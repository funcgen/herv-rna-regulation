#!/usr/bin/env python3
"""
Reconstruct genomic coordinates for primary FIMO motif calls.

FIMO was run with --no-pgc, so it reports positions relative to each scanned
sequence and does not attempt to decode the sequence name. This parser does that
decoding and nothing else: no filtering, no deduplication, no recomputation of
scores, p-values or q-values. Every raw call becomes exactly one parsed record.

Sequence-name grammar
---------------------
    <subfamily>-int_pos_<chrom>_<internal_start>_<internal_end>_strand_<+|->[...]

matched with a regex search, so trailing annotation such as
"_divergence_6.3_1" is tolerated. internal_start / internal_end are BED-style
(0-based, end-exclusive), as produced by bedtools getfasta.

Coordinate reconstruction
-------------------------
FIMO start/stop are 1-based inclusive relative to the scanned sequence. The
sense FASTA is transcript-oriented, so a locus on the genomic minus strand was
scanned as the reverse complement of the genome and its coordinates must be
mirrored:

    internal_strand '+'   (forward mapping)
        genomic_start = internal_start + (start - 1)
        genomic_end   = internal_start + stop

    internal_strand '-'   (reverse mapping)
        genomic_start = internal_end - stop
        genomic_end   = internal_end - (start - 1)

Both yield BED coordinates (0-based, end-exclusive) whose length equals the
motif span, stop - start + 1.

Motif strand
------------
FIMO ran with --norc, so it reports '+' for every call: the strand relative to
the scanned sequence. For a minus-strand locus the scanned sequence is the
reverse complement of the genome, so the genomic strand of the match is flipped.
match_strand therefore equals internal_strand for this sense run; it is computed
from the FIMO strand rather than copied, to preserve the historical logic for an
antisense run.

Output ordering
---------------
Sorted by (chrom, genomic_start, genomic_end) with a stable sort, so ties keep
their original FIMO order. The BED is written from that same ordering.

Behaviour is identical to the historical implementation. Only paths, error
handling, determinism guarantees and validation have changed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import sys
import tempfile
from pathlib import Path

SEQUENCE_NAME_RE = re.compile(
    r"_(?P<chrom>(chr)?[^_]+)_(?P<start>\d+)_(?P<end>\d+)_strand_(?P<strand>[+-])"
)

RAW_COLUMNS = [
    "motif_id", "motif_alt_id", "sequence_name",
    "start", "stop", "strand",
    "score", "p-value", "q-value", "matched_sequence",
]
DERIVED_COLUMNS = [
    "chrom", "internal_start", "internal_end", "internal_strand", "match_strand",
]
TSV_COLUMNS = RAW_COLUMNS + DERIVED_COLUMNS


def die(message: str) -> None:
    print(f"[ERROR] {message}", file=sys.stderr)
    sys.exit(1)


def flip_strand(strand: str) -> str:
    return {"+": "-", "-": "+"}.get(strand, strand)


def decode_sequence_name(name: str):
    """Return (chrom, internal_start, internal_end, internal_strand)."""
    match = SEQUENCE_NAME_RE.search(name)
    if not match:
        raise ValueError(f"sequence name does not match the expected grammar: {name!r}")
    return (
        match.group("chrom"),
        int(match.group("start")),
        int(match.group("end")),
        match.group("strand"),
    )


def read_fimo(path: Path):
    """Yield raw FIMO data rows as lists of strings, preserving every value verbatim.

    Header, blank lines and '#' comment lines (the FIMO trailer) are excluded.
    Raw values are never parsed into floats, so score, p-value and q-value reach
    the output exactly as FIMO wrote them.
    """
    if not path.is_file():
        die(f"input FIMO TSV not found: {path}")

    with path.open(encoding="utf-8") as handle:
        header_line = handle.readline().rstrip("\n")
        header = header_line.split("\t")
        if header != RAW_COLUMNS:
            die(
                "unexpected FIMO header.\n"
                f"  expected: {RAW_COLUMNS}\n"
                f"  observed: {header}"
            )
        for line_number, raw_line in enumerate(handle, start=2):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != len(RAW_COLUMNS):
                die(
                    f"{path}:{line_number}: expected {len(RAW_COLUMNS)} columns, "
                    f"found {len(fields)}"
                )
            yield line_number, fields


def scientific_row_sha256(path: Path) -> str:
    """Hash of the raw FIMO data rows, matching the Stage 3 regression definition."""
    digest = hashlib.sha256()
    for _, fields in read_fimo(path):
        digest.update(("\t".join(fields) + "\n").encode("utf-8"))
    return digest.hexdigest()


def parse_records(path: Path):
    """Decode every raw call into a parsed record. Fails loudly; skips nothing."""
    records = []
    for line_number, fields in read_fimo(path):
        row = dict(zip(RAW_COLUMNS, fields))
        try:
            chrom, internal_start, internal_end, internal_strand = decode_sequence_name(
                row["sequence_name"]
            )
        except ValueError as exc:
            die(f"{path}:{line_number}: {exc}")

        try:
            fimo_start = int(row["start"])
            fimo_stop = int(row["stop"])
        except ValueError:
            die(f"{path}:{line_number}: non-integer FIMO start/stop")

        if internal_strand == "+":
            genomic_start = internal_start + (fimo_start - 1)
            genomic_end = internal_start + fimo_stop
            match_strand = row["strand"]
        else:
            genomic_start = internal_end - fimo_stop
            genomic_end = internal_end - (fimo_start - 1)
            match_strand = flip_strand(row["strand"])

        if genomic_end <= genomic_start:
            die(
                f"{path}:{line_number}: genomic_end <= genomic_start "
                f"({genomic_end} <= {genomic_start}); coordinate convention broken"
            )
        if genomic_start < 0:
            die(f"{path}:{line_number}: negative genomic_start {genomic_start}")

        records.append(
            {
                **row,
                "chrom": chrom,
                "internal_start": internal_start,
                "internal_end": internal_end,
                "internal_strand": internal_strand,
                "match_strand": match_strand,
                "_genomic_start": genomic_start,
                "_genomic_end": genomic_end,
            }
        )
    return records


def write_atomic(path: Path, write_fn) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(path.parent), prefix=".tmp_", suffix=".part")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with tmp.open("w", encoding="utf-8", newline="") as handle:
            write_fn(handle)
        # mkstemp creates 0600; release artifacts must be world-readable.
        os.chmod(tmp, 0o644)
        tmp.replace(path)
    finally:
        if tmp.exists():
            tmp.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reconstruct genomic coordinates for primary FIMO motif calls."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--out-bed", required=True, type=Path)
    parser.add_argument(
        "--expect-input-sha256",
        default=None,
        help="Scientific-row checksum the input must have (Stage 3 regression reference).",
    )
    parser.add_argument("--expect-rows", type=int, default=None)
    args = parser.parse_args()

    if args.expect_input_sha256:
        observed = scientific_row_sha256(args.input)
        if observed != args.expect_input_sha256:
            die(
                "input scientific-row checksum mismatch; refusing to parse\n"
                f"  expected {args.expect_input_sha256}\n"
                f"  observed {observed}"
            )
        print(f"[info] input scientific-row sha256 verified: {observed}")

    records = parse_records(args.input)
    n = len(records)
    print(f"[info] parsed records        : {n}")

    if args.expect_rows is not None and n != args.expect_rows:
        die(f"parsed {n} records, expected {args.expect_rows}")

    # Stable sort: ties preserve original FIMO order.
    records.sort(key=lambda r: (r["chrom"], r["_genomic_start"], r["_genomic_end"]))

    n_loci = len({r["sequence_name"] for r in records})
    print(f"[info] unique motif-positive loci: {n_loci}")

    def write_tsv(handle):
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for r in records:
            writer.writerow([r[c] for c in TSV_COLUMNS])

    def write_bed(handle):
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for r in records:
            bed_name = f"{r['sequence_name']}_{r['motif_id']}_{r['motif_alt_id']}"
            writer.writerow(
                [
                    r["chrom"],
                    r["_genomic_start"],
                    r["_genomic_end"],
                    bed_name,
                    r["score"],
                    r["match_strand"],
                ]
            )

    write_atomic(args.out_tsv, write_tsv)
    write_atomic(args.out_bed, write_bed)

    print(f"[ok] wrote {args.out_tsv}")
    print(f"[ok] wrote {args.out_bed}")


if __name__ == "__main__":
    main()
