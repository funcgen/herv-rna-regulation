#!/usr/bin/env python3
"""
03_transcript_context step 1 — extract BED assets from the GENCODE GTF.

Refactor of scripts/context/01_build_gencode_feature_beds.py. Every schema and
semantic rule is reproduced exactly; only paths, argument handling, validation
and memory behaviour changed.

OUTPUTS — five gzipped BED6 files:
    genes/<prefix>.genes.bed.gz
    transcripts/<prefix>.transcripts.bed.gz
    features/<prefix>.exons.bed.gz
    features/<prefix>.cds.bed.gz
    features/<prefix>.utr.bed.gz

SEMANTICS PRESERVED FROM THE HISTORICAL PRODUCER
------------------------------------------------
* COORDINATES. GTF is 1-based inclusive, BED 0-based half-open:
      start0 = start1 - 1   (clamped at 0)
      end    = end1
* BED6 columns: chrom, start0, end1, name, score (always 0), strand.
* NAME FIELDS, pipe-delimited, per feature class:
      gene        gene_id|gene_name|gene_type
      transcript  transcript_id|gene_id|gene_name|transcript_type
      exon        transcript_id|gene_id|gene_name|exon|gene_type|transcript_type
      CDS         transcript_id|gene_id|gene_name|CDS|gene_type|transcript_type
      UTR         transcript_id|gene_id|gene_name|<FEATURE>|gene_type|transcript_type
  For UTRs the LITERAL GTF feature name is kept, so `UTR`, `five_prime_UTR` and
  `three_prime_UTR` stay distinguishable in one combined file. The 5'/3' split
  is step 4's job, not this one.
* IDENTIFIERS ARE VERBATIM. Ensembl version suffixes are NOT stripped.
  gene_name defaults to "NA"; gene_type falls back to gene_biotype then "NA";
  transcript_type falls back to transcript_biotype then "NA".
* NO CHROMOSOME FILTERING. Every sequence in the GTF is kept.
* NO SORTING. Records are emitted in GTF order; genome-order sorting is step 3.
* Comment lines are skipped. Lines without 9 fields, or with non-integer
  coordinates, or missing the identifiers their class requires, are counted as
  skipped rather than aborting.

CHANGED — memory only, not content
----------------------------------
The historical version accumulated all five row lists in memory and wrote at the
end; on this 1.8 GB GTF that is several GB of Python objects. This version
STREAMS: every record is written as it is read, so memory is constant. Row order
and bytes are unaffected. The tqdm progress bar and its preliminary
line-counting pass are dropped with it -- they required a second full read of the
GTF purely for a progress display.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import sys
from pathlib import Path

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')
UTR_FEATURES = ("UTR", "five_prime_UTR", "three_prime_UTR")


def die(msg: str) -> "None":
    sys.stderr.write(f"[error] {msg}\n")
    sys.exit(1)


def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_attrs(attr_field: str) -> "dict[str, str]":
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(attr_field)}


def bed_line(chrom: str, start1: int, end1: int, name: str, strand: str) -> str:
    start0 = start1 - 1
    if start0 < 0:
        start0 = 0
    return f"{chrom}\t{start0}\t{end1}\t{name}\t0\t{strand}\n"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (plain or .gz)")
    ap.add_argument("--outdir", required=True, help="output base directory")
    ap.add_argument("--prefix", default="gencode.v48")
    ap.add_argument("--expect-gtf-sha256", default=None,
                    help="refuse to run unless the GTF matches this checksum")
    args = ap.parse_args()

    gtf = Path(args.gtf)
    if not gtf.is_file():
        die(f"GTF not found: {gtf}")

    if args.expect_gtf_sha256:
        got = sha256_file(gtf)
        if got != args.expect_gtf_sha256:
            die(f"GTF checksum mismatch\n       expected {args.expect_gtf_sha256}\n"
                f"       observed {got}\n     Refusing to build from an unverified annotation.")
        print(f"[ok] GTF sha256 verified: {got}")

    outdir = Path(args.outdir)
    prefix = args.prefix
    targets = {
        "genes":       outdir / "genes" / f"{prefix}.genes.bed.gz",
        "transcripts": outdir / "transcripts" / f"{prefix}.transcripts.bed.gz",
        "exons":       outdir / "features" / f"{prefix}.exons.bed.gz",
        "cds":         outdir / "features" / f"{prefix}.cds.bed.gz",
        "utr":         outdir / "features" / f"{prefix}.utr.bed.gz",
    }
    for p in targets.values():
        p.parent.mkdir(parents=True, exist_ok=True)

    counts = {k: 0 for k in targets}
    skipped = 0
    chroms: "set[str]" = set()
    utr_kinds: "dict[str, int]" = {}

    # Stage to temporary files and publish only on success.
    tmp = {k: v.with_suffix(v.suffix + ".tmp") for k, v in targets.items()}
    handles = {k: gzip.open(p, "wt") for k, p in tmp.items()}

    try:
        with open_maybe_gzip(str(gtf), "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    skipped += 1
                    continue

                chrom, _source, feature, start, end, _score, strand, _frame, attrs = parts
                if feature not in ("gene", "transcript", "exon", "CDS") \
                        and feature not in UTR_FEATURES:
                    continue                      # not an emitted class

                try:
                    start1 = int(start)
                    end1 = int(end)
                except ValueError:
                    skipped += 1
                    continue

                a = parse_attrs(attrs)
                gene_id = a.get("gene_id")
                gene_name = a.get("gene_name", "NA")
                gene_type = a.get("gene_type", a.get("gene_biotype", "NA"))
                transcript_id = a.get("transcript_id")
                transcript_type = a.get("transcript_type",
                                        a.get("transcript_biotype", "NA"))

                if feature == "gene":
                    if not gene_id:
                        skipped += 1
                        continue
                    name = f"{gene_id}|{gene_name}|{gene_type}"
                    handles["genes"].write(bed_line(chrom, start1, end1, name, strand))
                    counts["genes"] += 1

                elif feature == "transcript":
                    if not (gene_id and transcript_id):
                        skipped += 1
                        continue
                    name = f"{transcript_id}|{gene_id}|{gene_name}|{transcript_type}"
                    handles["transcripts"].write(bed_line(chrom, start1, end1, name, strand))
                    counts["transcripts"] += 1

                elif feature == "exon":
                    if not (gene_id and transcript_id):
                        skipped += 1
                        continue
                    name = (f"{transcript_id}|{gene_id}|{gene_name}|exon|"
                            f"{gene_type}|{transcript_type}")
                    handles["exons"].write(bed_line(chrom, start1, end1, name, strand))
                    counts["exons"] += 1

                elif feature == "CDS":
                    if not (gene_id and transcript_id):
                        skipped += 1
                        continue
                    name = (f"{transcript_id}|{gene_id}|{gene_name}|CDS|"
                            f"{gene_type}|{transcript_type}")
                    handles["cds"].write(bed_line(chrom, start1, end1, name, strand))
                    counts["cds"] += 1

                else:                              # a UTR class
                    if not (gene_id and transcript_id):
                        skipped += 1
                        continue
                    name = (f"{transcript_id}|{gene_id}|{gene_name}|{feature}|"
                            f"{gene_type}|{transcript_type}")
                    handles["utr"].write(bed_line(chrom, start1, end1, name, strand))
                    counts["utr"] += 1
                    utr_kinds[feature] = utr_kinds.get(feature, 0) + 1

                chroms.add(chrom)
    except Exception:
        for h in handles.values():
            h.close()
        for p in tmp.values():
            if p.exists():
                p.unlink()
        raise

    for h in handles.values():
        h.close()

    if any(v == 0 for v in counts.values()):
        for p in tmp.values():
            if p.exists():
                p.unlink()
        die(f"a feature class produced zero records: {counts}. Refusing to publish.")

    for k, p in tmp.items():
        p.chmod(0o644)
        p.replace(targets[k])

    print(f"[ok] genes       {counts['genes']:>9,}  {targets['genes']}")
    print(f"[ok] transcripts {counts['transcripts']:>9,}  {targets['transcripts']}")
    print(f"[ok] exons       {counts['exons']:>9,}  {targets['exons']}")
    print(f"[ok] CDS         {counts['cds']:>9,}  {targets['cds']}")
    print(f"[ok] UTR         {counts['utr']:>9,}  {targets['utr']}")
    print(f"[info] UTR feature kinds: {dict(sorted(utr_kinds.items()))}")
    print(f"[info] sequences: {len(chroms)}")
    print(f"[info] skipped lines: {skipped:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
