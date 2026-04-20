#!/usr/bin/env python3
"""
Build BED assets from a GENCODE GTF with a progress bar.

Outputs (BED6, gzipped):
- context/assets/genes/gencode.v48.genes.bed.gz
- context/assets/transcripts/gencode.v48.transcripts.bed.gz
- context/assets/features/gencode.v48.exons.bed.gz
- context/assets/features/gencode.v48.cds.bed.gz
- context/assets/features/gencode.v48.utr.bed.gz

Notes:
- GTF is 1-based inclusive; BED is 0-based half-open:
  start0 = start1 - 1; end = end1
- Uses a 2-pass scan to obtain total line count for tqdm.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Iterator, TextIO, Tuple

try:
    from tqdm import tqdm
except ImportError:
    sys.stderr.write(
        "[ERROR] tqdm is not installed in this environment.\n"
        "Install it with: pip install tqdm\n"
        "or add it to your conda env.\n"
    )
    raise

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def open_maybe_gzip(path: str, mode: str = "rt") -> TextIO:
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_attrs(attr_field: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(attr_field)}


def bed_fields(chrom: str, start1: int, end1: int, name: str, strand: str) -> Tuple[str, int, int, str, int, str]:
    start0 = start1 - 1
    if start0 < 0:
        start0 = 0
    return (chrom, start0, end1, name, 0, strand)


def iter_gtf_lines(gtf_path: str) -> Iterator[str]:
    with open_maybe_gzip(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            yield line.rstrip("\n")


def count_noncomment_lines(gtf_path: str) -> int:
    n = 0
    with open_maybe_gzip(gtf_path, "rt") as fh:
        for line in fh:
            if line and not line.startswith("#"):
                n += 1
    return n


def write_bed_gz(out_path: Path, rows) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as out:
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (plain or .gz)")
    ap.add_argument("--outdir", required=True, help="Output base dir (e.g. context/assets)")
    ap.add_argument("--prefix", default="gencode.v48", help="Output filename prefix (default: gencode.v48)")
    ap.add_argument("--log", default=None, help="Optional log file")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    prefix = args.prefix

    genes_out = outdir / "genes" / f"{prefix}.genes.bed.gz"
    tx_out    = outdir / "transcripts" / f"{prefix}.transcripts.bed.gz"
    exon_out  = outdir / "features" / f"{prefix}.exons.bed.gz"
    cds_out   = outdir / "features" / f"{prefix}.cds.bed.gz"
    utr_out   = outdir / "features" / f"{prefix}.utr.bed.gz"

    def log(msg: str):
        if args.log:
            Path(args.log).parent.mkdir(parents=True, exist_ok=True)
            with open(args.log, "a") as lf:
                lf.write(msg + "\n")
        else:
            print(msg, file=sys.stderr)

    log(f"[INFO] Counting lines in {args.gtf} (for progress bar)...")
    total = count_noncomment_lines(args.gtf)
    log(f"[INFO] Non-comment lines: {total:,}")

    genes_rows = []
    tx_rows = []
    exon_rows = []
    cds_rows = []
    utr_rows = []

    n_gene = n_tx = n_exon = n_cds = n_utr = 0
    n_skipped = 0

    pbar = tqdm(total=total, unit="lines", desc="Parsing GTF", smoothing=0.05)
    try:
        for line in iter_gtf_lines(args.gtf):
            pbar.update(1)

            parts = line.split("\t")
            if len(parts) != 9:
                n_skipped += 1
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            try:
                start1 = int(start)
                end1 = int(end)
            except ValueError:
                n_skipped += 1
                continue

            a = parse_attrs(attrs)

            gene_id = a.get("gene_id")
            gene_name = a.get("gene_name", "NA")
            gene_type = a.get("gene_type", a.get("gene_biotype", "NA"))

            transcript_id = a.get("transcript_id")
            transcript_type = a.get("transcript_type", a.get("transcript_biotype", "NA"))

            if feature == "gene":
                if not gene_id:
                    n_skipped += 1
                    continue
                name = f"{gene_id}|{gene_name}|{gene_type}"
                genes_rows.append(bed_fields(chrom, start1, end1, name, strand))
                n_gene += 1

            elif feature == "transcript":
                if not (gene_id and transcript_id):
                    n_skipped += 1
                    continue
                name = f"{transcript_id}|{gene_id}|{gene_name}|{transcript_type}"
                tx_rows.append(bed_fields(chrom, start1, end1, name, strand))
                n_tx += 1

            elif feature == "exon":
                if not (gene_id and transcript_id):
                    n_skipped += 1
                    continue
                name = f"{transcript_id}|{gene_id}|{gene_name}|exon|{gene_type}|{transcript_type}"
                exon_rows.append(bed_fields(chrom, start1, end1, name, strand))
                n_exon += 1

            elif feature == "CDS":
                if not (gene_id and transcript_id):
                    n_skipped += 1
                    continue
                name = f"{transcript_id}|{gene_id}|{gene_name}|CDS|{gene_type}|{transcript_type}"
                cds_rows.append(bed_fields(chrom, start1, end1, name, strand))
                n_cds += 1

            elif feature in ("UTR", "five_prime_UTR", "three_prime_UTR"):
                if not (gene_id and transcript_id):
                    n_skipped += 1
                    continue
                name = f"{transcript_id}|{gene_id}|{gene_name}|{feature}|{gene_type}|{transcript_type}"
                utr_rows.append(bed_fields(chrom, start1, end1, name, strand))
                n_utr += 1

    finally:
        pbar.close()

    # Write outputs
    log("[INFO] Writing BED assets...")
    write_bed_gz(genes_out, genes_rows)
    write_bed_gz(tx_out, tx_rows)
    write_bed_gz(exon_out, exon_rows)
    write_bed_gz(cds_out, cds_rows)
    write_bed_gz(utr_out, utr_rows)

    log(f"[OK] Wrote genes:        {genes_out} ({n_gene:,})")
    log(f"[OK] Wrote transcripts:  {tx_out} ({n_tx:,})")
    log(f"[OK] Wrote exons:        {exon_out} ({n_exon:,})")
    log(f"[OK] Wrote CDS:          {cds_out} ({n_cds:,})")
    log(f"[OK] Wrote UTR:          {utr_out} ({n_utr:,})")
    log(f"[INFO] Skipped lines:    {n_skipped:,}")


if __name__ == "__main__":
    main()
