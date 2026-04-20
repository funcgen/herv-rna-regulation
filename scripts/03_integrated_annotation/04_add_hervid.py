#!/usr/bin/env python3
"""
Add HERV_id (from GTF gene_id) into internal_integrated.tsv by coordinate matching.

Assumptions:
- internal_integrated.tsv uses 0-based start and 1-based end (BED-like half-open): [start, end)
- GTF uses 1-based inclusive coordinates: [start, end]
- Therefore: tsv_start + 1 == gtf_start AND tsv_end == gtf_end
- Match key: chrom, start, end, strand (and geneRegion=="internal" in GTF)

Usage:
  ./add_herv_id_from_gtf.py \
    --tsv internal_integrated.tsv \
    --gtf hervs_full_ltr_int_v2.gtf \
    --out internal_integrated.with_herv_id.tsv

Optional:
  --tsv-start-col start --tsv-end-col end --tsv-strand-col strand --tsv-chrom-col chrom
"""

import argparse
import gzip
import re
import sys
import pandas as pd


def open_maybe_gz(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


_attr_re = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_gtf_attributes(attr: str) -> dict:
    # robust enough for typical GTF attributes: key "value";
    # returns dict like {"gene_id": "...", "transcript_id": "...", ...}
    d = {}
    for m in _attr_re.finditer(attr):
        d[m.group(1)] = m.group(2)
    return d


def read_internal_gtf(gtf_path: str) -> pd.DataFrame:
    rows = []
    with open_maybe_gz(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            ad = parse_gtf_attributes(attrs)

            # keep only geneRegion == internal
            if ad.get("geneRegion") != "internal":
                continue

            gene_id = ad.get("gene_id")
            if gene_id is None:
                continue

            # GTF start/end are 1-based inclusive in file
            start_i = int(start)
            end_i = int(end)

            rows.append((chrom, start_i, end_i, strand, gene_id))

    gtf_df = pd.DataFrame(rows, columns=["chrom", "gtf_start_1based", "gtf_end_1based", "strand", "HERV_id"])

    # Check duplicates (same coords mapping to multiple gene_ids) – should not happen, but guard anyway
    dup = gtf_df.duplicated(subset=["chrom", "gtf_start_1based", "gtf_end_1based", "strand"], keep=False)
    if dup.any():
        bad = gtf_df.loc[dup].sort_values(["chrom", "gtf_start_1based", "gtf_end_1based", "strand"])
        raise SystemExit(
            "ERROR: duplicate coordinate keys in GTF (same chrom/start/end/strand but different gene_id).\n"
            f"First few duplicates:\n{bad.head(20).to_string(index=False)}"
        )

    return gtf_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, help="internal_integrated.tsv (tab-separated with header)")
    ap.add_argument("--gtf", required=True, help="hervs_full_ltr_int_v2.gtf (can be .gz)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--tsv-chrom-col", default="chrom")
    ap.add_argument("--tsv-start-col", default="start")
    ap.add_argument("--tsv-end-col", default="end")
    ap.add_argument("--tsv-strand-col", default="strand")
    ap.add_argument("--keep-gtf-cols", action="store_true",
                    help="If set, keep gtf_start_1based/gtf_end_1based after merge (debug).")
    args = ap.parse_args()

    # Read TSV
    tsv = pd.read_csv(args.tsv, sep="\t", dtype=str)
    for c in (args.tsv_chrom_col, args.tsv_start_col, args.tsv_end_col, args.tsv_strand_col):
        if c not in tsv.columns:
            raise SystemExit(f"ERROR: TSV is missing required column: {c}")

    # Coerce coords
    tsv[args.tsv_start_col] = tsv[args.tsv_start_col].astype(int)
    tsv[args.tsv_end_col] = tsv[args.tsv_end_col].astype(int)

    # Read internal-only GTF mapping
    gtf_df = read_internal_gtf(args.gtf)

    # Build merge keys:
    # TSV: start is 0-based; convert to 1-based for matching GTF start
    tsv = tsv.copy()
    tsv["__start_1based"] = tsv[args.tsv_start_col] + 1
    tsv["__end_1based"] = tsv[args.tsv_end_col]

    # Merge
    merged = tsv.merge(
        gtf_df,
        how="left",
        left_on=[args.tsv_chrom_col, "__start_1based", "__end_1based", args.tsv_strand_col],
        right_on=["chrom", "gtf_start_1based", "gtf_end_1based", "strand"],
        suffixes=("", "_gtf"),
    )

    # Summary
    n_total = len(merged)
    n_hit = merged["HERV_id"].notna().sum()
    n_miss = n_total - n_hit
    sys.stderr.write(f"[INFO] Total TSV rows: {n_total}\n")
    sys.stderr.write(f"[INFO] Mapped to GTF internal coords: {n_hit} ({(n_hit/n_total)*100:.2f}%)\n")
    sys.stderr.write(f"[INFO] Unmapped: {n_miss} ({(n_miss/n_total)*100:.2f}%)\n")

    if n_miss > 0:
        # show a few examples
        ex = merged.loc[merged["HERV_id"].isna(), [args.tsv_chrom_col, args.tsv_start_col, args.tsv_end_col, args.tsv_strand_col]].head(15)
        sys.stderr.write("[WARN] Example unmapped rows (chrom/start/end/strand):\n")
        sys.stderr.write(ex.to_string(index=False) + "\n")

    # Clean up columns
    drop_cols = ["__start_1based", "__end_1based", "chrom_gtf", "strand_gtf"]
    if not args.keep_gtf_cols:
        drop_cols += ["gtf_start_1based", "gtf_end_1based"]
    # Only drop those that exist
    drop_cols = [c for c in drop_cols if c in merged.columns]
    merged.drop(columns=drop_cols, inplace=True)

    # Put HERV_id as the first column (if present)
    if "HERV_id" in merged.columns:
        cols = ["HERV_id"] + [c for c in merged.columns if c != "HERV_id"]
        merged = merged.loc[:, cols]

    # Write output
    merged.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[INFO] Wrote: {args.out}\n")


if __name__ == "__main__":
    main()
