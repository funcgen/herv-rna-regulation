#!/usr/bin/env python3
import argparse
import gzip
import os
import re
from typing import List, Tuple

import pandas as pd


def make_domain_hit_key_from_overlap_aname(a_name: str) -> str:
    if a_name is None or pd.isna(a_name):
        return "."
    s = str(a_name).strip()
    m = re.match(r'^(.+?_strand_[+-])_(\d+\|\d+-\d+)', s)
    if not m:
        return s
    return f"{m.group(1)}_{m.group(2)}"


def make_domain_hit_key_from_domains_row(query: str, locus: str) -> str:
    if (query is None or pd.isna(query)) or (locus is None or pd.isna(locus)):
        return "."
    q = str(query).strip()
    l = str(locus).strip()
    m = re.search(r'_(\d+\|\d+-\d+)$', q)
    if not m:
        return q
    return f"{l}_{m.group(1)}"


def parse_bname_tx_feature(bname: str):
    if bname is None or pd.isna(bname):
        return ".", ".", ".", ".", ".", "."
    s = str(bname).strip()
    parts = s.split("|")
    tx_id = parts[0] if len(parts) >= 1 else "."
    gene_id = parts[1] if len(parts) >= 2 else "."
    gene_name = parts[2] if len(parts) >= 3 else "."
    feature = parts[3] if len(parts) >= 4 else "."
    gene_type = parts[4] if len(parts) >= 5 else "."
    tx_type = parts[5] if len(parts) >= 6 else "."
    return tx_id, gene_id, gene_name, feature, gene_type, tx_type


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Return union of half-open intervals [start, end)."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:  # overlap/touch
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def union_overlap_bp_from_group(g: pd.DataFrame) -> int:
    """
    Compute unique overlapped bases for a group using ovl_start/ovl_end columns.
    This prevents overlap fractions > 1 caused by double-counting.
    """
    # Drop NA/invalid intervals
    sub = g[["ovl_start", "ovl_end"]].dropna()
    if sub.empty:
        return 0
    intervals = []
    for s, e in sub.itertuples(index=False, name=None):
        s = int(s); e = int(e)
        if e > s:
            intervals.append((s, e))
    if not intervals:
        return 0
    merged = _merge_intervals(intervals)
    return int(sum(e - s for s, e in merged))


def read_ov(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")

    expected = ["a_chrom","a_start","a_end","a_name","a_score","a_strand",
                "b_chrom","b_start","b_end","b_name","b_score","b_strand","overlap_bp"]
    if df.shape[1] != 13:
        raise ValueError(
            f"{path}: expected 13 columns, got {df.shape[1]}. "
            f"Likely malformed TSV (missing tabs)."
        )
    df.columns = expected

    for c in ["a_start","a_end","b_start","b_end","overlap_bp"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["domain_hit_key"] = df["a_name"].apply(make_domain_hit_key_from_overlap_aname)
    df["domain_len_bp"] = (df["a_end"] - df["a_start"]).astype("float")

    # Compute exact overlapped genomic interval
    df["ovl_start"] = df[["a_start", "b_start"]].max(axis=1)
    df["ovl_end"] = df[["a_end", "b_end"]].min(axis=1)
    df["ovl_len"] = (df["ovl_end"] - df["ovl_start"]).clip(lower=0)

    # This per-row frac is fine (but we won't SUM these)
    df["overlap_frac_of_domain"] = (df["ovl_len"] / df["domain_len_bp"]).fillna(0.0)

    parsed = df["b_name"].apply(parse_bname_tx_feature)
    df["tx_id"] = parsed.apply(lambda x: x[0])
    df["gene_id"] = parsed.apply(lambda x: x[1])
    df["gene_name"] = parsed.apply(lambda x: x[2])
    df["feature"] = parsed.apply(lambda x: x[3])
    df["gene_type"] = parsed.apply(lambda x: x[4])
    df["tx_type"] = parsed.apply(lambda x: x[5])

    return df


def load_tx_bounds_from_exons_bed_gz(exons_bed_gz: str):
    bounds = {}  # tx -> [min_start, max_end, strand]
    opener = gzip.open if exons_bed_gz.endswith(".gz") else open
    with opener(exons_bed_gz, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
            start = int(start); end = int(end)
            tx_id, *_ = name.split("|")
            if tx_id not in bounds:
                bounds[tx_id] = [start, end, strand]
            else:
                bounds[tx_id][0] = min(bounds[tx_id][0], start)
                bounds[tx_id][1] = max(bounds[tx_id][1], end)
    return bounds


def is_terminal_exon_overlap(row, tx_bounds):
    tx = row["tx_id"]
    if tx not in tx_bounds:
        return False
    tx_min, tx_max, tx_strand = tx_bounds[tx]
    if tx_strand == "+":
        return int(row["b_end"]) == int(tx_max)
    elif tx_strand == "-":
        return int(row["b_start"]) == int(tx_min)
    return False


def main():
    ap = argparse.ArgumentParser(
        description="Transcript-aware domain overlaps vs exon/CDS/UTR (5'/3'), keeping ALL matches."
    )
    ap.add_argument("--domains-tsv", required=True, help="ERV_GyDB_v6_filtered.tsv")
    ap.add_argument("--ov-exon", required=True, help="domains_vs_gencode_tx_exons.same_strand.wo.tsv")
    ap.add_argument("--ov-cds", required=True, help="domains_vs_gencode_cds.same_strand.wo.tsv")
    ap.add_argument("--ov-utr5p3p", required=True, help="domains_vs_gencode_utr_5p3p.same_strand.wo.tsv")
    ap.add_argument("--tx-exons-bed", required=True, help="bed/gencode.v48.exons.sorted.bed.gz (full exons bed)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-domain-overlap", type=float, default=0.80,
                    help="Min fraction of DOMAIN length overlapped by the feature(s) within a transcript.")
    ap.add_argument("--min-coverage", type=float, default=0.50,
                    help="Min domain Coverage from domains.tsv.")
    ap.add_argument("--write-wide", action="store_true",
                    help="Also write a wide table per (domain_hit_key, tx_id) with feature columns.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    dom = pd.read_csv(args.domains_tsv, sep="\t", dtype=str)
    required_dom = ["Query", "DomainType", "Coverage", "Gene", "Domain", "Subfamily", "Locus"]
    missing = [c for c in required_dom if c not in dom.columns]
    if missing:
        raise ValueError(f"Domains TSV missing columns: {missing}")

    dom["domain_hit_key"] = dom.apply(
        lambda r: make_domain_hit_key_from_domains_row(r["Query"], r["Locus"]), axis=1
    )
    dom["Coverage_num"] = pd.to_numeric(dom["Coverage"], errors="coerce")

    ov_exon = read_ov(args.ov_exon)
    ov_cds = read_ov(args.ov_cds)
    ov_utr = read_ov(args.ov_utr5p3p)

    ov_exon["feature_class"] = "exon"
    ov_cds["feature_class"] = "CDS"
    ov_utr["feature_class"] = ov_utr["feature"].fillna("UTR")

    # terminal exon flag (computed from FULL exons bed)
    tx_bounds = load_tx_bounds_from_exons_bed_gz(args.tx_exons_bed)
    ov_exon["is_terminal_exon"] = ov_exon.apply(lambda r: is_terminal_exon_overlap(r, tx_bounds), axis=1)
    ov_cds["is_terminal_exon"] = False
    ov_utr["is_terminal_exon"] = False

    ov_all = pd.concat([ov_exon, ov_cds, ov_utr], ignore_index=True)

    # --- collapse within (domain_hit_key, tx_id, feature_class) using UNION overlap bp
    agg = (ov_all
           .groupby(["domain_hit_key", "tx_id", "feature_class"], as_index=False)
           .agg(
               chrom=("a_chrom", "first"),
               start=("a_start", "min"),
               end=("a_end", "max"),
               strand=("a_strand", "first"),
               domain_len_bp=("domain_len_bp", "first"),
               overlap_bp_total=("ovl_len", "sum"),  # placeholder, overwritten below
               gene_id=("gene_id", "first"),
               gene_name=("gene_name", "first"),
               gene_type=("gene_type", "first"),
               tx_type=("tx_type", "first"),
               any_terminal_exon_overlap=("is_terminal_exon", "max"),
           )
          )

    # Compute true unique overlap bp per group (prevents >1 fractions)
    # We need access to per-row ovl_start/ovl_end, so compute separately and merge.
    # Compute union overlap bp per group (version-safe)
    union_bp = (
        ov_all
        .groupby(["domain_hit_key", "tx_id", "feature_class"])
        .apply(union_overlap_bp_from_group)
    )

    # union_bp is a Series with a MultiIndex -> convert to DataFrame safely
    union_bp = union_bp.reset_index()
    union_bp = union_bp.rename(columns={0: "overlap_bp_union"})
    agg = agg.merge(union_bp, on=["domain_hit_key", "tx_id", "feature_class"], how="left")
    agg["overlap_bp_total"] = pd.to_numeric(agg["overlap_bp_union"], errors="coerce").fillna(0).astype(int)
    agg.drop(columns=["overlap_bp_union"], inplace=True)

    agg["overlap_frac_of_domain"] = (agg["overlap_bp_total"] / agg["domain_len_bp"]).fillna(0.0)
    # Clamp as a final safety net (should be unnecessary after union)
    agg["overlap_frac_of_domain"] = agg["overlap_frac_of_domain"].clip(lower=0.0, upper=1.0)
    # Human-readable genomic coordinates for the domain hit
    agg["domain_coord"] = (
        agg["chrom"].astype(str)
        + ":"
        + agg["start"].astype(int).astype(str)
        + "-"
        + agg["end"].astype(int).astype(str)
    )

    merged = agg.merge(dom, how="left", on="domain_hit_key")
    n_match = merged["Coverage_num"].notna().sum()
    print(f"[DEBUG] domain key matches with domains.tsv: {n_match}/{merged.shape[0]}", flush=True)

    kept = merged[
        (merged["overlap_frac_of_domain"] >= args.min_domain_overlap) &
        (merged["Coverage_num"].fillna(0) >= args.min_coverage)
    ].copy()

    if kept.empty:
        print("[WARN] No rows passed filters. Try lowering --min-domain-overlap or --min-coverage.")
        return

    out_long = os.path.join(args.outdir, "domain_tx_feature_overlaps.long.tsv")
    out_cols = [
        "domain_hit_key",
        "domain_coord",
        "tx_id", "gene_id", "gene_name", "gene_type", "tx_type",
        "feature_class",
        "any_terminal_exon_overlap",
        "chrom", "start", "end", "strand",
        "domain_len_bp", "overlap_bp_total", "overlap_frac_of_domain",
        "DomainType", "Coverage", "Coverage_num",
        "Gene", "Domain", "Subfamily", "Locus"
    ]
    out_cols = [c for c in out_cols if c in kept.columns]
    kept[out_cols].sort_values(
        ["DomainType", "feature_class", "overlap_frac_of_domain", "Coverage_num"],
        ascending=[True, True, False, False]
    ).to_csv(out_long, sep="\t", index=False)

    if args.write_wide:
        wide = kept.pivot_table(
            index=[
                "domain_hit_key", "domain_coord",
                "tx_id", "gene_id", "gene_name", "gene_type", "tx_type",
                "DomainType", "Coverage", "Coverage_num",
                "Gene", "Domain", "Subfamily", "Locus"
            ],
            columns="feature_class",
            values="overlap_frac_of_domain",
            aggfunc="max",
            fill_value=0.0
        ).reset_index()
        term = (kept.groupby(["domain_hit_key", "tx_id"], as_index=False)["any_terminal_exon_overlap"].max())
        wide = wide.merge(term, on=["domain_hit_key", "tx_id"], how="left")
        out_wide = os.path.join(args.outdir, "domain_tx_feature_overlaps.wide.tsv")
        wide.to_csv(out_wide, sep="\t", index=False)

    print(f"[OK] Wrote: {out_long}")
    if args.write_wide:
        print(f"[OK] Wrote wide table too.")
    print(f"[OK] Kept {kept.shape[0]} (domain_hit_key, tx_id, feature_class) rows after filters.")


if __name__ == "__main__":
    main()
