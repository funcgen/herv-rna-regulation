#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
import pandas as pd
from typing import Dict, List, Optional


def read_bed(path: str) -> pd.DataFrame:
    """
    Reads BED6. Supports .gz.
    Columns: chrom, start, end, name, score, strand
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        df = pd.read_csv(
            f, sep="\t", header=None,
            names=["chrom", "start", "end", "locid", "score", "strand"],
            dtype={"chrom": str, "start": int, "end": int, "locid": str, "score": str, "strand": str},
        )
    return df


def safe_read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def parse_bname_gene_fields(bname: str) -> Dict[str, str]:
    """
    Parse bedtools 'b_name' field to extract gene_id, gene_symbol, gene_type.
    Handles formats you showed, e.g.:
      - ENSG...|MGMT|protein_coding
      - ENSG...|ENSG...|lncRNA
      - ENST...|ENSG...|TAF12|CDS|protein_coding|protein_coding
      - ENSG...|ENSG...|lncRNA|intron
    Returns dict with keys: gene_id, gene_symbol, gene_type
    """
    if bname is None or pd.isna(bname):
        return {"gene_id": ".", "gene_symbol": ".", "gene_type": "."}

    s = str(bname).strip()
    if s in ["", ".", "NA", "NaN"]:
        return {"gene_id": ".", "gene_symbol": ".", "gene_type": "."}

    parts = [p.strip() for p in s.split("|") if p.strip()]
    gene_id = "."
    gene_symbol = "."
    gene_type = "."

    # gene_id: first ENSG token anywhere
    for p in parts:
        if p.startswith("ENSG"):
            gene_id = p
            break

    # gene_symbol
    if len(parts) >= 2 and not parts[1].startswith("ENSG") and not parts[1].startswith("ENST"):
        gene_symbol = parts[1]
    else:
        gene_symbol = gene_id if gene_id != "." else "."

    common_biotypes = {
        "protein_coding", "lncRNA", "pseudogene",
        "processed_pseudogene", "unprocessed_pseudogene",
        "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
        "unitary_pseudogene", "miRNA", "snRNA", "snoRNA", "rRNA",
        "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
        "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
        "TEC", "misc_RNA"
    }

    for p in parts:
        if p in common_biotypes:
            gene_type = p
            break

    if gene_type == "." and len(parts) >= 1:
        gene_type = parts[-1]

    return {"gene_id": gene_id, "gene_symbol": gene_symbol, "gene_type": gene_type}


def parse_bname_list_gene_types(bname_list: str) -> str:
    if bname_list is None or pd.isna(bname_list):
        return "."
    s = str(bname_list).strip()
    if s in ["", ".", "NA", "NaN"]:
        return "."

    toks = [t.strip() for t in s.split(",") if t.strip()]
    types = []
    for t in toks:
        types.append(parse_bname_gene_fields(t)["gene_type"])
    types = sorted(set([x for x in types if x not in [".", "", "NA"]]))
    return "|".join(types) if types else "."


def collapse_domains(dom: pd.DataFrame) -> pd.DataFrame:
    required = ["Locus", "Gene", "DomainType", "Coverage"]
    missing = [c for c in required if c not in dom.columns]
    if missing:
        raise ValueError(f"Domain TSV missing required columns: {missing}")

    d = dom.copy()
    d["Coverage"] = pd.to_numeric(d["Coverage"], errors="coerce")

    for c in ["Gene", "DomainType", "Domain"]:
        if c in d.columns:
            d[c] = d[c].fillna(".").astype(str).str.strip()

    def fmt_gene(r):
        gene = r["Gene"]
        cov = r["Coverage"]
        return f"{gene}:{'NA' if pd.isna(cov) else f'{float(cov):.3f}'}"

    def fmt_type(r):
        gene = r["Gene"]
        dtype = r["DomainType"]
        cov = r["Coverage"]
        return f"{gene}|{dtype}:{'NA' if pd.isna(cov) else f'{float(cov):.3f}'}"

    def fmt_profile(r):
        gene = r["Gene"]
        dtype = r["DomainType"]
        domname = r["Domain"] if "Domain" in r.index else "."
        cov = r["Coverage"]
        return f"{gene}|{dtype}|{domname}:{'NA' if pd.isna(cov) else f'{float(cov):.3f}'}"

    d["entry_gene"] = d.apply(fmt_gene, axis=1)
    d["entry_type"] = d.apply(fmt_type, axis=1)
    d["entry_profile"] = d.apply(fmt_profile, axis=1) if "Domain" in d.columns else d["entry_type"]

    agg = d.groupby("Locus", as_index=False).agg(
        domains_gene=("entry_gene", lambda x: ";".join(sorted(set(x)))),
        domains_type=("entry_type", lambda x: ";".join(sorted(set(x)))),
        domains_profile=("entry_profile", lambda x: ";".join(sorted(set(x)))),
        domain_count=("entry_type", "count"),
        max_domain_cov=("Coverage", "max"),
        domain_genes_unique=("Gene", lambda x: pd.Series(x).nunique()),
        domain_types_unique=("DomainType", lambda x: pd.Series(x).nunique()),
    )
    return agg.rename(columns={"Locus": "locid"})


def parse_map_has_ltrs(mapdf: pd.DataFrame) -> pd.DataFrame:
    m = mapdf.copy()

    for c in ["ltr5_name", "ltr3_name"]:
        if c in m.columns:
            m[c] = m[c].fillna(".")
            m[c] = m[c].replace({"": ".", "NA": ".", "NaN": "."})

    def has_ltr(x):
        return 0 if (pd.isna(x) or str(x).strip() in ["", ".", "NA"]) else 1

    m["has_ltr5"] = m["ltr5_name"].apply(has_ltr) if "ltr5_name" in m.columns else 0
    m["has_ltr3"] = m["ltr3_name"].apply(has_ltr) if "ltr3_name" in m.columns else 0

    out_cols = ["locid", "ltr5_name", "ltr3_name", "has_ltr5", "has_ltr3"]
    for c in out_cols:
        if c not in m.columns:
            m[c] = "."
    return m[out_cols]


def extract_ltr_metrics(ltr: pd.DataFrame) -> pd.DataFrame:
    l = ltr.copy()

    for c in ["sequence_name", "LTR_length", "n_motifs", "status"]:
        if c not in l.columns:
            raise ValueError(f"LTR annotated TSV missing required column: {c}")

    l = coerce_numeric(l, ["LTR_length", "n_motifs"])
    l["conf_reconstructed"] = (l["status"].astype(str) == "OK").astype(int)

    return l[["sequence_name", "LTR_length", "n_motifs", "status", "conf_reconstructed"]].rename(
        columns={
            "sequence_name": "ltr_name",
            "LTR_length": "ltr_length",
            "n_motifs": "ltr_tfbm_burden",
            "status": "ltr_status",
        }
    )


def join_ltr_side(base: pd.DataFrame, ltr_metrics: pd.DataFrame, side: str) -> pd.DataFrame:
    key = f"ltr{side}_name"
    tmp = ltr_metrics.copy().rename(
        columns={
            "ltr_length": f"ltr{side}_length",
            "ltr_tfbm_burden": f"ltr{side}_tfbm_burden",
            "ltr_status": f"ltr{side}_status",
            "conf_reconstructed": f"ltr{side}_confident",
        }
    )
    merged = base.merge(tmp, how="left", left_on=key, right_on="ltr_name")
    merged = merged.drop(columns=["ltr_name"], errors="ignore")

    for c in [f"ltr{side}_length", f"ltr{side}_tfbm_burden", f"ltr{side}_confident"]:
        if c in merged.columns:
            merged[c] = merged[c].fillna(0).astype(int)
    for c in [f"ltr{side}_status"]:
        if c in merged.columns:
            merged[c] = merged[c].fillna(".")

    return merged


def collapse_rbp(fimo: pd.DataFrame, qval_thresh: float) -> pd.DataFrame:
    f = fimo.copy()
    if "sequence_name" not in f.columns:
        raise ValueError("RBP fimo.tsv missing 'sequence_name'")
    if "q-value" in f.columns:
        f["q-value"] = pd.to_numeric(f["q-value"], errors="coerce")
        f = f[(f["q-value"].isna()) | (f["q-value"] <= qval_thresh)]

    agg = f.groupby("sequence_name", as_index=False).agg(
        rbp_burden=("motif_alt_id", "count") if "motif_alt_id" in f.columns else ("motif_id", "count"),
        rbp_unique=("motif_alt_id", lambda x: pd.Series(x).nunique()) if "motif_alt_id" in f.columns else ("motif_id", lambda x: pd.Series(x).nunique()),
    )
    return agg.rename(columns={"sequence_name": "locid"})



def parse_closest_table(path: str, out_prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, comment=None)
    df = df.rename(columns={c: c.lstrip("#") for c in df.columns})

    required = ["a_name", "b_name", "distance_bp"]
    for r in required:
        if r not in df.columns:
            raise ValueError(f"Closest TSV {path} missing column {r}")

    df["distance_bp"] = pd.to_numeric(df["distance_bp"], errors="coerce")

    min_dist = df.groupby("a_name", as_index=False)["distance_bp"].min().rename(columns={"distance_bp": "min_dist"})
    df2 = df.merge(min_dist, on="a_name", how="left")
    df2 = df2[df2["distance_bp"] == df2["min_dist"]].copy()

    parsed = df2["b_name"].apply(parse_bname_gene_fields)
    df2[f"{out_prefix}_gene_id"] = parsed.apply(lambda x: x["gene_id"])
    df2[f"{out_prefix}_gene_symbol"] = parsed.apply(lambda x: x["gene_symbol"])
    df2[f"{out_prefix}_gene_type"] = parsed.apply(lambda x: x["gene_type"])

    agg = df2.groupby("a_name", as_index=False).agg(
        **{
            f"{out_prefix}_name": ("b_name", lambda x: "|".join(sorted(set(map(str, x))))),
            f"{out_prefix}_distance_bp": ("min_dist", "first"),
            f"{out_prefix}_gene_id": (f"{out_prefix}_gene_id", lambda x: "|".join(sorted(set(map(str, x))))),
            f"{out_prefix}_gene_symbol": (f"{out_prefix}_gene_symbol", lambda x: "|".join(sorted(set(map(str, x))))),
            f"{out_prefix}_gene_type": (f"{out_prefix}_gene_type", lambda x: "|".join(sorted(set(map(str, x))))),
        }
    ).rename(columns={"a_name": "locid"})
    return agg


def collapse_overlaps(overlap_files: Dict[str, str]) -> pd.DataFrame:
    """
    Build overlap flags per locid and:
      - ov_gene, ov_exon, ov_intron, ov_cds, ov_utr (0/1)
      - ov_gene_names (pipe-separated b_name strings from gene overlaps)
      - ov_gene_types (union of biotypes from gene overlaps)
      - feature_overlap_per_gene: "ENSG...:intron;ENSG...:exon|intron"
      - feature_overlap: "intron;exon|intron"
        If no overlaps -> "intergenic"
    """
    feats = ["gene", "exon", "intron", "cds", "utr"]
    frames = []

    for feat in feats:
        path = (overlap_files.get(feat, "") or "").strip()
        if not path:
            continue
        if not os.path.exists(path):
            print(f"[WARN] overlap file not found for {feat}: {path}", file=sys.stderr)
            continue

        df = pd.read_csv(path, sep="\t", dtype=str)
        df = df.rename(columns={c: c.lstrip("#") for c in df.columns})

        if "a_name" not in df.columns or "b_name" not in df.columns:
            raise ValueError(f"Overlap TSV {path} missing a_name/b_name. Columns={df.columns.tolist()[:20]}")

        tmp = df[["a_name", "b_name"]].copy()
        tmp = tmp.rename(columns={"a_name": "locid", "b_name": "b_name"})
        tmp["feat"] = feat
        frames.append(tmp)

    if not frames:
        out = pd.DataFrame({"locid": []})
        for feat in feats:
            out[f"ov_{feat}"] = []
        out["ov_gene_names"] = []
        out["ov_gene_types"] = []
        out["feature_overlap_per_gene"] = []
        out["feature_overlap"] = []
        return out

    allx = pd.concat(frames, ignore_index=True)

    flags = (
        allx.groupby(["locid", "feat"], as_index=False)
            .size()
            .assign(val=1)
            .pivot(index="locid", columns="feat", values="val")
            .fillna(0)
            .astype(int)
            .reset_index()
    )

    for feat in feats:
        if feat not in flags.columns:
            flags[feat] = 0
    flags = flags.rename(columns={feat: f"ov_{feat}" for feat in feats})

    gene_names = (
        allx[allx["feat"] == "gene"]
        .groupby("locid", as_index=False)
        .agg(ov_gene_names=("b_name", lambda x: "|".join(sorted(set(map(str, x))))))
    )

    out = flags.merge(gene_names, how="left", on="locid")
    out["ov_gene_names"] = out["ov_gene_names"].fillna(".")

    gene_types = (
        allx[allx["feat"] == "gene"]
        .assign(gene_type=lambda d: d["b_name"].apply(lambda x: parse_bname_gene_fields(x)["gene_type"]))
        .groupby("locid", as_index=False)
        .agg(ov_gene_types=("gene_type", lambda x: "|".join(sorted(set([str(v) for v in x if str(v) not in [".", "", "NA"]])))))
    )
    out = out.merge(gene_types, how="left", on="locid")
    out["ov_gene_types"] = out["ov_gene_types"].fillna(".")

    def extract_gene_ids(feat: str, bname: str) -> List[str]:
        if bname is None or pd.isna(bname):
            return []
        s = str(bname)

        if feat == "exon":
            toks = [t.strip() for t in s.split(",") if t.strip()]
            return sorted(set([t for t in toks if t.startswith("ENSG")]))

        parts = s.split("|")
        genes = [p for p in parts if p.startswith("ENSG")]
        if not genes and s.startswith("ENSG"):
            genes = [s.split("|")[0]]
        return sorted(set(genes))

    gene_feat_rows = []
    for locid, feat, bname in allx[["locid", "feat", "b_name"]].itertuples(index=False, name=None):
        for gid in extract_gene_ids(feat, bname):
            gene_feat_rows.append((locid, gid, feat))

    def fallback_class_row(r):
        if r.get("ov_cds", 0) == 1: return "cds"
        if r.get("ov_utr", 0) == 1: return "utr"
        if r.get("ov_exon", 0) == 1: return "exon"
        if r.get("ov_intron", 0) == 1: return "intron"
        if r.get("ov_gene", 0) == 1: return "gene"
        return "intergenic"

    fallback_feature_overlap = out.apply(fallback_class_row, axis=1)

    if gene_feat_rows:
        gdf = pd.DataFrame(gene_feat_rows, columns=["locid", "gene_id", "feat"]).drop_duplicates()

        feat_order = {"cds": 0, "utr": 1, "exon": 2, "intron": 3, "gene": 4}

        per_gene_featset = (
            gdf.groupby(["locid", "gene_id"])["feat"]
               .apply(lambda x: "|".join(sorted(set(x), key=lambda f: feat_order.get(f, 99))))
               .reset_index(name="featset")
        )

        per_gene = (
            per_gene_featset.sort_values(["locid", "gene_id"])
            .groupby("locid")
            .apply(lambda blk: ";".join(
                f"{gid}:{fs}" for gid, fs in blk[["gene_id", "featset"]].itertuples(index=False, name=None)
            ))
            .rename("feature_overlap_per_gene")
            .reset_index()
        )

        per_locid = (
            per_gene_featset.sort_values(["locid", "gene_id"])
            .groupby("locid")
            .apply(lambda blk: ";".join(pd.unique(blk["featset"])))
            .rename("feature_overlap")
            .reset_index()
        )

        out = out.merge(per_gene, how="left", on="locid")
        out = out.merge(per_locid, how="left", on="locid")

        out["feature_overlap_per_gene"] = out["feature_overlap_per_gene"].fillna(".")
        out["feature_overlap"] = out["feature_overlap"].fillna(fallback_feature_overlap).fillna("intergenic")
    else:
        out["feature_overlap_per_gene"] = "."
        out["feature_overlap"] = fallback_feature_overlap.fillna("intergenic")

    if "feature_overlap" not in out.columns:
        out["feature_overlap"] = "intergenic"

    return out



def only_same_strand_overlap_columns(ov: pd.DataFrame) -> pd.DataFrame:
    """
    Take a collapse_overlaps() output and keep only:
      - feature_overlap_same_strand
      - feature_overlap_per_gene_same_strand
      - ov_gene_types_same_strand
    """
    if ov is None or ov.shape[0] == 0:
        return pd.DataFrame(columns=[
            "locid",
            "feature_overlap_same_strand",
            "feature_overlap_per_gene_same_strand",
            "ov_gene_types_same_strand",
        ])

    keep = ["locid", "feature_overlap", "feature_overlap_per_gene", "ov_gene_types"]
    keep = [c for c in keep if c in ov.columns]
    out = ov[keep].copy()

    out = out.rename(columns={
        "feature_overlap": "feature_overlap_same_strand",
        "feature_overlap_per_gene": "feature_overlap_per_gene_same_strand",
        "ov_gene_types": "ov_gene_types_same_strand",
    })

    # defaults
    if "feature_overlap_same_strand" in out.columns:
        out["feature_overlap_same_strand"] = out["feature_overlap_same_strand"].fillna("intergenic")
    else:
        out["feature_overlap_same_strand"] = "intergenic"

    if "feature_overlap_per_gene_same_strand" in out.columns:
        out["feature_overlap_per_gene_same_strand"] = out["feature_overlap_per_gene_same_strand"].fillna(".")
    else:
        out["feature_overlap_per_gene_same_strand"] = "."

    if "ov_gene_types_same_strand" in out.columns:
        out["ov_gene_types_same_strand"] = out["ov_gene_types_same_strand"].fillna(".")
    else:
        out["ov_gene_types_same_strand"] = "."

    return out


def main():
    ap = argparse.ArgumentParser(description="Build integrated ncRNA-layer table for internal HERV regions.")
    ap.add_argument("--internal-bed", required=True, help="internal_regions.sorted.bed.gz (BED6)")
    ap.add_argument("--domains-tsv", required=True, help="ERV_GyDB_v6_filtered.tsv")
    ap.add_argument("--map-tsv", required=True, help="ERV_full_plus_components.map.tsv")
    ap.add_argument("--ltr-annot-tsv", required=True, help="LTR_fully_annotated.tsv")
    ap.add_argument("--rbp-fimo", required=True, help="rbp/fimo/internal/fimo.tsv")
    ap.add_argument("--closest-gene-abs", required=True, help="internal_closest_gene.tsv")
    ap.add_argument("--closest-gene-same", required=True, help="internal_closest_gene_same_strand.tsv")
    ap.add_argument("--closest-tss-abs", required=True, help="internal_closest_tss.tsv")
    ap.add_argument("--closest-tss-same", required=True, help="internal_closest_tss_same_strand.tsv")

    # existing (any-strand) overlaps
    ap.add_argument("--overlap-genes", required=False, default="", help="internal_vs_genes.tsv")
    ap.add_argument("--overlap-exons", required=False, default="", help="internal_vs_gene_exons.tsv")
    ap.add_argument("--overlap-introns", required=False, default="", help="internal_vs_gene_introns.tsv")
    ap.add_argument("--overlap-cds", required=False, default="", help="internal_vs_cds.tsv")
    ap.add_argument("--overlap-utr", required=False, default="", help="internal_vs_utr.tsv")

    # NEW: same-strand overlaps (generated with bedtools intersect -s)
    ap.add_argument("--overlap-genes-same", required=False, default="", help="internal_vs_genes.same_strand.tsv")
    ap.add_argument("--overlap-exons-same", required=False, default="", help="internal_vs_gene_exons.same_strand.tsv")
    ap.add_argument("--overlap-introns-same", required=False, default="", help="internal_vs_gene_introns.same_strand.tsv")
    ap.add_argument("--overlap-cds-same", required=False, default="", help="internal_vs_cds.same_strand.tsv")
    ap.add_argument("--overlap-utr-same", required=False, default="", help="internal_vs_utr.same_strand.tsv")

    ap.add_argument("--rbp-qval", type=float, default=1,
                    help="Optional FIMO q-value filter; set to 1 (default) to disable.")
    ap.add_argument("--out", required=True, help="Output TSV path")

    args = ap.parse_args()

    # ---- base internal bed
    base = read_bed(args.internal_bed)
    base["length"] = base["end"] - base["start"]

    # ---- domains
    dom = safe_read_tsv(args.domains_tsv)
    dom_agg = collapse_domains(dom)
    base = base.merge(dom_agg, how="left", on="locid")
    for c in ["domains_gene", "domains_type", "domains_profile"]:
        if c in base.columns:
            base[c] = base[c].fillna(".")
    base["domain_count"] = base["domain_count"].fillna(0).astype(int)

    # ---- LTR map + metrics
    mapdf = safe_read_tsv(args.map_tsv)
    map_sel = parse_map_has_ltrs(mapdf)
    base = base.merge(map_sel, how="left", on="locid")
    for c in ["has_ltr5", "has_ltr3"]:
        base[c] = base[c].fillna(0).astype(int)
    for c in ["ltr5_name", "ltr3_name"]:
        base[c] = base[c].fillna(".")

    ltr = safe_read_tsv(args.ltr_annot_tsv)
    ltr_metrics = extract_ltr_metrics(ltr)
    base = join_ltr_side(base, ltr_metrics, "5")
    base = join_ltr_side(base, ltr_metrics, "3")

    # ---- RBP burden
    rbp = safe_read_tsv(args.rbp_fimo)
    rbp_agg = collapse_rbp(rbp, args.rbp_qval)
    base = base.merge(rbp_agg, how="left", on="locid")
    base["rbp_burden"] = base["rbp_burden"].fillna(0).astype(int)
    base["rbp_unique"] = base["rbp_unique"].fillna(0).astype(int)


    # ---- closest gene/tss
    cg_abs = parse_closest_table(args.closest_gene_abs, "closest_gene_abs")
    cg_same = parse_closest_table(args.closest_gene_same, "closest_gene_same_strand")
    ct_abs = parse_closest_table(args.closest_tss_abs, "closest_tss_abs")
    ct_same = parse_closest_table(args.closest_tss_same, "closest_tss_same_strand")

    for d in [cg_abs, cg_same, ct_abs, ct_same]:
        base = base.merge(d, how="left", on="locid")

    # ---- overlaps (any-strand)
    overlap_files = {
        "gene": args.overlap_genes,
        "exon": args.overlap_exons,
        "intron": args.overlap_introns,
        "cds": args.overlap_cds,
        "utr": args.overlap_utr,
    }
    ov = collapse_overlaps(overlap_files)
    base = base.merge(ov, how="left", on="locid")

    for c in ["ov_gene", "ov_exon", "ov_intron", "ov_cds", "ov_utr"]:
        if c in base.columns:
            base[c] = base[c].fillna(0).astype(int)

    if "feature_overlap" in base.columns:
        base["feature_overlap"] = base["feature_overlap"].fillna("intergenic")
    else:
        base["feature_overlap"] = "intergenic"

    if "feature_overlap_per_gene" in base.columns:
        base["feature_overlap_per_gene"] = base["feature_overlap_per_gene"].fillna(".")
    else:
        base["feature_overlap_per_gene"] = "."

    if "ov_gene_types" in base.columns:
        base["ov_gene_types"] = base["ov_gene_types"].fillna(".")
    else:
        base["ov_gene_types"] = "."

    # ---- overlaps (SAME STRAND) -> NEW requested columns
    overlap_files_same = {
        "gene": args.overlap_genes_same,
        "exon": args.overlap_exons_same,
        "intron": args.overlap_introns_same,
        "cds": args.overlap_cds_same,
        "utr": args.overlap_utr_same,
    }
    ov_same_full = collapse_overlaps(overlap_files_same)
    ov_same = only_same_strand_overlap_columns(ov_same_full)
    base = base.merge(ov_same, how="left", on="locid")

    # defaults
    base["feature_overlap_same_strand"] = base.get("feature_overlap_same_strand", "intergenic")
    base["feature_overlap_same_strand"] = base["feature_overlap_same_strand"].fillna("intergenic")

    base["feature_overlap_per_gene_same_strand"] = base.get("feature_overlap_per_gene_same_strand", ".")
    base["feature_overlap_per_gene_same_strand"] = base["feature_overlap_per_gene_same_strand"].fillna(".")

    base["ov_gene_types_same_strand"] = base.get("ov_gene_types_same_strand", ".")
    base["ov_gene_types_same_strand"] = base["ov_gene_types_same_strand"].fillna(".")

    # ---- select output columns
    out_cols = [
        "chrom", "start", "end", "strand", "length", "locid",
        "domains_gene", "domains_type", "domains_profile", "domain_count", "max_domain_cov",
        "has_ltr5", "has_ltr3",
        "ltr5_name", "ltr5_length", "ltr5_tfbm_burden", "ltr5_status",
        "ltr3_name", "ltr3_length", "ltr3_tfbm_burden", "ltr3_status",
        "rbp_burden", "rbp_unique",
        "closest_gene_abs_name", "closest_gene_abs_gene_type",
        "closest_gene_abs_distance_bp",
        "closest_gene_same_strand_name", "closest_gene_same_strand_gene_type",
        "closest_gene_same_strand_distance_bp",
        "closest_tss_abs_name", "closest_tss_abs_gene_type",
        "closest_tss_abs_distance_bp",
        "closest_tss_same_strand_gene_type",
        "closest_tss_same_strand_name",
        "closest_tss_same_strand_distance_bp",

        # ANY-strand overlaps
        "feature_overlap",
        "feature_overlap_per_gene",
        "ov_gene_types",

        # SAME-strand overlaps
        "feature_overlap_same_strand",
        "feature_overlap_per_gene_same_strand",
        "ov_gene_types_same_strand",
    ]
    out_cols = [c for c in out_cols if c in base.columns]
    out = base[out_cols].copy()

    for c in ["has_ltr5", "has_ltr3", "ltr5_confident", "ltr3_confident"]:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype(int)

    out.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] Wrote: {args.out} ({out.shape[0]} rows, {out.shape[1]} cols)", file=sys.stderr)


if __name__ == "__main__":
    main()
