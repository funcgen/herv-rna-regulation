#!/usr/bin/env python3
"""
03_transcript_context step 9 (frozen global order 16, effective order 17):
the integrated locus-level table -- one row per HERV internal region, 41 columns.

Refactor of scripts/integration/01_build_integrated_internal_table_simplify.py as
invoked by run_build_internal_table.sh (the full 19-argument invocation, D010).

What is refactored and what is not
  * Every aggregation function (read_bed ... only_same_strand_overlap_columns) is
    carried over VERBATIM from the historical module, between the markers below.
    They are the definition of the table's semantics -- delimiters, sorting,
    fill values, the per-gene key quirks -- and downstream parsers depend on them.
    Block sha256 (historical, identical): 79e75b01d5916f414ddad46190a2e1dd926d98dd689fd840eef468b9d186e443
  * main() is new: explicit canonical inputs, provenance guards, and a
    one-row-per-locus invariant after every join. Merge order, fill values and
    the 41-column selection are the historical main()'s, unchanged.

Canonical inputs (accepted release artifacts only)
  --internal-bed     03 step 5  results/herv_internal/internal_regions.sorted.bed.gz  (53,176 loci)
  --domains-tsv      SRC06      01_sources/frozen/herv/HERV_loci_annotated_domains_v5.tsv
  --map-tsv          SRC10      01_sources/frozen/ltr/HERV_internal_LTR_map.tsv
  --ltr-projection   07_sparcs  results/ltr_structure/herv_ltr_minimal_projection.tsv.gz (D017)
                     -- replaces the unpublished LTR_fully_annotated.tsv (SRC11, refused)
  --rbp-fimo         02_rbp stage 3  results/fimo_primary/fimo.tsv -- the RAW FIMO output (D018);
                     a stage-4 parsed table is refused; --rbp-qval 1 as historically
  --closest-dir      03 step 6  results/closest/                  (D023 TSS used there)
  --overlap-dir      03 step 5  results/overlap/                  (D022 per-gene exons/introns)

Output: one TSV, historical header, 53,176 rows x 41 columns, base (sorted BED)
order. No HERV_id (effective order 18). No field from the older February atlas.
"""
import argparse
import gzip
import hashlib
import os
import sys
import tempfile
import pandas as pd
from typing import Dict, List, Optional

N_LOCI = 53176
PINNED = {  # accepted identities of the inputs that live outside 03_transcript_context
    "domains_tsv": ("md5", "162cba73027d1e4500dc1dbb2671c84a"),      # SRC06 (01_sources)
    "map_tsv": ("md5", "a0aa5d8cdcac58bf719d0df4bf07f3fe"),          # SRC10 (01_sources)
    "rbp_fimo": ("sha256", "46b4685eb160786946dd9bb80cfecd9d1519e5807538c9e4608905404c6a513b"),  # 02_rbp stage 3
    "ltr_projection": ("sha256", "2e06a23a5a1382a82a18e81a2f59fbe7ac6eefdc5d16572fa747ca28bf5d9c42"),  # 07_sparcs, D017
}
RAW_FIMO_HEADER = "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence"

# ============================================================================
# BEGIN verbatim historical functions (01_build_integrated_internal_table_simplify.py)
# ============================================================================
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
# ============================================================================
# END verbatim historical functions
# ============================================================================


def die(msg):
    sys.exit(f"[error] {msg}")


def one_row_per_locus(df: pd.DataFrame, after: str) -> pd.DataFrame:
    if len(df) != N_LOCI or df["locid"].nunique() != N_LOCI:
        die(f"one-row-per-locus invariant broken after {after}: {len(df)} rows, {df['locid'].nunique()} unique loci")
    return df


def digest(path: str, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def verify_pinned(key: str, path: str):
    algo, want = PINNED[key]
    got = digest(path, algo)
    if got != want:
        die(f"{key}: {path} {algo} {got} != accepted {want}")
    print(f"[ok] {key}: {algo} {got}", file=sys.stderr)


def verify_in_manifest(path: str):
    """03 step outputs: verify against the checksums.sha256 written next to them."""
    d, name = os.path.split(os.path.abspath(path))
    want = None
    with open(os.path.join(d, "checksums.sha256")) as fh:
        for line in fh:
            h, f = line.rstrip("\n").split(None, 1)
            if f.lstrip("*").removeprefix("./") == name:
                want = h
    if want is None:
        die(f"{path} is not listed in {d}/checksums.sha256")
    got = digest(path, "sha256")
    if got != want:
        die(f"{path} sha256 {got} != manifest {want}")
    print(f"[ok] {name}: sha256 {got} (manifest)", file=sys.stderr)


def guard_raw_fimo(path: str):
    """D018: the integration input must be the RAW stage-3 fimo.tsv."""
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
    if first != RAW_FIMO_HEADER:
        die(f"--rbp-fimo is not a raw FIMO fimo.tsv (header {first[:80]!r}); D018 requires the raw "
            "stage-3 output, not the stage-4 parsed table")
    with open(path, "rb") as fh:
        fh.seek(max(0, os.path.getsize(path) - 4096))
        if b"# FIMO (Find Individual Motif Occurrences)" not in fh.read():
            die("--rbp-fimo lacks the FIMO trailer; D018 requires the unmodified raw output")


def main():
    ap = argparse.ArgumentParser(description="integrated HERV-locus table (effective order 17)")
    ap.add_argument("--internal-bed", required=True)
    ap.add_argument("--domains-tsv", required=True)
    ap.add_argument("--map-tsv", required=True)
    ap.add_argument("--ltr-projection", required=True)
    ap.add_argument("--rbp-fimo", required=True)
    ap.add_argument("--closest-dir", required=True)
    ap.add_argument("--overlap-dir", required=True)
    ap.add_argument("--rbp-qval", type=float, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    for p in (args.internal_bed, args.domains_tsv, args.map_tsv, args.ltr_projection, args.rbp_fimo,
              args.closest_dir, args.overlap_dir):
        if "LTR_fully_annotated" in p:
            die("LTR_fully_annotated.tsv (SRC11) is excluded (D014); use the D017 minimal projection")
        if "historical_reproduction" in p:
            die("historical-reproduction (forensic) inputs are refused")
    guard_raw_fimo(args.rbp_fimo)
    for key in PINNED:
        verify_pinned(key, getattr(args, key))
    C, O = args.closest_dir, args.overlap_dir
    verify_in_manifest(args.internal_bed)
    for fn in ("internal_closest_gene.tsv", "internal_closest_gene_same_strand.tsv",
               "internal_closest_tss.tsv", "internal_closest_tss_same_strand.tsv"):
        verify_in_manifest(f"{C}/{fn}")
    for t in ("genes", "gene_exons", "gene_introns", "cds", "utr"):
        for sfx in (".tsv", ".same_strand.tsv"):
            verify_in_manifest(f"{O}/internal_vs_{t}{sfx}")

    # ---- base internal bed (universe and order)
    base = read_bed(args.internal_bed)
    base["length"] = base["end"] - base["start"]
    one_row_per_locus(base, "base")

    # ---- domains (SRC06)
    dom = safe_read_tsv(args.domains_tsv)
    dom_agg = collapse_domains(dom)
    base = one_row_per_locus(base.merge(dom_agg, how="left", on="locid"), "domains")
    for c in ["domains_gene", "domains_type", "domains_profile"]:
        if c in base.columns:
            base[c] = base[c].fillna(".")
    base["domain_count"] = base["domain_count"].fillna(0).astype(int)

    # ---- LTR map (SRC10) + minimal LTR projection (D017)
    mapdf = safe_read_tsv(args.map_tsv)
    map_sel = parse_map_has_ltrs(mapdf)
    base = one_row_per_locus(base.merge(map_sel, how="left", on="locid"), "LTR map")
    for c in ["has_ltr5", "has_ltr3"]:
        base[c] = base[c].fillna(0).astype(int)
    for c in ["ltr5_name", "ltr3_name"]:
        base[c] = base[c].fillna(".")
    ltr = safe_read_tsv(args.ltr_projection)
    ltr_metrics = extract_ltr_metrics(ltr)
    base = one_row_per_locus(join_ltr_side(base, ltr_metrics, "5"), "5' LTR metrics")
    base = one_row_per_locus(join_ltr_side(base, ltr_metrics, "3"), "3' LTR metrics")

    # ---- RBP burden (raw stage-3 FIMO, D018)
    rbp = safe_read_tsv(args.rbp_fimo)
    rbp_agg = collapse_rbp(rbp, args.rbp_qval)
    base = one_row_per_locus(base.merge(rbp_agg, how="left", on="locid"), "RBP")
    base["rbp_burden"] = base["rbp_burden"].fillna(0).astype(int)
    base["rbp_unique"] = base["rbp_unique"].fillna(0).astype(int)

    # ---- closest gene / TSS (step 6)
    for fn, prefix in (("internal_closest_gene.tsv", "closest_gene_abs"),
                       ("internal_closest_gene_same_strand.tsv", "closest_gene_same_strand"),
                       ("internal_closest_tss.tsv", "closest_tss_abs"),
                       ("internal_closest_tss_same_strand.tsv", "closest_tss_same_strand")):
        base = one_row_per_locus(base.merge(parse_closest_table(f"{C}/{fn}", prefix), how="left", on="locid"), prefix)

    # ---- overlaps, any strand (step 5)
    overlap_files = {"gene": f"{O}/internal_vs_genes.tsv", "exon": f"{O}/internal_vs_gene_exons.tsv",
                     "intron": f"{O}/internal_vs_gene_introns.tsv", "cds": f"{O}/internal_vs_cds.tsv",
                     "utr": f"{O}/internal_vs_utr.tsv"}
    for p in overlap_files.values():
        if not os.path.exists(p):
            die(f"missing overlap table {p}")
    ov = collapse_overlaps(overlap_files)
    base = one_row_per_locus(base.merge(ov, how="left", on="locid"), "any-strand overlaps")
    for c in ["ov_gene", "ov_exon", "ov_intron", "ov_cds", "ov_utr"]:
        if c in base.columns:
            base[c] = base[c].fillna(0).astype(int)
    base["feature_overlap"] = base["feature_overlap"].fillna("intergenic") if "feature_overlap" in base.columns else "intergenic"
    base["feature_overlap_per_gene"] = base["feature_overlap_per_gene"].fillna(".") if "feature_overlap_per_gene" in base.columns else "."
    base["ov_gene_types"] = base["ov_gene_types"].fillna(".") if "ov_gene_types" in base.columns else "."

    # ---- overlaps, same strand
    overlap_files_same = {k: v.replace(".tsv", ".same_strand.tsv") for k, v in overlap_files.items()}
    for p in overlap_files_same.values():
        if not os.path.exists(p):
            die(f"missing overlap table {p}")
    ov_same = only_same_strand_overlap_columns(collapse_overlaps(overlap_files_same))
    base = one_row_per_locus(base.merge(ov_same, how="left", on="locid"), "same-strand overlaps")
    base["feature_overlap_same_strand"] = base.get("feature_overlap_same_strand", "intergenic")
    base["feature_overlap_same_strand"] = base["feature_overlap_same_strand"].fillna("intergenic")
    base["feature_overlap_per_gene_same_strand"] = base.get("feature_overlap_per_gene_same_strand", ".")
    base["feature_overlap_per_gene_same_strand"] = base["feature_overlap_per_gene_same_strand"].fillna(".")
    base["ov_gene_types_same_strand"] = base.get("ov_gene_types_same_strand", ".")
    base["ov_gene_types_same_strand"] = base["ov_gene_types_same_strand"].fillna(".")

    # ---- the historical 41-column selection (no HERV_id; no February-atlas fields)
    out_cols = [
        "chrom", "start", "end", "strand", "length", "locid",
        "domains_gene", "domains_type", "domains_profile", "domain_count", "max_domain_cov",
        "has_ltr5", "has_ltr3",
        "ltr5_name", "ltr5_length", "ltr5_tfbm_burden", "ltr5_status",
        "ltr3_name", "ltr3_length", "ltr3_tfbm_burden", "ltr3_status",
        "rbp_burden", "rbp_unique",
        "closest_gene_abs_name", "closest_gene_abs_gene_type", "closest_gene_abs_distance_bp",
        "closest_gene_same_strand_name", "closest_gene_same_strand_gene_type", "closest_gene_same_strand_distance_bp",
        "closest_tss_abs_name", "closest_tss_abs_gene_type", "closest_tss_abs_distance_bp",
        "closest_tss_same_strand_gene_type", "closest_tss_same_strand_name", "closest_tss_same_strand_distance_bp",
        "feature_overlap", "feature_overlap_per_gene", "ov_gene_types",
        "feature_overlap_same_strand", "feature_overlap_per_gene_same_strand", "ov_gene_types_same_strand",
    ]
    missing = [c for c in out_cols if c not in base.columns]
    if missing:
        die(f"columns missing from the integrated table: {missing}")
    out = base[out_cols].copy()
    for c in ["has_ltr5", "has_ltr3"]:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype(int)
    one_row_per_locus(out, "column selection")

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".stage_integrated.", dir=outdir)
    os.close(fd)
    try:
        out.to_csv(tmp, sep="\t", index=False)
        os.chmod(tmp, 0o644)
        os.replace(tmp, args.out)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    print(f"[OK] Wrote: {args.out} ({out.shape[0]} rows, {out.shape[1]} cols)", file=sys.stderr)


if __name__ == "__main__":
    main()
