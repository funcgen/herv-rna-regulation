#!/usr/bin/env python3
"""
07_sparcs step 1 (frozen global order 42, effective order 43), CNAG: antisense terminal-LTR candidates.

REFACTOR_SCRIPT of ncRNA_layer/scripts/integration/05_ltr_features_overlap_tables.py (sha256 7abb1e3f...).
The scientific rules are the historical ones; parse_bname_tx_feature, _merge_intervals,
union_overlap_bp_from_group and the aggregation in main() are carried over verbatim:

  input          the five accepted order-15 tables (03 results/overlap_ltrs/), each
                 `bedtools intersect -sorted -g -S -wo` of the SRC07 LTRs against one 03 feature BED
  antisense      LTR strand != strand of the transcript feature (the -S of order 15; nothing is
                 re-filtered here, a guard only checks it)
  feature class  exon (transcript exons), CDS, the 5'/3' label of the UTR split, three_prime_UTR
                 (utr3 table), last_exon (per-transcript terminal exon of 03 step 4)
  overlap        overlap = min(ends) - max(starts), BED half-open; kept when >= sparcs.candidate_min_overlap_bp
                 (1: any positive overlap; no LTR or feature fraction threshold)
  row unit       LTR x transcript x feature_class; overlap_bp_total = union of the overlapped intervals
                 (so the utr3 and utr_5p3p copies of one 3' UTR relationship count once, and an LTR
                 spanning two exons of one transcript is one row); overlap_frac_of_ltr = total / LTR length
  candidates     feature_class in {three_prime_UTR, last_exon}, every gene/transcript biotype kept
  pc/lnc subset  candidates with gene_type in {protein_coding, lncRNA}
  branches       (historically applied by the consumer E_sparcs/01_analysis_sparcs.R; written here as an
                 explicit table) coding_3UTR = gene_type protein_coding & feature_class three_prime_UTR;
                 lncRNA_last_exon = gene_type lncRNA (every such candidate row is a last_exon row: guarded)

Release changes (declared before production):
  D033  read_ov() parsed each table with pd.read_csv(comment="#"): the '#a_chrom ...' header was
        discarded as a comment, the first DATA row became the header and was overwritten, so every
        table lost one overlap row. The header is now parsed as the header (and must be exactly the
        historical one); no data row is dropped.
  guards inputs must match their accepted manifests (never the historical tree); every LTR must be an
        SRC07 record with identical coordinates and strand; every transcript and gene must be a GENCODE
        v48 record of 03 step 1 with identical attributes and strand; valid half-open coordinates,
        overlap_bp recomputed, strands in {+,-} and opposite, no duplicate input rows.
  order historical sort keys first, then a_chrom, a_start, a_end, a_name, tx_id: a total order, so ties
        no longer depend on the sort algorithm (the historical quicksort left them unspecified).
  outputs staged and published atomically, with checksums.sha256, a branch table and a QC summary.
"""

import argparse
import gzip
import hashlib
import os
import re
import shutil
import sys
import tempfile
from typing import List, Tuple

import pandas as pd


EXPECTED = ["a_chrom", "a_start", "a_end", "a_name", "a_score", "a_strand",
            "b_chrom", "b_start", "b_end", "b_name", "b_score", "b_strand", "overlap_bp"]
TABLES = {   # argument -> (order-15 file, allowed FEATURE field(s) of b_name)
    "exon": ("ltrs_vs_gencode_tx_exons.antisense.wo.tsv", {"exon"}),
    "cds": ("ltrs_vs_gencode_cds.antisense.wo.tsv", {"CDS"}),
    "utr5p3p": ("ltrs_vs_gencode_utr_5p3p.antisense.wo.tsv", {"five_prime_UTR", "three_prime_UTR"}),
    "utr3": ("ltrs_vs_gencode_utr3.antisense.wo.tsv", {"three_prime_UTR"}),
    "last_exon": ("ltrs_vs_gencode_last_exon.antisense.wo.tsv", {"last_exon"}),
}
TX_RE = re.compile(r"^ENST\d{11}\.\d+(_PAR_Y)?$")
GENE_RE = re.compile(r"^ENSG\d{11}\.\d+(_PAR_Y)?$")
BRANCHES = ("coding_3UTR", "lncRNA_last_exon")
TIE = ["a_chrom", "a_start", "a_end", "a_name", "tx_id"]


# ---------------------------------------------------------------- historical functions (verbatim)
def parse_bname_tx_feature(bname: str):
    """
    Expected b_name format:
      ENST|ENSG|GENE|FEATURE|GENE_TYPE|TX_TYPE
    (same as your pipeline)
    """
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
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def union_overlap_bp_from_group(g: pd.DataFrame) -> int:
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
    # RELEASE (D033): the header is parsed as the header. Historically comment="#" discarded it and the
    # first data row was consumed as the header. keep_default_na=False: no field is ever read as NaN.
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    expected = ["a_chrom","a_start","a_end","a_name","a_score","a_strand",
                "b_chrom","b_start","b_end","b_name","b_score","b_strand","overlap_bp"]
    if df.shape[1] != 13:
        raise ValueError(f"{path}: expected 13 columns, got {df.shape[1]}")
    if list(df.columns) != ["#" + expected[0]] + expected[1:]:
        raise ValueError(f"{path}: header is not the historical 13-column '#a_chrom ... overlap_bp' header")
    df.columns = expected

    for c in ["a_start","a_end","b_start","b_end","overlap_bp"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Overlapped genomic interval (half-open)
    df["ovl_start"] = df[["a_start", "b_start"]].max(axis=1)
    df["ovl_end"]   = df[["a_end", "b_end"]].min(axis=1)
    df["ovl_len"]   = (df["ovl_end"] - df["ovl_start"]).clip(lower=0)

    # Parse feature fields
    parsed = df["b_name"].apply(parse_bname_tx_feature)
    df["tx_id"] = parsed.apply(lambda x: x[0])
    df["gene_id"] = parsed.apply(lambda x: x[1])
    df["gene_name"] = parsed.apply(lambda x: x[2])
    df["feature"] = parsed.apply(lambda x: x[3])
    df["gene_type"] = parsed.apply(lambda x: x[4])
    df["tx_type"] = parsed.apply(lambda x: x[5])

    # LTR length
    df["ltr_len_bp"] = (df["a_end"] - df["a_start"]).astype("float")

    return df


# ---------------------------------------------------------------- release: parameters, manifests, guards
def read_min_overlap(path: str) -> int:
    """sparcs.candidate_min_overlap_bp of config/parameters.yaml (no YAML module on CNAG)."""
    text = open(path, encoding="utf-8").read()
    m = re.search(r"^sparcs:\s*\n(.*?)(?=^\S|\Z)", text, re.M | re.S)
    k = re.search(r"^\s+candidate_min_overlap_bp:\s*([0-9]+)\s*(#.*)?$", m.group(1), re.M) if m else None
    if not k:
        raise ValueError("parameter sparcs.candidate_min_overlap_bp not found")
    v = int(k.group(1))
    if v < 1:
        raise ValueError("sparcs.candidate_min_overlap_bp must be >= 1 (a positive overlap)")
    return v


def verify_in_manifest(path: str) -> None:
    p = os.path.abspath(path)
    if "hervs_annotation4" in p or "ncRNA_layer" in p or "HERVs4_v2" in p:
        raise ValueError(f"{path}: historical tree; only the accepted release artifacts are inputs")
    d = os.path.dirname(p)
    while not os.path.exists(os.path.join(d, "checksums.sha256")):
        if d == os.path.dirname(d):
            raise ValueError(f"no checksums.sha256 above {path}")
        d = os.path.dirname(d)
    rel = os.path.relpath(p, d)
    want = None
    for line in open(os.path.join(d, "checksums.sha256")):
        h, f = line.rstrip("\n").split(None, 1)
        if f in (rel, "./" + rel):
            want = h
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if want is None or got != want:
        raise ValueError(f"{path} does not match its accepted manifest")
    print(f"[ok] {rel} sha256 {got} (manifest)", flush=True)


def load_bed_names(path: str, n_name_fields: int):
    """name -> (chrom, start, end, strand) of a gzipped BED6; names must be unique."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if n_name_fields:
                key = c[3].split("|")[0]
                if len(c[3].split("|")) != n_name_fields:
                    raise ValueError(f"{path}: unexpected name field {c[3]}")
            else:
                key = c[3]
            if key in out:
                raise ValueError(f"{path}: duplicate identifier {key}")
            out[key] = (c[0], int(c[1]), int(c[2]), c[5], c[3])
    return out


def guard_rows(df: pd.DataFrame, table: str, features: set, ltrs: dict, txs: dict, genes: dict) -> None:
    def refuse(mask, what):
        if mask.any():
            r = df.loc[mask].iloc[0]
            raise ValueError(f"{table}: {int(mask.sum())} row(s) {what}; first: {r['a_name']} x {r['b_name']}")
    coords = ["a_start", "a_end", "b_start", "b_end", "overlap_bp"]
    refuse(df[coords].isna().any(axis=1), "with a non-integer coordinate or overlap")
    refuse((df["a_start"] < 0) | (df["b_start"] < 0) | (df["a_end"] <= df["a_start"]) | (df["b_end"] <= df["b_start"]),
           "with an invalid half-open interval (start < 0 or end <= start)")
    refuse(~df["a_strand"].isin(["+", "-"]) | ~df["b_strand"].isin(["+", "-"]), "with a missing or invalid strand")
    refuse(df["a_strand"] == df["b_strand"], "on the same strand (the order-15 tables are -S, antisense only)")
    refuse(df["a_chrom"] != df["b_chrom"], "on different chromosomes")
    refuse((df["overlap_bp"] != df["ovl_len"]) | (df["ovl_len"] <= 0), "whose overlap_bp is not min(ends) - max(starts) > 0")
    refuse(df.duplicated(EXPECTED, keep=False), "duplicated")
    refuse(~df["feature"].isin(features), f"whose feature field is not {sorted(features)}")
    refuse(df["b_name"].str.count(r"\|") != 5, "whose b_name is not ENST|ENSG|GENE|FEATURE|GENE_TYPE|TX_TYPE")
    refuse(~df["tx_id"].str.match(TX_RE) | ~df["gene_id"].str.match(GENE_RE), "with a malformed transcript or gene id")
    a = list(zip(df["a_chrom"], df["a_start"].astype(int), df["a_end"].astype(int), df["a_strand"]))
    refuse(pd.Series([ltrs.get(n, (None,))[:4] != t for n, t in zip(df["a_name"], a)], index=df.index),
           "whose LTR is not an SRC07 record with these coordinates and strand")
    tx_ok = [(t in txs and txs[t][3] == s and txs[t][4] == f"{t}|{g}|{gn}|{tt}" and b[0] == txs[t][0]
              and txs[t][1] <= b[1] and b[2] <= txs[t][2])
             for t, g, gn, tt, s, b in zip(df["tx_id"], df["gene_id"], df["gene_name"], df["tx_type"], df["b_strand"],
                                           zip(df["b_chrom"], df["b_start"].astype(int), df["b_end"].astype(int)))]
    refuse(pd.Series([not x for x in tx_ok], index=df.index),
           "whose transcript is not a GENCODE v48 transcript with this gene, name, type, strand and span")
    gene_ok = [g in genes and genes[g][4] == f"{g}|{gn}|{gt}" and genes[g][3] == s
               for g, gn, gt, s in zip(df["gene_id"], df["gene_name"], df["gene_type"], df["b_strand"])]
    refuse(pd.Series([not x for x in gene_ok], index=df.index),
           "whose gene is not a GENCODE v48 gene with this name, type and strand")


def qc_rows(long_t: pd.DataFrame, cand: pd.DataFrame, pc_lnc: pd.DataFrame, br: pd.DataFrame, n_in: dict, n_ltr: int):
    rows = [("inputs", f"rows_read:{k}", v) for k, v in n_in.items()]
    rows.append(("inputs", "src07_ltrs", n_ltr))

    def unit(scope, t):
        rows.extend([
            (scope, "rows", len(t)),
            (scope, "unique_ltrs", t["a_name"].nunique()),
            (scope, "unique_ltr_x_transcript", len(t[["a_name", "tx_id"]].drop_duplicates())),
            (scope, "unique_ltr_x_gene", len(t[["a_name", "gene_id"]].drop_duplicates())),
            (scope, "unique_genes_gene_id", t["gene_id"].nunique()),
            (scope, "unique_genes_gene_name", t["gene_name"].nunique()),
            (scope, "unique_transcripts", t["tx_id"].nunique()),
            (scope, "ltrs_with_more_than_one_gene", int((t.groupby("a_name")["gene_id"].nunique() > 1).sum())),
            (scope, "ltrs_with_more_than_one_transcript", int((t.groupby("a_name")["tx_id"].nunique() > 1).sum())),
            (scope, "genes_with_more_than_one_ltr", int((t.groupby("gene_id")["a_name"].nunique() > 1).sum())),
        ])
    unit("long_table", long_t)
    unit("candidates", cand)
    unit("candidates_protein_coding_or_lncRNA", pc_lnc)
    for b in BRANCHES:
        unit(b, br[br["branch"] == b])
    for (gt, fc), n in cand.groupby(["gene_type", "feature_class"]).size().sort_index().items():
        rows.append(("candidates", f"rows:{gt}:{fc}", int(n)))
    for fc, n in long_t.groupby("feature_class").size().sort_index().items():
        rows.append(("long_table", f"rows:{fc}", int(n)))
    for b in BRANCHES:
        for tt, n in br[br["branch"] == b].groupby("tx_type").size().sort_index().items():
            rows.append((b, f"rows:tx_type:{tt}", int(n)))
    c, l = br[br["branch"] == BRANCHES[0]], br[br["branch"] == BRANCHES[1]]
    rows.extend([
        ("branch_intersection", "shared_ltrs", len(set(c["a_name"]) & set(l["a_name"]))),
        ("branch_intersection", "shared_genes_gene_id", len(set(c["gene_id"]) & set(l["gene_id"]))),
        ("branch_intersection", "shared_ltr_x_gene_events", len(set(zip(c["a_name"], c["gene_id"])) & set(zip(l["a_name"], l["gene_id"])))),
        ("branch_intersection", "rows_in_both_branches", len(set(zip(c["a_name"], c["tx_id"], c["feature_class"])) & set(zip(l["a_name"], l["tx_id"], l["feature_class"])))),
        ("branches", "union_unique_ltrs", len(set(c["a_name"]) | set(l["a_name"]))),
    ])
    return pd.DataFrame(rows, columns=["scope", "metric", "value"])


def main():
    ap = argparse.ArgumentParser(description="Antisense terminal-LTR candidates (07_sparcs step 1, effective order 43).")
    ap.add_argument("--overlap-dir", required=True, help="03 results/overlap_ltrs (the five order-15 tables)")
    ap.add_argument("--ltr-bed", required=True, help="03 results/herv_ltrs/ltrs.sorted.bed.gz (SRC07, sorted)")
    ap.add_argument("--transcripts-bed", required=True, help="03 gencode_features transcripts/sorted/gencode.v48.transcripts.sorted.bed.gz")
    ap.add_argument("--genes-bed", required=True, help="03 gencode_features genes/sorted/gencode.v48.genes.sorted.bed.gz")
    ap.add_argument("--params", required=True, help="config/parameters.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    min_overlap_bp = read_min_overlap(args.params)
    paths = {k: os.path.join(args.overlap_dir, f) for k, (f, _) in TABLES.items()}
    for p in [*paths.values(), args.ltr_bed, args.transcripts_bed, args.genes_bed]:
        verify_in_manifest(p)
    for p in (args.transcripts_bed, args.genes_bed):
        if ".v48." not in os.path.basename(p):
            raise ValueError(f"{p}: GENCODE release 48 assets required")
    ltrs = load_bed_names(args.ltr_bed, 0)
    txs = load_bed_names(args.transcripts_bed, 4)
    genes = load_bed_names(args.genes_bed, 3)
    print(f"[INFO] {len(ltrs):,} SRC07 LTRs; {len(txs):,} transcripts; {len(genes):,} genes", flush=True)

    frames, n_in = {}, {}
    for k, (f, features) in TABLES.items():
        df = read_ov(paths[k])
        guard_rows(df, f, features, ltrs, txs, genes)
        frames[k], n_in[k] = df, len(df)
        print(f"[ok] {f}: {len(df):,} rows pass the input guards", flush=True)

    ov_exon = frames["exon"];       ov_exon["feature_class"] = "exon"
    ov_cds  = frames["cds"];        ov_cds["feature_class"]  = "CDS"
    ov_utr  = frames["utr5p3p"];    ov_utr["feature_class"]  = ov_utr["feature"].fillna("UTR")
    ov_utr3 = frames["utr3"];       ov_utr3["feature_class"] = "three_prime_UTR"
    ov_last = frames["last_exon"];  ov_last["feature_class"] = "last_exon"

    ov_all = pd.concat([ov_exon, ov_cds, ov_utr, ov_utr3, ov_last], ignore_index=True)

    # Basic filter: real overlap
    ov_all = ov_all[ov_all["ovl_len"].fillna(0).astype(int) >= min_overlap_bp].copy()

    # Collapse within (LTR, transcript, feature_class) using union overlap bp
    group_cols = ["a_name", "a_chrom", "a_start", "a_end", "a_strand",
                  "tx_id", "gene_id", "gene_name", "gene_type", "tx_type",
                  "feature_class"]
    agg = (ov_all.groupby(group_cols, as_index=False)
           .agg(
               ltr_len_bp=("ltr_len_bp", "first"),
               overlap_bp_sum=("ovl_len", "sum"),  # placeholder
           ))

    union_bp = (ov_all.groupby(group_cols).apply(union_overlap_bp_from_group)).reset_index()
    union_bp = union_bp.rename(columns={0: "overlap_bp_union"})
    agg = agg.merge(union_bp, on=group_cols, how="left")
    agg["overlap_bp_total"] = pd.to_numeric(agg["overlap_bp_union"], errors="coerce").fillna(0).astype(int)
    agg.drop(columns=["overlap_bp_sum", "overlap_bp_union"], inplace=True)

    agg["overlap_frac_of_ltr"] = (agg["overlap_bp_total"] / agg["ltr_len_bp"]).fillna(0.0).clip(0, 1)

    # Candidate-focused: last exon OR 3'UTR (strong SPARCS-like contexts)
    cand = agg[agg["feature_class"].isin(["three_prime_UTR", "last_exon"])].copy()
    # Convenience: protein-coding + lncRNA only (but keep full table too)
    cand_pc_lnc = cand[cand["gene_type"].isin(["protein_coding", "lncRNA"])].copy()

    # RELEASE: the two downstream branches, exactly as E_sparcs/01_analysis_sparcs.R selects them
    coding = cand[(cand["gene_type"] == "protein_coding") & (cand["feature_class"] == "three_prime_UTR")].copy()
    lnc = cand[cand["gene_type"] == "lncRNA"].copy()
    if (lnc["feature_class"] != "last_exon").any():
        raise ValueError("a lncRNA candidate is not a last_exon row: the lncRNA_last_exon label would be ambiguous")
    coding.insert(0, "branch", BRANCHES[0])
    lnc.insert(0, "branch", BRANCHES[1])
    branches = pd.concat([coding, lnc], ignore_index=True)

    # output guards
    key = ["a_name", "tx_id", "feature_class"]
    for name, t in (("long", agg), ("candidates", cand), ("branches", branches)):
        k = key + (["branch"] if name == "branches" else [])
        if t.duplicated(k).any():
            raise ValueError(f"{name}: duplicate LTR x transcript x feature_class key")
        if ((t["a_end"] <= t["a_start"]) | (t["a_start"] < 0) | ~t["a_strand"].isin(["+", "-"])).any():
            raise ValueError(f"{name}: invalid LTR coordinates or strand")
        if ((t["overlap_bp_total"] < min_overlap_bp) | (t["overlap_bp_total"] > t["ltr_len_bp"])).any():
            raise ValueError(f"{name}: overlap_bp_total outside [min_overlap_bp, LTR length]")
        if not (set(t["a_name"]) <= set(ltrs) and set(t["tx_id"]) <= set(txs) and set(t["gene_id"]) <= set(genes)):
            raise ValueError(f"{name}: unknown LTR, transcript or gene identifier")
    if not set(branches["branch"]) <= set(BRANCHES):
        raise ValueError("impossible branch label")
    if len(branches) != len(coding) + len(lnc) or (set(coding["gene_type"]) & set(lnc["gene_type"])):
        raise ValueError("a row belongs to both branches")

    qc = qc_rows(agg, cand, cand_pc_lnc, branches, n_in, len(ltrs))
    qc.loc[len(qc)] = ("parameters", "candidate_min_overlap_bp", min_overlap_bp)

    # staged output: created only after every computation and guard has succeeded
    final_outdir = args.outdir
    parent = os.path.dirname(os.path.abspath(final_outdir)) or "."
    os.makedirs(parent, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".stage_candidates.", dir=parent)
    try:
        outputs = [
            ("ltr_tx_feature_overlaps.long.tsv", agg, ["feature_class", "overlap_bp_total", "overlap_frac_of_ltr"], [True, False, False]),
            ("ltr_candidates_antisense_last_exon_or_utr3.tsv", cand, ["feature_class", "overlap_bp_total"], [True, False]),
            ("ltr_candidates_antisense_last_exon_or_utr3.protein_coding_or_lncRNA.tsv", cand_pc_lnc,
             ["gene_type", "feature_class", "overlap_bp_total"], [True, True, False]),
            ("sparcs_candidate_branches.tsv", branches, ["branch"], [True]),
        ]
        for name, t, keys, asc in outputs:
            # historical keys first; the release tie-breakers make the order total
            t.sort_values(keys + TIE, ascending=asc + [True] * len(TIE), kind="mergesort") \
             .to_csv(os.path.join(stage, name), sep="\t", index=False)
        qc.to_csv(os.path.join(stage, "candidate_qc_summary.tsv"), sep="\t", index=False)
        with open(os.path.join(stage, "checksums.sha256"), "w") as fh:
            for name in sorted(n for n in os.listdir(stage) if n.endswith(".tsv")):
                fh.write(f"{hashlib.sha256(open(os.path.join(stage, name), 'rb').read()).hexdigest()}  {name}\n")
        for name in os.listdir(stage):
            os.chmod(os.path.join(stage, name), 0o644)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    if os.path.exists(final_outdir):
        shutil.rmtree(final_outdir)
    os.rename(stage, final_outdir)
    os.chmod(final_outdir, 0o755)

    q = {(s, m): v for s, m, v in qc.itertuples(index=False)}
    print(f"[OK] Long rows: {agg.shape[0]}")
    print(f"[OK] Candidate rows (utr3/last_exon): {cand.shape[0]}")
    print(f"[OK] Candidate rows (protein_coding|lncRNA): {cand_pc_lnc.shape[0]}")
    for b in BRANCHES:
        print(f"[OK] {b}: {q[(b, 'rows')]} rows, {q[(b, 'unique_ltrs')]} LTRs, {q[(b, 'unique_genes_gene_id')]} genes (gene_id)")
    print(f"[OK] wrote {final_outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
