#!/usr/bin/env python3
"""
04_domains_terminal_exons step 5 (frozen global order 33, effective order 34), CNAG:
terminal-exon architecture of transcripts and genes carrying conserved HERV-derived domains.

REFACTOR_SCRIPT of ncRNA_layer/scripts/integration/09_compare_terminal_exon_architecture.py
(sha256 b91fca85...). Scientific definitions are the historical ones, and the grouping,
summary, test, BH and gene-table functions are carried over verbatim:

  conserved association  order-30 all_coverage qualifying association (>= 80% exonic overlap)
                         with Coverage >= domains.min_hmm_profile_coverage
  terminal association   terminal_under_stringent_rule (>= 80% of the domain in the terminal exon)
  transcript groups      (priority order) terminal_conserved_domain > nonterminal_conserved_domain >
                         other_transcript_in_conserved_domain_gene > background_gene_without_conserved_domain
  unit                   MULTI-EXON transcripts (single-exon transcripts are excluded from the
                         summaries, tests and gene table; they stay in the transcript table)
  exon metrics           exons merged per transcript; terminal exon = 3'-most block;
                         length = end - start (BED); fraction = terminal / sum of merged blocks
  gene table             multi-exon transcripts per gene; group by the highest-ranked transcript
                         (terminal 2 > non-terminal 1 > other/background 0); medians and max
  tests                  scipy mannwhitneyu two-sided, method="asymptotic"; Cliff's delta =
                         2U/(n_a n_b) - 1 (> 0: group a larger); BH over each whole test table
                         (transcript and gene tables are separate families)

Refactor (declared before production):
  * the terminal exon comes from the shared 03 step-4 asset last_exon_per_tx (identical to the
    historical merged-block terminal for all 387,954 transcripts; a guard checks it per transcript);
  * coverage and exonic-overlap thresholds come from config/parameters.yaml (domains:);
  * inputs must match their accepted manifests; outputs are staged and published atomically;
  * no figures: plotting belongs to 08_figures_supplement (the historical PNGs were not used in
    the published figure), so the qc metric figures_created is 0 (historically 6).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import math
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


# ---------------------------------------------------------------- historical constants (verbatim)
TRANSCRIPT_GROUP_ORDER = [
    "terminal_conserved_domain",
    "nonterminal_conserved_domain",
    "other_transcript_in_conserved_domain_gene",
    "background_gene_without_conserved_domain",
]

PRIMARY_PLOT_GROUPS = [
    "terminal_conserved_domain",
    "nonterminal_conserved_domain",
    "background_gene_without_conserved_domain",
]

GROUP_LABELS = {
    "terminal_conserved_domain": "Terminal conserved\ndomain",
    "nonterminal_conserved_domain": "Non-terminal conserved\ndomain",
    "other_transcript_in_conserved_domain_gene": "Other transcript in\nhost gene",
    "background_gene_without_conserved_domain": "Background gene\n(no conserved domain)",
    "terminal_conserved_domain_gene": "Terminal conserved\ndomain gene",
    "nonterminal_only_conserved_domain_gene": "Non-terminal-only\nconserved domain gene",
    "background_gene_without_conserved_domain_gene": "Background gene\n(no conserved domain)",
}

TRANSCRIPT_METRICS = [
    "terminal_exon_length_bp",
    "total_exonic_length_bp",
    "terminal_exon_fraction",
    "n_exons",
]

GENE_METRICS = [
    "median_terminal_exon_length_bp",
    "max_terminal_exon_length_bp",
    "median_total_exonic_length_bp",
    "median_terminal_exon_fraction",
    "median_n_exons",
]



# ---------------------------------------------------------------- release: parameters and guards
def read_domain_params(path: str) -> Dict[str, float]:
    """domains: block of config/parameters.yaml (no YAML module on CNAG)."""
    text = open(path, encoding="utf-8").read()
    m = re.search(r"^domains:\s*\n(.*?)(?=^\S|\Z)", text, re.M | re.S)
    if not m:
        raise ValueError("no domains: block in parameters.yaml")
    out = {}
    for key in ("min_hmm_profile_coverage", "min_exonic_overlap_fraction"):
        k = re.search(r"^\s+" + re.escape(key) + r":\s*([0-9.eE+-]+)", m.group(1), re.M)
        if not k:
            raise ValueError(f"parameter domains.{key} not found")
        out[key] = float(k.group(1))
    return out


def verify_in_manifest(path: str) -> None:
    p = os.path.abspath(path)
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


def guard_inputs(args: argparse.Namespace) -> None:
    a = os.path.abspath(args.associations)
    if os.path.basename(a) != "domain_centered_qualifying_domain_transcript_associations.tsv" \
            or os.path.basename(os.path.dirname(a)) != "all_coverage":
        raise ValueError("--associations must be the order-30 all_coverage association table (D026)")
    for p in (args.tx_exons_bed, args.last_exons_bed):
        if "gene_exons_merged" in p or "gene_introns" in p:
            raise ValueError("transcript-resolved exons required (D022 gene-level assets are not an input)")
    for p in (args.associations, args.tx_exons_bed, args.last_exons_bed):
        verify_in_manifest(p)


# ---------------------------------------------------------------- historical helpers (verbatim)
def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    clean = sorted((int(start), int(end)) for start, end in intervals if int(end) > int(start))
    if not clean:
        return tuple()
    merged: List[List[int]] = [[clean[0][0], clean[0][1]]]
    for start, end in clean[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def parse_exon_name(name: str, line_number: int) -> Tuple[str, str, str, str, str, str]:
    parts = name.split("|")
    if len(parts) < 6:
        raise ValueError(
            f"Could not parse transcript metadata from exon BED line {line_number}: {name}"
        )
    return tuple(parts[:6])  # type: ignore[return-value]



# ---------------------------------------------------------------- transcript architecture
def load_terminal_exons(path: str) -> Dict[str, Tuple[int, int]]:
    """Shared 03 step-4 asset: one terminal exon per transcript."""
    terminal: Dict[str, Tuple[int, int]] = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            tx_id = fields[3].split("|", 1)[0]
            if tx_id in terminal:
                raise ValueError(f"last-exon asset lists {tx_id} twice")
            terminal[tx_id] = (int(fields[1]), int(fields[2]))
    return terminal


def build_transcript_architecture(exon_path: str, terminal_exons: Dict[str, Tuple[int, int]]) -> pd.DataFrame:
    """Historical build_transcript_architecture(), with the terminal exon taken from the shared
    step-4 asset; it must be the 3'-most merged block of the transcript (historical definition)."""
    raw: Dict[str, dict] = {}
    with open_text(exon_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Malformed exon BED line {line_number}: fewer than six columns")
            chrom, start, end, name, _score, strand = fields[:6]
            tx_id, gene_id, gene_name, _feature, gene_type, tx_type = parse_exon_name(name, line_number)
            if tx_id not in raw:
                raw[tx_id] = {"tx_id": tx_id, "gene_id": gene_id, "gene_name": gene_name, "gene_type": gene_type,
                              "tx_type": tx_type, "chrom": chrom, "strand": strand, "exons": []}
            else:
                previous = raw[tx_id]
                expected = (previous["gene_id"], previous["gene_name"], previous["gene_type"],
                            previous["tx_type"], previous["chrom"], previous["strand"])
                if expected != (gene_id, gene_name, gene_type, tx_type, chrom, strand):
                    raise ValueError(f"Inconsistent annotation for transcript {tx_id}")
            raw[tx_id]["exons"].append((int(start), int(end)))

    rows = []
    for tx_id, record in raw.items():
        exons = merge_intervals(record["exons"])
        if not exons:
            continue
        terminal_exon = terminal_exons.get(tx_id)
        if terminal_exon != (exons[-1] if record["strand"] == "+" else exons[0]):
            raise ValueError(f"step-4 terminal exon of {tx_id} is not its 3'-most merged exon block")
        terminal_length = terminal_exon[1] - terminal_exon[0]
        total_exonic_length = sum(end - start for start, end in exons)
        n_exons = len(exons)
        rows.append({
            "tx_id": tx_id, "gene_id": record["gene_id"], "gene_name": record["gene_name"],
            "gene_type": record["gene_type"], "tx_type": record["tx_type"], "chrom": record["chrom"],
            "strand": record["strand"], "transcript_start": exons[0][0], "transcript_end": exons[-1][1],
            "n_exons": n_exons, "is_single_exon_transcript": n_exons == 1,
            "terminal_exon_start": terminal_exon[0], "terminal_exon_end": terminal_exon[1],
            "terminal_exon_length_bp": terminal_length, "total_exonic_length_bp": total_exonic_length,
            "terminal_exon_fraction": terminal_length / total_exonic_length,
        })
    transcript_table = pd.DataFrame(rows)
    if transcript_table.empty:
        raise RuntimeError("No transcript models were loaded")
    if transcript_table["tx_id"].duplicated().any():
        raise ValueError("Transcript architecture table contains duplicate transcript IDs")
    return transcript_table


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    normalized = series.astype(str).str.strip().str.lower()
    return normalized.isin({"true", "t", "1", "yes"})



def load_conserved_associations(path: str, coverage_threshold: float, min_exonic_overlap: float) -> pd.DataFrame:
    """Historical load_conserved_associations(); the qualifying-overlap guard reads its threshold
    from parameters.yaml instead of a literal 0.80."""
    associations = pd.read_csv(path, sep="\t", dtype={"tx_id": str, "gene_id": str, "domain_uid": str})
    required = {"domain_uid", "tx_id", "gene_id", "Coverage", "terminal_under_stringent_rule", "exonic_overlap_fraction"}
    missing = sorted(required - set(associations.columns))
    if missing:
        raise ValueError("Association table is missing columns: " + ", ".join(missing))
    associations["Coverage"] = pd.to_numeric(associations["Coverage"], errors="raise")
    associations["exonic_overlap_fraction"] = pd.to_numeric(associations["exonic_overlap_fraction"], errors="raise")
    associations["terminal_under_stringent_rule"] = as_bool(associations["terminal_under_stringent_rule"])
    if (associations["exonic_overlap_fraction"] < min_exonic_overlap - 1e-12).any():
        raise ValueError("The input contains associations below the qualifying exonic overlap; use the order-30 table.")
    conserved = associations.loc[associations["Coverage"] >= coverage_threshold].copy()
    if conserved.empty:
        raise RuntimeError("No conserved associations remain at the selected threshold")
    return conserved


# ---------------------------------------------------------------- historical grouping, summaries, tests,
# ---------------------------------------------------------------- BH, gene table, counts (verbatim)
def annotate_transcript_groups(
    transcript_table: pd.DataFrame,
    conserved: pd.DataFrame,
) -> pd.DataFrame:
    per_tx = (
        conserved.groupby("tx_id", sort=False)
        .agg(
            n_conserved_domain_associations=("domain_uid", "size"),
            n_unique_conserved_domains=("domain_uid", "nunique"),
            n_terminal_conserved_domains=(
                "terminal_under_stringent_rule",
                lambda x: int(np.sum(x.astype(bool))),
            ),
        )
        .reset_index()
    )
    per_tx["n_nonterminal_conserved_domains"] = (
        per_tx["n_conserved_domain_associations"] - per_tx["n_terminal_conserved_domains"]
    )

    host_genes = set(conserved["gene_id"].astype(str))
    terminal_tx = set(
        conserved.loc[conserved["terminal_under_stringent_rule"], "tx_id"].astype(str)
    )
    conserved_tx = set(conserved["tx_id"].astype(str))

    table = transcript_table.merge(per_tx, on="tx_id", how="left", validate="one_to_one")
    count_columns = [
        "n_conserved_domain_associations",
        "n_unique_conserved_domains",
        "n_terminal_conserved_domains",
        "n_nonterminal_conserved_domains",
    ]
    table[count_columns] = table[count_columns].fillna(0).astype(int)
    table["gene_has_qualifying_conserved_domain"] = table["gene_id"].isin(host_genes)
    table["transcript_has_qualifying_conserved_domain"] = table["tx_id"].isin(conserved_tx)
    table["transcript_has_terminal_conserved_domain"] = table["tx_id"].isin(terminal_tx)

    conditions = [
        table["transcript_has_terminal_conserved_domain"],
        table["transcript_has_qualifying_conserved_domain"],
        table["gene_has_qualifying_conserved_domain"],
    ]
    choices = [
        "terminal_conserved_domain",
        "nonterminal_conserved_domain",
        "other_transcript_in_conserved_domain_gene",
    ]
    table["transcript_domain_group"] = np.select(
        conditions,
        choices,
        default="background_gene_without_conserved_domain",
    )
    return table


def analysis_scopes(table: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    scopes = {"all_gene_types": table}
    for gene_type in ["lncRNA", "protein_coding"]:
        subset = table.loc[table["gene_type"] == gene_type].copy()
        if not subset.empty:
            scopes[gene_type] = subset
    return scopes


def summarize_groups(
    table: pd.DataFrame,
    group_column: str,
    metrics: Sequence[str],
    level: str,
) -> pd.DataFrame:
    rows = []
    for scope, scoped in analysis_scopes(table).items():
        for group, grouped in scoped.groupby(group_column, sort=False):
            for metric in metrics:
                values = pd.to_numeric(grouped[metric], errors="coerce").dropna().to_numpy(float)
                if values.size == 0:
                    continue
                rows.append(
                    {
                        "level": level,
                        "scope": scope,
                        "group": group,
                        "metric": metric,
                        "n": values.size,
                        "mean": float(np.mean(values)),
                        "sd": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
                        "median": float(np.median(values)),
                        "q1": float(np.quantile(values, 0.25)),
                        "q3": float(np.quantile(values, 0.75)),
                        "minimum": float(np.min(values)),
                        "maximum": float(np.max(values)),
                    }
                )
    return pd.DataFrame(rows)


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not np.any(valid):
        return adjusted
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = ranked.size
    corrected = ranked * n / np.arange(1, n + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    corrected = np.clip(corrected, 0, 1)
    restored = np.empty(n, dtype=float)
    restored[order] = corrected
    adjusted[valid] = restored
    return adjusted


def pairwise_tests(
    table: pd.DataFrame,
    group_column: str,
    metrics: Sequence[str],
    comparisons: Sequence[Tuple[str, str]],
    level: str,
    min_group_size: int,
) -> pd.DataFrame:
    rows = []
    for scope, scoped in analysis_scopes(table).items():
        for group_a, group_b in comparisons:
            for metric in metrics:
                x = pd.to_numeric(
                    scoped.loc[scoped[group_column] == group_a, metric], errors="coerce"
                ).dropna().to_numpy(float)
                y = pd.to_numeric(
                    scoped.loc[scoped[group_column] == group_b, metric], errors="coerce"
                ).dropna().to_numpy(float)

                if x.size < min_group_size or y.size < min_group_size:
                    rows.append(
                        {
                            "level": level,
                            "scope": scope,
                            "metric": metric,
                            "group_a": group_a,
                            "group_b": group_b,
                            "n_a": x.size,
                            "n_b": y.size,
                            "median_a": float(np.median(x)) if x.size else np.nan,
                            "median_b": float(np.median(y)) if y.size else np.nan,
                            "median_difference_a_minus_b": (
                                float(np.median(x) - np.median(y))
                                if x.size and y.size
                                else np.nan
                            ),
                            "median_ratio_a_over_b": (
                                float(np.median(x) / np.median(y))
                                if x.size and y.size and np.median(y) != 0
                                else np.nan
                            ),
                            "mann_whitney_u": np.nan,
                            "cliffs_delta": np.nan,
                            "p_two_sided": np.nan,
                            "test_status": "insufficient_group_size",
                        }
                    )
                    continue

                result = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
                u_stat = float(result.statistic)
                cliffs_delta = 2.0 * u_stat / (x.size * y.size) - 1.0
                rows.append(
                    {
                        "level": level,
                        "scope": scope,
                        "metric": metric,
                        "group_a": group_a,
                        "group_b": group_b,
                        "n_a": x.size,
                        "n_b": y.size,
                        "median_a": float(np.median(x)),
                        "median_b": float(np.median(y)),
                        "median_difference_a_minus_b": float(np.median(x) - np.median(y)),
                        "median_ratio_a_over_b": (
                            float(np.median(x) / np.median(y))
                            if np.median(y) != 0
                            else np.nan
                        ),
                        "mann_whitney_u": u_stat,
                        "cliffs_delta": float(cliffs_delta),
                        "p_two_sided": float(result.pvalue),
                        "test_status": "tested",
                    }
                )

    output = pd.DataFrame(rows)
    if not output.empty:
        output["p_adjust_bh"] = benjamini_hochberg(output["p_two_sided"])
    return output


def build_gene_table(transcript_table: pd.DataFrame) -> pd.DataFrame:
    multi = transcript_table.loc[~transcript_table["is_single_exon_transcript"]].copy()
    multi["is_terminal_domain_transcript"] = (
        multi["transcript_domain_group"] == "terminal_conserved_domain"
    ).astype(int)
    multi["is_nonterminal_domain_transcript"] = (
        multi["transcript_domain_group"] == "nonterminal_conserved_domain"
    ).astype(int)
    group_rank = {
        "background_gene_without_conserved_domain": 0,
        "other_transcript_in_conserved_domain_gene": 0,
        "nonterminal_conserved_domain": 1,
        "terminal_conserved_domain": 2,
    }
    multi["gene_group_rank"] = multi["transcript_domain_group"].map(group_rank).astype(int)

    gene_table = (
        multi.groupby("gene_id", sort=False, observed=True)
        .agg(
            gene_name=("gene_name", "first"),
            gene_type=("gene_type", "first"),
            chrom=("chrom", "first"),
            strand=("strand", "first"),
            gene_group_rank=("gene_group_rank", "max"),
            n_multi_exon_transcripts=("tx_id", "nunique"),
            n_terminal_domain_transcripts=("is_terminal_domain_transcript", "sum"),
            n_nonterminal_domain_transcripts=("is_nonterminal_domain_transcript", "sum"),
            median_terminal_exon_length_bp=("terminal_exon_length_bp", "median"),
            max_terminal_exon_length_bp=("terminal_exon_length_bp", "max"),
            median_total_exonic_length_bp=("total_exonic_length_bp", "median"),
            median_terminal_exon_fraction=("terminal_exon_fraction", "median"),
            median_n_exons=("n_exons", "median"),
        )
        .reset_index()
    )

    gene_table["gene_domain_group"] = gene_table["gene_group_rank"].map(
        {
            2: "terminal_conserved_domain_gene",
            1: "nonterminal_only_conserved_domain_gene",
            0: "background_gene_without_conserved_domain_gene",
        }
    )
    gene_table = gene_table.drop(columns=["gene_group_rank"])
    return gene_table


def make_group_counts(
    table: pd.DataFrame,
    group_column: str,
    level: str,
) -> pd.DataFrame:
    rows = []
    for scope, scoped in analysis_scopes(table).items():
        counts = scoped[group_column].value_counts()
        total = len(scoped)
        for group, n in counts.items():
            rows.append(
                {
                    "level": level,
                    "scope": scope,
                    "group": group,
                    "n": int(n),
                    "percent_of_scope": 100.0 * n / total if total else np.nan,
                }
            )
    return pd.DataFrame(rows)



# ---------------------------------------------------------------- main
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Terminal-exon architecture (effective order 34).")
    parser.add_argument("--associations", required=True, help="order-30 all_coverage association table")
    parser.add_argument("--tx-exons-bed", required=True, help="03 step-3 sorted transcript exon BED")
    parser.add_argument("--last-exons-bed", required=True, help="03 step-4 last_exon_per_tx BED")
    parser.add_argument("--params", required=True, help="config/parameters.yaml")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-group-size", type=int, default=5,
                        help="Minimum observations per group required for a statistical comparison (historical default)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_group_size < 2:
        raise ValueError("--min-group-size must be at least 2")
    params = read_domain_params(args.params)
    coverage_threshold = params["min_hmm_profile_coverage"]
    guard_inputs(args)

    transcripts = build_transcript_architecture(args.tx_exons_bed, load_terminal_exons(args.last_exons_bed))
    print(f"[INFO] Loaded {len(transcripts):,} transcript models", flush=True)
    conserved = load_conserved_associations(args.associations, coverage_threshold, params["min_exonic_overlap_fraction"])
    print(f"[INFO] {len(conserved):,} conserved domain-transcript associations at coverage >= {coverage_threshold:g}", flush=True)
    missing_transcripts = sorted(set(conserved["tx_id"]) - set(transcripts["tx_id"]))
    if missing_transcripts:
        raise ValueError("Conserved associations reference transcripts absent from the exon BED: " + ", ".join(missing_transcripts[:10]))

    transcript_table = annotate_transcript_groups(transcripts, conserved)
    multi = transcript_table.loc[~transcript_table["is_single_exon_transcript"]].copy()
    transcript_group_counts = make_group_counts(multi, "transcript_domain_group", "transcript")
    transcript_summary = summarize_groups(multi, "transcript_domain_group", TRANSCRIPT_METRICS, "transcript")
    transcript_comparisons = [
        ("terminal_conserved_domain", "background_gene_without_conserved_domain"),
        ("nonterminal_conserved_domain", "background_gene_without_conserved_domain"),
        ("terminal_conserved_domain", "nonterminal_conserved_domain"),
        ("terminal_conserved_domain", "other_transcript_in_conserved_domain_gene"),
    ]
    transcript_tests = pairwise_tests(multi, "transcript_domain_group", TRANSCRIPT_METRICS, transcript_comparisons,
                                      "transcript", args.min_group_size)
    gene_table = build_gene_table(transcript_table)
    gene_group_counts = make_group_counts(gene_table, "gene_domain_group", "gene")
    gene_summary = summarize_groups(gene_table, "gene_domain_group", GENE_METRICS, "gene")
    gene_comparisons = [
        ("terminal_conserved_domain_gene", "background_gene_without_conserved_domain_gene"),
        ("nonterminal_only_conserved_domain_gene", "background_gene_without_conserved_domain_gene"),
        ("terminal_conserved_domain_gene", "nonterminal_only_conserved_domain_gene"),
    ]
    gene_tests = pairwise_tests(gene_table, "gene_domain_group", GENE_METRICS, gene_comparisons, "gene", args.min_group_size)

    # staged output: created only after every computation and guard has succeeded
    final_outdir = args.outdir
    parent = os.path.dirname(os.path.abspath(final_outdir)) or "."
    os.makedirs(parent, exist_ok=True)
    outdir = tempfile.mkdtemp(prefix=".stage_architecture.", dir=parent)
    tables = {
        "transcript_architecture_groups.tsv": transcript_table,
        "transcript_architecture_group_counts.tsv": transcript_group_counts,
        "transcript_architecture_summary.tsv": transcript_summary,
        "transcript_architecture_pairwise_tests.tsv": transcript_tests,
        "gene_architecture_groups.tsv": gene_table,
        "gene_architecture_group_counts.tsv": gene_group_counts,
        "gene_architecture_summary.tsv": gene_summary,
        "gene_architecture_pairwise_tests.tsv": gene_tests,
    }
    for name, table in tables.items():
        table.to_csv(os.path.join(outdir, name), sep="\t", index=False)
    qc = pd.DataFrame({
        "metric": ["coverage_threshold", "all_gencode_transcripts", "all_gencode_genes", "multi_exon_transcripts",
                   "multi_exon_genes", "conserved_domain_transcript_associations", "unique_conserved_domains_in_associations",
                   "unique_conserved_domain_host_genes", "terminal_conserved_domain_transcripts_multi_exon",
                   "nonterminal_conserved_domain_transcripts_multi_exon", "other_transcripts_in_conserved_domain_genes_multi_exon",
                   "background_transcripts_multi_exon", "terminal_conserved_domain_genes_multi_exon",
                   "nonterminal_only_conserved_domain_genes_multi_exon", "background_genes_multi_exon", "figures_created"],
        "value": [coverage_threshold, len(transcript_table), transcript_table["gene_id"].nunique(), len(multi),
                  multi["gene_id"].nunique(), len(conserved), conserved["domain_uid"].nunique(), conserved["gene_id"].nunique(),
                  int((multi["transcript_domain_group"] == "terminal_conserved_domain").sum()),
                  int((multi["transcript_domain_group"] == "nonterminal_conserved_domain").sum()),
                  int((multi["transcript_domain_group"] == "other_transcript_in_conserved_domain_gene").sum()),
                  int((multi["transcript_domain_group"] == "background_gene_without_conserved_domain").sum()),
                  int((gene_table["gene_domain_group"] == "terminal_conserved_domain_gene").sum()),
                  int((gene_table["gene_domain_group"] == "nonterminal_only_conserved_domain_gene").sum()),
                  int((gene_table["gene_domain_group"] == "background_gene_without_conserved_domain_gene").sum()),
                  0],   # REFACTOR: no figures rendered here (08_figures_supplement)
    })
    qc.to_csv(os.path.join(outdir, "terminal_exon_architecture_qc.tsv"), sep="\t", index=False)

    for name in os.listdir(outdir):
        os.chmod(os.path.join(outdir, name), 0o644)
    if os.path.exists(final_outdir):
        shutil.rmtree(final_outdir)
    os.rename(outdir, final_outdir)
    os.chmod(final_outdir, 0o755)

    central = transcript_tests.loc[(transcript_tests["metric"] == "terminal_exon_length_bp")
                                   & (transcript_tests["group_a"] == "terminal_conserved_domain")
                                   & (transcript_tests["group_b"] == "background_gene_without_conserved_domain")]
    print("\n=== PRIMARY TERMINAL-EXON LENGTH COMPARISON ===")
    print(central[["scope", "n_a", "n_b", "median_a", "median_b", "cliffs_delta", "p_two_sided", "p_adjust_bh"]].to_string(index=False))
    print(f"[OK] wrote {final_outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
