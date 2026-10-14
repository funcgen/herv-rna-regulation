#!/usr/bin/env python3
# ==============================================================================
# 04_domains_terminal_exons step 2 (frozen global order 30, effective order 31), CNAG
# Transcript-aware positional null. COPY_RENAME of the historical
#   ncRNA_layer/scripts/integration/06_domain_centered_terminal_exon_null_v2.py
#   (sha256 d65f044c...). v1 (06_domain_centered_terminal_exon_null.py) is superseded
#   and NOT migrated.
# The null model, RNG, statistics, CLI and defaults are the historical v2 code,
# unchanged. Release-only additions, all marked "RELEASE:":
#   * input guards: the three inputs must match their accepted manifests, and the
#     domain / association tables must be the order-30 conserved profile;
#   * outputs are written to a staging directory and published atomically;
#   * the permutation table's gzip header carries mtime 0 (pandas otherwise writes the
#     run time), so re-runs are byte-identical; its decompressed content is unchanged.
# ==============================================================================
"""
Domain-centred, transcript-architecture-aware terminal-exon null model (v2).

Biological question
-------------------
Among unique conserved HERV-derived domains that overlap exonic sequence
stringently (>= min_exon_overlap of the domain), is terminal-exon overlap more
frequent than expected from domain length, exon length, and transcript
architecture?

Primary null: observed_qualifying_transcripts
---------------------------------------------
For each domain, preserve the exact transcript IDs in which that domain met the
observed >= min_exon_overlap criterion. Enumerate every genomic start position
within those transcript models at which the same-length domain would still meet
the exonic-overlap criterion. This is the closest match to the original
selection rule and directly controls for transcript-specific exon lengths.

Sensitivity null: all_isoforms_of_associated_genes
--------------------------------------------------
For each domain, use all same-chromosome, same-strand transcript isoforms of the
associated observed gene(s). This is a broader gene-matched null and reproduces
the assumption used by the first version of this script.

For every domain and scope, the script calculates exact position-based null
probabilities and then reports both:
  * exact Poisson-binomial tail probabilities for the total number of terminal
    domains, assuming domains are the analysis units; and
  * Monte Carlo empirical probabilities from n_permutations simulated datasets.

Notes
-----
* The expected count is exact under the stated domain-level null.
* Multiple domains from the same HERV locus are treated as separate analysis
  units, matching the requested domain-centred analysis. This independence
  assumption affects null variance, not the expected count.
* "terminal in any transcript" is permissive in isoform-rich genes. Therefore,
  "terminal in all qualifying transcripts", a stringent terminal-overlap rule,
  and multi-exon-only analyses are reported as sensitivities.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import math
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TranscriptModel:
    tx_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: Tuple[Tuple[int, int], ...]
    terminal_exon: Tuple[int, int]

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Domain-centred terminal-exon null model with transcript-matched "
            "primary and gene-isoform sensitivity scopes."
        )
    )
    parser.add_argument("--domains", required=True)
    parser.add_argument("--associations", required=True)
    parser.add_argument("--tx-exons-bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-exon-overlap", type=float, default=0.80)
    parser.add_argument("--n-permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument(
        "--scope",
        choices=("both", "observed", "gene_all"),
        default="both",
        help=(
            "Null transcript scope. 'observed' preserves observed qualifying "
            "transcript IDs; 'gene_all' uses all isoforms of associated genes."
        ),
    )
    return parser.parse_args()


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    clean = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    if not clean:
        return tuple()
    merged: List[List[int]] = [[clean[0][0], clean[0][1]]]
    for start, end in clean[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((s, e) for s, e in merged)


def load_transcripts(
    exon_path: str,
    relevant_genes: set[str],
) -> Tuple[Dict[str, TranscriptModel], Dict[str, List[str]]]:
    raw: Dict[str, dict] = {}

    with open_text(exon_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(
                    f"Malformed exon BED line {line_number}: fewer than 6 columns"
                )

            chrom, start, end, name, _score, strand = fields[:6]
            parts = name.split("|")
            if len(parts) < 2:
                raise ValueError(
                    "Could not parse transcript/gene identifiers from exon BED "
                    f"line {line_number}"
                )
            tx_id, gene_id = parts[0], parts[1]
            if gene_id not in relevant_genes:
                continue

            start_i, end_i = int(start), int(end)
            if tx_id not in raw:
                raw[tx_id] = {
                    "gene_id": gene_id,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                }
            else:
                previous = raw[tx_id]
                if (
                    previous["gene_id"] != gene_id
                    or previous["chrom"] != chrom
                    or previous["strand"] != strand
                ):
                    raise ValueError(f"Inconsistent annotation for transcript {tx_id}")
            raw[tx_id]["exons"].append((start_i, end_i))

    transcripts: Dict[str, TranscriptModel] = {}
    gene_to_transcripts: Dict[str, List[str]] = defaultdict(list)

    for tx_id, item in raw.items():
        exons = merge_intervals(item["exons"])
        if not exons:
            continue
        terminal = exons[-1] if item["strand"] == "+" else exons[0]
        model = TranscriptModel(
            tx_id=tx_id,
            gene_id=item["gene_id"],
            chrom=item["chrom"],
            strand=item["strand"],
            exons=exons,
            terminal_exon=terminal,
        )
        transcripts[tx_id] = model
        gene_to_transcripts[model.gene_id].append(tx_id)

    missing = sorted(relevant_genes - set(gene_to_transcripts))
    if missing:
        raise ValueError(
            "No transcript exons found for relevant genes: " + ", ".join(missing[:10])
        )

    return transcripts, gene_to_transcripts


def transcript_window_overlap(
    model: TranscriptModel,
    search_start: int,
    search_end_exclusive: int,
    domain_length: int,
    n_starts: int,
) -> np.ndarray:
    """Exonic bp overlap for every possible domain start in one transcript."""
    sequence_length = search_end_exclusive - search_start
    diff = np.zeros(sequence_length + 1, dtype=np.int16)

    for exon_start, exon_end in model.exons:
        clipped_start = max(exon_start, search_start)
        clipped_end = min(exon_end, search_end_exclusive)
        if clipped_end <= clipped_start:
            continue
        diff[clipped_start - search_start] += 1
        diff[clipped_end - search_start] -= 1

    covered = (np.cumsum(diff[:-1]) > 0).astype(np.int8)
    cumulative = np.concatenate(([0], np.cumsum(covered, dtype=np.int32)))
    return cumulative[domain_length : domain_length + n_starts] - cumulative[:n_starts]


def terminal_overlap_bp(
    starts: np.ndarray,
    domain_length: int,
    terminal_exon: Tuple[int, int],
) -> np.ndarray:
    terminal_start, terminal_end = terminal_exon
    return np.maximum(
        0,
        np.minimum(starts + domain_length, terminal_end)
        - np.maximum(starts, terminal_start),
    )


def filter_compatible_transcripts(
    tx_ids: Sequence[str],
    transcripts: Dict[str, TranscriptModel],
    chrom: str,
    strand: str,
) -> List[str]:
    return sorted(
        {
            tx_id
            for tx_id in tx_ids
            if tx_id in transcripts
            and transcripts[tx_id].chrom == chrom
            and transcripts[tx_id].strand == strand
        }
    )


def enumerate_domain_null(
    domain_row: pd.Series,
    tx_ids: Sequence[str],
    transcripts: Dict[str, TranscriptModel],
    min_exon_overlap: float,
    scope: str,
) -> dict:
    domain_length = int(domain_row["domain_length_bp"])
    minimum_bp = int(math.ceil(min_exon_overlap * domain_length))
    chrom = str(domain_row["chrom"])
    strand = str(domain_row["strand"])

    tx_ids = filter_compatible_transcripts(tx_ids, transcripts, chrom, strand)
    if not tx_ids:
        raise ValueError(
            f"No compatible transcript models for {domain_row['domain_uid']} in {scope}"
        )

    # A domain may extend outside an exon/transcript by at most L-minimum_bp while
    # still satisfying the >= min_exon_overlap criterion.
    search_start = min(transcripts[tx_id].start for tx_id in tx_ids) - (
        domain_length - minimum_bp
    )
    final_start = max(transcripts[tx_id].end for tx_id in tx_ids) - minimum_bp
    if final_start < search_start:
        raise ValueError(f"Empty candidate range for {domain_row['domain_uid']}")

    starts = np.arange(search_start, final_start + 1, dtype=np.int64)
    n_starts = starts.size
    search_end_exclusive = final_start + domain_length

    qualifying_count = np.zeros(n_starts, dtype=np.uint16)
    terminal_count = np.zeros(n_starts, dtype=np.uint16)
    stringent_terminal_count = np.zeros(n_starts, dtype=np.uint16)

    multi_qualifying_count = np.zeros(n_starts, dtype=np.uint16)
    multi_terminal_count = np.zeros(n_starts, dtype=np.uint16)
    multi_stringent_terminal_count = np.zeros(n_starts, dtype=np.uint16)

    for tx_id in tx_ids:
        model = transcripts[tx_id]
        exonic_bp = transcript_window_overlap(
            model=model,
            search_start=search_start,
            search_end_exclusive=search_end_exclusive,
            domain_length=domain_length,
            n_starts=n_starts,
        )
        qualifies = exonic_bp >= minimum_bp
        if not np.any(qualifies):
            continue

        terminal_bp = terminal_overlap_bp(starts, domain_length, model.terminal_exon)
        is_terminal = terminal_bp > 0
        is_stringently_terminal = terminal_bp >= minimum_bp

        qualifying_count += qualifies
        terminal_count += qualifies & is_terminal
        stringent_terminal_count += qualifies & is_stringently_terminal

        if model.n_exons > 1:
            multi_qualifying_count += qualifies
            multi_terminal_count += qualifies & is_terminal
            multi_stringent_terminal_count += qualifies & is_stringently_terminal

    valid = qualifying_count > 0
    valid_multi = multi_qualifying_count > 0
    if not np.any(valid):
        raise ValueError(f"No valid randomized positions for {domain_row['domain_uid']}")

    terminal_any = valid & (terminal_count > 0)
    terminal_all = valid & (terminal_count == qualifying_count)
    stringent_terminal_any = valid & (stringent_terminal_count > 0)
    stringent_terminal_all = valid & (stringent_terminal_count == qualifying_count)

    terminal_any_multi = valid_multi & (multi_terminal_count > 0)
    terminal_all_multi = valid_multi & (multi_terminal_count == multi_qualifying_count)
    stringent_terminal_any_multi = valid_multi & (multi_stringent_terminal_count > 0)

    n_valid = int(valid.sum())
    n_valid_multi = int(valid_multi.sum())

    def safe_fraction(num: int, den: int) -> float:
        return float(num / den) if den > 0 else np.nan

    return {
        "scope": scope,
        "domain_uid": domain_row["domain_uid"],
        "n_transcript_models": len(tx_ids),
        "n_single_exon_models": sum(transcripts[t].n_exons == 1 for t in tx_ids),
        "n_multi_exon_models": sum(transcripts[t].n_exons > 1 for t in tx_ids),
        "domain_length_bp": domain_length,
        "minimum_exonic_bp": minimum_bp,
        "candidate_start_min": search_start,
        "candidate_start_max": final_start,
        "n_candidate_starts": n_starts,
        "n_valid_exonic_starts": n_valid,
        "n_terminal_any_starts": int(terminal_any.sum()),
        "n_terminal_all_starts": int(terminal_all.sum()),
        "n_stringent_terminal_any_starts": int(stringent_terminal_any.sum()),
        "n_stringent_terminal_all_starts": int(stringent_terminal_all.sum()),
        "p_terminal_any": safe_fraction(int(terminal_any.sum()), n_valid),
        "p_terminal_all": safe_fraction(int(terminal_all.sum()), n_valid),
        "p_stringent_terminal_any": safe_fraction(
            int(stringent_terminal_any.sum()), n_valid
        ),
        "p_stringent_terminal_all": safe_fraction(
            int(stringent_terminal_all.sum()), n_valid
        ),
        "n_valid_multi_exon_starts": n_valid_multi,
        "n_terminal_any_multi_exon_starts": int(terminal_any_multi.sum()),
        "n_terminal_all_multi_exon_starts": int(terminal_all_multi.sum()),
        "n_stringent_terminal_any_multi_exon_starts": int(
            stringent_terminal_any_multi.sum()
        ),
        "p_terminal_any_multi_exon": safe_fraction(
            int(terminal_any_multi.sum()), n_valid_multi
        ),
        "p_terminal_all_multi_exon": safe_fraction(
            int(terminal_all_multi.sum()), n_valid_multi
        ),
        "p_stringent_terminal_any_multi_exon": safe_fraction(
            int(stringent_terminal_any_multi.sum()), n_valid_multi
        ),
    }


def simulate_counts(
    probabilities: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
    chunk_size: int = 1000,
) -> np.ndarray:
    probabilities = np.asarray(probabilities, dtype=float)
    if probabilities.size == 0 or np.isnan(probabilities).any():
        raise ValueError("Invalid probability vector for simulation")

    output = np.empty(n_permutations, dtype=np.int32)
    offset = 0
    while offset < n_permutations:
        n_chunk = min(chunk_size, n_permutations - offset)
        draws = rng.random((n_chunk, probabilities.size)) < probabilities
        output[offset : offset + n_chunk] = draws.sum(axis=1)
        offset += n_chunk
    return output


def poisson_binomial_pmf(probabilities: np.ndarray) -> np.ndarray:
    """Exact PMF of a sum of independent Bernoulli variables."""
    pmf = np.array([1.0], dtype=float)
    for probability in np.asarray(probabilities, dtype=float):
        pmf = np.convolve(pmf, np.array([1.0 - probability, probability]))
    # Guard against tiny floating-point drift.
    pmf = np.clip(pmf, 0.0, 1.0)
    pmf /= pmf.sum()
    return pmf


def discrete_interval_from_pmf(pmf: np.ndarray, lower_q=0.025, upper_q=0.975):
    cdf = np.cumsum(pmf)
    lower = int(np.searchsorted(cdf, lower_q, side="left"))
    upper = int(np.searchsorted(cdf, upper_q, side="left"))
    return lower, upper


def summarize_null(
    scope: str,
    scope_role: str,
    analysis: str,
    observed: int,
    probabilities: np.ndarray,
    simulated: np.ndarray,
) -> dict:
    probabilities = np.asarray(probabilities, dtype=float)
    expected = float(probabilities.sum())

    exact_pmf = poisson_binomial_pmf(probabilities)
    exact_lower, exact_upper = discrete_interval_from_pmf(exact_pmf)
    exact_p_enrichment = float(exact_pmf[observed:].sum())
    exact_p_depletion = float(exact_pmf[: observed + 1].sum())
    exact_p_two_sided = min(1.0, 2.0 * min(exact_p_enrichment, exact_p_depletion))

    mc_lower, mc_upper = np.quantile(simulated, [0.025, 0.975])
    mc_p_enrichment = float((1 + np.sum(simulated >= observed)) / (len(simulated) + 1))
    mc_p_depletion = float((1 + np.sum(simulated <= observed)) / (len(simulated) + 1))
    mc_p_two_sided = min(1.0, 2.0 * min(mc_p_enrichment, mc_p_depletion))

    return {
        "scope": scope,
        "scope_role": scope_role,
        "analysis": analysis,
        "n_domains": int(probabilities.size),
        "observed_terminal_domains": int(observed),
        "observed_fraction": float(observed / probabilities.size),
        "expected_terminal_domains": expected,
        "expected_fraction": float(expected / probabilities.size),
        "observed_expected_ratio": float(observed / expected) if expected > 0 else np.nan,
        "median_domain_probability": float(np.median(probabilities)),
        "domains_with_probability_1": int(np.sum(probabilities == 1.0)),
        "domains_with_probability_ge_0_9": int(np.sum(probabilities >= 0.9)),
        "exact_null_mode": int(np.argmax(exact_pmf)),
        "exact_95pct_lower": exact_lower,
        "exact_95pct_upper": exact_upper,
        "exact_p_enrichment": exact_p_enrichment,
        "exact_p_depletion": exact_p_depletion,
        "exact_p_two_sided": exact_p_two_sided,
        "mc_null_median": float(np.median(simulated)),
        "mc_95pct_lower": float(mc_lower),
        "mc_95pct_upper": float(mc_upper),
        "mc_p_enrichment": mc_p_enrichment,
        "mc_p_depletion": mc_p_depletion,
        "mc_p_two_sided": mc_p_two_sided,
    }


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    normalized = series.astype(str).str.strip().str.lower()
    return normalized.isin({"true", "t", "1", "yes"})


# RELEASE: accepted-input guards ------------------------------------------------
def _release_verify_in_manifest(path: str) -> None:
    p = os.path.abspath(path)
    d = os.path.dirname(p)
    while not os.path.exists(os.path.join(d, "checksums.sha256")):
        if d == os.path.dirname(d):
            raise ValueError(f"RELEASE: no checksums.sha256 above {path}")
        d = os.path.dirname(d)
    rel = os.path.relpath(p, d)
    want = None
    for line in open(os.path.join(d, "checksums.sha256")):
        h, f = line.rstrip("\n").split(None, 1)
        if f in (rel, "./" + rel):
            want = h
    if want is None:
        raise ValueError(f"RELEASE: {rel} is not listed in {d}/checksums.sha256")
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if got != want:
        raise ValueError(f"RELEASE: {path} does not match its accepted digest")
    print(f"[RELEASE] {rel} sha256 {got} (manifest)", flush=True)


def _release_guard_inputs(args: argparse.Namespace) -> None:
    for p in (args.domains, args.associations):
        if os.path.basename(os.path.dirname(os.path.abspath(p))) != "conserved":
            raise ValueError("RELEASE: the positional null uses the order-30 conserved profile only")
    if "gene_exons_merged" in args.tx_exons_bed or "gene_introns" in args.tx_exons_bed:
        raise ValueError("RELEASE: transcript-resolved exons required (D022 gene-level assets are not an input)")
    for p in (args.domains, args.associations, args.tx_exons_bed):
        _release_verify_in_manifest(p)


def main() -> int:
    args = parse_args()
    _release_guard_inputs(args)  # RELEASE
    final_outdir = args.outdir  # RELEASE: publish atomically at the end
    os.makedirs(os.path.dirname(os.path.abspath(final_outdir)) or ".", exist_ok=True)
    args.outdir = tempfile.mkdtemp(prefix=".stage_positional_null.",
                                   dir=os.path.dirname(os.path.abspath(final_outdir)))
    if not 0 < args.min_exon_overlap <= 1:
        raise ValueError("--min-exon-overlap must be within (0, 1]")
    if args.n_permutations < 1:
        raise ValueError("--n-permutations must be positive")

    os.makedirs(args.outdir, exist_ok=True)

    domains = pd.read_csv(args.domains, sep="\t")
    associations = pd.read_csv(
        args.associations,
        sep="\t",
        dtype={"gene_id": str, "tx_id": str},
    )

    required_domain_columns = {
        "domain_uid",
        "chrom",
        "strand",
        "domain_length_bp",
        "Locus",
        "n_multi_exon_transcripts",
        "multi_exon_classification",
        "terminal_in_any_qualifying_transcript",
        "terminal_in_all_qualifying_transcripts",
    }
    missing = sorted(required_domain_columns - set(domains.columns))
    if missing:
        raise ValueError("Domain table is missing columns: " + ", ".join(missing))

    required_association_columns = {"domain_uid", "gene_id", "tx_id", "n_exons"}
    missing = sorted(required_association_columns - set(associations.columns))
    if missing:
        raise ValueError("Association table is missing columns: " + ", ".join(missing))

    if domains["domain_uid"].duplicated().any():
        raise ValueError("Domain table contains duplicated domain_uid values")

    domains["terminal_in_any_qualifying_transcript"] = as_bool(
        domains["terminal_in_any_qualifying_transcript"]
    )
    domains["terminal_in_all_qualifying_transcripts"] = as_bool(
        domains["terminal_in_all_qualifying_transcripts"]
    )

    domain_to_genes = (
        associations.groupby("domain_uid")["gene_id"]
        .apply(lambda x: sorted(set(x.dropna().astype(str))))
        .to_dict()
    )
    domain_to_observed_txs = (
        associations.groupby("domain_uid")["tx_id"]
        .apply(lambda x: sorted(set(x.dropna().astype(str))))
        .to_dict()
    )
    domain_to_observed_multi_txs = (
        associations.loc[pd.to_numeric(associations["n_exons"], errors="coerce") > 1]
        .groupby("domain_uid")["tx_id"]
        .apply(lambda x: sorted(set(x.dropna().astype(str))))
        .to_dict()
    )

    missing_gene_links = sorted(set(domains["domain_uid"]) - set(domain_to_genes))
    missing_tx_links = sorted(set(domains["domain_uid"]) - set(domain_to_observed_txs))
    if missing_gene_links or missing_tx_links:
        raise ValueError(
            "Domains missing associations: "
            + ", ".join((missing_gene_links + missing_tx_links)[:10])
        )

    relevant_genes = set(associations["gene_id"].dropna().astype(str))
    transcripts, gene_to_transcripts = load_transcripts(
        args.tx_exons_bed, relevant_genes
    )
    print(
        f"[INFO] Loaded {len(transcripts)} transcripts from "
        f"{len(gene_to_transcripts)} genes",
        flush=True,
    )

    scopes = []
    if args.scope in {"both", "observed"}:
        scopes.append(("observed_qualifying_transcripts", "primary"))
    if args.scope in {"both", "gene_all"}:
        scopes.append(("all_isoforms_of_associated_genes", "sensitivity"))

    per_domain_rows = []
    for scope, _role in scopes:
        print(f"[INFO] Enumerating scope: {scope}", flush=True)
        for index, domain_row in domains.iterrows():
            uid = domain_row["domain_uid"]
            if scope == "observed_qualifying_transcripts":
                tx_ids = domain_to_observed_txs[uid]
            else:
                tx_ids = []
                for gene_id in domain_to_genes[uid]:
                    tx_ids.extend(gene_to_transcripts.get(gene_id, []))

            result = enumerate_domain_null(
                domain_row=domain_row,
                tx_ids=tx_ids,
                transcripts=transcripts,
                min_exon_overlap=args.min_exon_overlap,
                scope=scope,
            )
            result["associated_gene_ids"] = ";".join(domain_to_genes[uid])
            result["observed_qualifying_transcript_ids"] = ";".join(
                domain_to_observed_txs[uid]
            )
            per_domain_rows.append(result)

            if (index + 1) % 25 == 0 or index + 1 == len(domains):
                print(
                    f"[INFO] {scope}: processed {index + 1}/{len(domains)} domains",
                    flush=True,
                )

    per_domain = pd.DataFrame(per_domain_rows)
    per_domain_path = os.path.join(
        args.outdir, "domain_centered_terminal_exon_null_probabilities_v2.tsv"
    )
    per_domain.to_csv(per_domain_path, sep="\t", index=False)

    rng = np.random.default_rng(args.seed)
    summaries = []
    permutation_tables = []

    def run_analysis(
        scope: str,
        role: str,
        name: str,
        domain_subset: pd.DataFrame,
        probability_column: str,
        observed: int,
    ):
        scope_prob = per_domain.loc[
            per_domain["scope"] == scope,
            ["domain_uid", probability_column],
        ]
        merged = domain_subset[["domain_uid"]].merge(
            scope_prob, on="domain_uid", how="left", validate="one_to_one"
        )
        probabilities = merged[probability_column].to_numpy(dtype=float)
        if np.isnan(probabilities).any():
            raise ValueError(f"Missing probabilities for {scope} / {name}")

        simulated = simulate_counts(probabilities, args.n_permutations, rng)
        summaries.append(
            summarize_null(
                scope=scope,
                scope_role=role,
                analysis=name,
                observed=observed,
                probabilities=probabilities,
                simulated=simulated,
            )
        )
        permutation_tables.append(
            pd.DataFrame(
                {
                    "permutation": np.arange(1, args.n_permutations + 1),
                    "scope": scope,
                    "analysis": name,
                    "terminal_domain_count": simulated,
                }
            )
        )

    for scope, role in scopes:
        run_analysis(
            scope,
            role,
            "all_domains_terminal_in_any_transcript",
            domains,
            "p_terminal_any",
            int(domains["terminal_in_any_qualifying_transcript"].sum()),
        )
        run_analysis(
            scope,
            role,
            "all_domains_terminal_in_all_transcripts",
            domains,
            "p_terminal_all",
            int(domains["terminal_in_all_qualifying_transcripts"].sum()),
        )
        run_analysis(
            scope,
            role,
            "all_domains_stringent_terminal_in_any_transcript",
            domains,
            "p_stringent_terminal_any",
            int(domains["terminal_in_any_qualifying_transcript"].sum()),
        )

        multi = domains.loc[domains["n_multi_exon_transcripts"] > 0].copy()
        observed_multi_any = int(
            (
                multi["multi_exon_classification"].astype(str)
                != "nonterminal_all_multi_exon"
            ).sum()
        )
        observed_multi_all = int(
            (
                multi["multi_exon_classification"].astype(str)
                == "terminal_all_multi_exon"
            ).sum()
        )

        # For the observed-transcript primary scope, recompute the multi-exon
        # probabilities using only observed qualifying multi-exon transcripts.
        if scope == "observed_qualifying_transcripts":
            multi_rows = []
            for _, domain_row in multi.iterrows():
                uid = domain_row["domain_uid"]
                tx_ids = domain_to_observed_multi_txs.get(uid, [])
                if not tx_ids:
                    raise ValueError(f"Missing observed multi-exon transcripts for {uid}")
                multi_result = enumerate_domain_null(
                    domain_row=domain_row,
                    tx_ids=tx_ids,
                    transcripts=transcripts,
                    min_exon_overlap=args.min_exon_overlap,
                    scope="observed_qualifying_multi_exon_transcripts",
                )
                multi_rows.append(multi_result)
            multi_prob = pd.DataFrame(multi_rows)

            for name, probability_column, observed_value in [
                (
                    "multi_exon_domains_terminal_in_any_transcript",
                    "p_terminal_any",
                    observed_multi_any,
                ),
                (
                    "multi_exon_domains_terminal_in_all_transcripts",
                    "p_terminal_all",
                    observed_multi_all,
                ),
                (
                    "multi_exon_domains_stringent_terminal_in_any_transcript",
                    "p_stringent_terminal_any",
                    observed_multi_any,
                ),
            ]:
                probabilities = multi_prob[probability_column].to_numpy(dtype=float)
                simulated = simulate_counts(probabilities, args.n_permutations, rng)
                summaries.append(
                    summarize_null(
                        scope="observed_qualifying_transcripts",
                        scope_role="primary",
                        analysis=name,
                        observed=observed_value,
                        probabilities=probabilities,
                        simulated=simulated,
                    )
                )
                permutation_tables.append(
                    pd.DataFrame(
                        {
                            "permutation": np.arange(1, args.n_permutations + 1),
                            "scope": scope,
                            "analysis": name,
                            "terminal_domain_count": simulated,
                        }
                    )
                )
        else:
            run_analysis(
                scope,
                role,
                "multi_exon_domains_terminal_in_any_transcript",
                multi,
                "p_terminal_any_multi_exon",
                observed_multi_any,
            )
            run_analysis(
                scope,
                role,
                "multi_exon_domains_terminal_in_all_transcripts",
                multi,
                "p_terminal_all_multi_exon",
                observed_multi_all,
            )
            run_analysis(
                scope,
                role,
                "multi_exon_domains_stringent_terminal_in_any_transcript",
                multi,
                "p_stringent_terminal_any_multi_exon",
                observed_multi_any,
            )

    summary = pd.DataFrame(summaries)
    summary_path = os.path.join(
        args.outdir, "domain_centered_terminal_exon_null_summary_v2.tsv"
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    permutations = pd.concat(permutation_tables, ignore_index=True)
    permutations_path = os.path.join(
        args.outdir, "domain_centered_terminal_exon_permutation_counts_v2.tsv.gz"
    )
    permutations.to_csv(permutations_path, sep="\t", index=False,
                        compression={"method": "gzip", "mtime": 0})  # RELEASE: fixed gzip mtime

    scope_diag = (
        per_domain.groupby("scope", as_index=False)
        .agg(
            n_domains=("domain_uid", "nunique"),
            minimum_valid_positions_per_domain=("n_valid_exonic_starts", "min"),
            median_valid_positions_per_domain=("n_valid_exonic_starts", "median"),
            median_transcript_models=("n_transcript_models", "median"),
        )
    )
    scope_diag_path = os.path.join(
        args.outdir, "domain_centered_terminal_exon_null_scope_diagnostics_v2.tsv"
    )
    scope_diag.to_csv(scope_diag_path, sep="\t", index=False)

    qc = pd.DataFrame(
        {
            "metric": [
                "input_domains",
                "unique_herv_loci",
                "multi_domain_herv_loci",
                "qualifying_associations",
                "transcript_models_loaded",
                "genes_loaded",
                "null_scopes",
                "n_permutations",
                "seed",
            ],
            "value": [
                len(domains),
                domains["Locus"].nunique(),
                int((domains.groupby("Locus").size() > 1).sum()),
                len(associations),
                len(transcripts),
                len(gene_to_transcripts),
                ";".join(scope for scope, _role in scopes),
                args.n_permutations,
                args.seed,
            ],
        }
    )
    qc_path = os.path.join(
        args.outdir, "domain_centered_terminal_exon_null_qc_v2.tsv"
    )
    qc.to_csv(qc_path, sep="\t", index=False)

    # RELEASE: atomic publication of the staged directory
    for name in os.listdir(args.outdir):
        os.chmod(os.path.join(args.outdir, name), 0o644)
    if os.path.exists(final_outdir):
        shutil.rmtree(final_outdir)
    os.rename(args.outdir, final_outdir)
    os.chmod(final_outdir, 0o755)
    staged = args.outdir
    args.outdir = final_outdir
    per_domain_path, summary_path, permutations_path, scope_diag_path, qc_path = (
        p.replace(staged, final_outdir) for p in
        (per_domain_path, summary_path, permutations_path, scope_diag_path, qc_path))

    print("\n=== TERMINAL-EXON NULL SUMMARY V2 ===")
    display_cols = [
        "scope_role",
        "scope",
        "analysis",
        "n_domains",
        "observed_terminal_domains",
        "observed_fraction",
        "expected_terminal_domains",
        "expected_fraction",
        "exact_95pct_lower",
        "exact_95pct_upper",
        "exact_p_enrichment",
    ]
    print(summary[display_cols].to_string(index=False))
    print("\nOutputs:")
    for path in [
        per_domain_path,
        summary_path,
        permutations_path,
        scope_diag_path,
        qc_path,
    ]:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
