#!/usr/bin/env python3
"""
06_lncrna_function step 1 (frozen global order 38, effective order 39), LOCAL:
"Retrieve PubMed abstracts" -- SNAPSHOT-ONLY in the release (D031). No network access.

REFACTOR_SCRIPT of analysis/scripts/scripts_categories/F_literature_mining/01_fetch_abstracts_by_year.py
(sha256 379775f3...). The historical script retrieved a time-dependent corpus from NCBI PubMed with
Biopython Bio.Entrez (version UNKNOWN):
  query     ("RNA, Long Noncoding"[Mesh] OR "lncRNA"[TIAB] OR "lncRNAs"[TIAB] OR "long noncoding RNA"[TIAB]
             OR "long non-coding RNA"[TIAB] OR "lincRNA"[TIAB] OR "lincRNAs"[TIAB] OR "ceRNA"[TIAB])
            -- one topic query; no per-gene search (lncRNA genes are found later by text mining)
  slicing   esearch per publication year (datetype pdat, mindate = maxdate = year), 2000 .. current
            year (2026 at retrieval), usehistory; count capped at 9,999 per year
  fetch     efetch XML in batches of 500, 3 attempts, a batch failing 3 times silently skipped, 1 s pause
  record    PMID, Year (= the query slice, not a parsed date), ArticleTitle, AbstractText parts
            joined by one space; records without abstract text are NOT written
  output    CSV (csv.DictWriter, QUOTE_MINIMAL), one row per (year slice, PMID) hit
Re-running it today would return a different corpus, so the release freezes the corpus the canonical
v5 extraction consumed (SRC16, 01_sources/frozen/pubmed/lncrna_abstracts_dataset_full_no_function_filter.csv,
retrieved 2026-03-11) and this step only VERIFIES and INDEXES it:
  * sockets are disabled for the whole run (no code path can reach NCBI);
  * the corpus must match its SRC16 registration and must not be the SUPERSEDED file;
  * structure: header PMID,Year,Title,Abstract; numeric PMID; 4-digit year; non-empty title and abstract;
  * outputs (staged, atomic): corpus_manifest.tsv, corpus_year_counts.tsv, corpus_record_index.tsv
    (per-row PMID, Year and sha256 fingerprints of Title and Abstract).
Consumers (effective order 40) read the frozen corpus from 01_sources, verified by its sha256.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import socket
import sys
import tempfile
from collections import defaultdict


# ---------------------------------------------------------------- network prohibition guard
def _no_network(*_args, **_kwargs):
    raise RuntimeError("network access is disabled in the release build (D031: SRC16 is a frozen snapshot)")


socket.socket = _no_network                    # type: ignore[assignment]
socket.create_connection = _no_network         # type: ignore[assignment]
socket.getaddrinfo = _no_network               # type: ignore[assignment]

HISTORICAL_QUERY = ('("RNA, Long Noncoding"[Mesh] OR "lncRNA"[TIAB] OR "lncRNAs"[TIAB] OR "long noncoding RNA"[TIAB] '
                    'OR "long non-coding RNA"[TIAB] OR "lincRNA"[TIAB] OR "lincRNAs"[TIAB] OR "ceRNA"[TIAB])')
HEADER = ["PMID", "Year", "Title", "Abstract"]


def digest(path: str, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify_registration(corpus: str, sums_path: str, snapshots_path: str) -> str:
    name = os.path.basename(corpus)
    with open(snapshots_path, newline="", encoding="utf-8") as fh:
        status = [r["status"] for r in csv.DictReader(fh, delimiter="\t") if r["source_id"] == "SRC16" and r["snapshot_name"] == name]
    if status != ["AVAILABLE"]:
        raise ValueError(f"{name} is not the AVAILABLE SRC16 snapshot (status {status}); the superseded corpus is refused")
    with open(sums_path, newline="", encoding="utf-8") as fh:
        want = [r["sha256"] for r in csv.DictReader(fh, delimiter="\t") if r["source_id"] == "SRC16" and r["filename"] == name]
    got = digest(corpus, "sha256")
    if want != [got]:
        raise ValueError(f"{corpus} does not match its SRC16 registration")
    return got


def load_corpus(path: str) -> list:
    csv.field_size_limit(sys.maxsize)
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if header != HEADER:
            raise ValueError(f"unexpected corpus header {header}")
        rows = [r for r in reader]
    for i, r in enumerate(rows, start=1):
        if len(r) != 4:
            raise ValueError(f"row {i}: {len(r)} fields")
        pmid, year, title, abstract = r
        if not pmid.isdigit() or not (len(year) == 4 and year.isdigit()) or not title or not abstract:
            raise ValueError(f"row {i}: malformed record (PMID {pmid!r}, Year {year!r})")
    return rows


def sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description="Verify and index the frozen SRC16 PubMed corpus (no network).")
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--snapshot-checksums", required=True, help="01_sources/manifests/manual_snapshot_checksums.tsv")
    ap.add_argument("--snapshots", required=True, help="01_sources/manifests/manual_snapshots.tsv")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    sha256 = verify_registration(a.corpus, a.snapshot_checksums, a.snapshots)
    print(f"[ok] SRC16 {os.path.basename(a.corpus)} sha256 {sha256} (registered, AVAILABLE)", flush=True)
    rows = load_corpus(a.corpus)
    years = defaultdict(list)
    slices = defaultdict(set)
    for pmid, year, _t, _a in rows:
        years[year].append(pmid)
        slices[pmid].add(year)
    multi = sum(1 for s in slices.values() if len(s) > 1)
    dup_same_year = len(rows) - sum(len(set(v)) for v in years.values())
    manifest = [
        ("source_id", "SRC16"), ("source", "NCBI PubMed"), ("file", os.path.basename(a.corpus)), ("sha256", sha256),
        ("md5", digest(a.corpus, "md5")), ("bytes", str(os.path.getsize(a.corpus))), ("rows", str(len(rows))),
        ("unique_pmids", str(len(slices))), ("pmids_in_more_than_one_year_slice", str(multi)),
        ("duplicate_rows_within_a_year_slice", str(dup_same_year)), ("year_min", min(years)), ("year_max", max(years)),
        ("rows_without_abstract", "0 (such records were not written historically)"),
        ("historical_query", HISTORICAL_QUERY), ("historical_slicing", "esearch per publication year, datetype pdat, 2000..2026"),
        ("historical_per_year_cap", "9999 (maximum observed rows in a year: %d)" % max(len(v) for v in years.values())),
        ("historical_fetch", "efetch retmode xml, batches of 500, 3 attempts, failed batch skipped silently"),
        ("historical_tool", "Biopython Bio.Entrez (version UNKNOWN); NCBI E-utilities"),
        ("retrieval_date", "2026-03-11 (APPROXIMATE: output written 15:32; start time UNKNOWN)"),
        ("raw_entrez_responses", "NOT_RECOVERED (never saved)"),
        ("snapshot_nature", "ORIGINAL historical output of 01_fetch_abstracts_by_year.py (not a reconstruction)"),
        ("network_used_by_this_step", "NO"),
    ]
    parent = os.path.dirname(os.path.abspath(a.outdir)) or "."
    os.makedirs(parent, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".stage_pubmed_corpus.", dir=parent)
    with open(os.path.join(stage, "corpus_manifest.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["field", "value"])
        w.writerows(manifest)
    with open(os.path.join(stage, "corpus_year_counts.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Year", "rows", "unique_pmids"])
        for y in sorted(years):
            w.writerow([y, len(years[y]), len(set(years[y]))])
    with open(os.path.join(stage, "corpus_record_index.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["row", "PMID", "Year", "title_sha256", "abstract_sha256", "title_chars", "abstract_chars"])
        for i, (pmid, year, title, abstract) in enumerate(rows, start=1):
            w.writerow([i, pmid, year, sha(title), sha(abstract), len(title), len(abstract)])
    for n in os.listdir(stage):
        os.chmod(os.path.join(stage, n), 0o644)
    if os.path.exists(a.outdir):
        shutil.rmtree(a.outdir)
    os.rename(stage, a.outdir)
    os.chmod(a.outdir, 0o755)
    print(f"[ok] {len(rows)} rows, {len(slices)} PMIDs, years {min(years)}-{max(years)}, "
          f"{multi} PMIDs in more than one year slice; network used: NO")
    print(f"[OK] wrote {a.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
