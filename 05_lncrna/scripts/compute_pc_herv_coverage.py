#!/usr/bin/env python3
"""
05_lncrna step 2 (frozen global order 35, effective order 36), LOCAL:
per-gene HERV base-pair coverage of the protein-coding genes that overlap a HERV internal
region on the same strand.

REFACTOR_SCRIPT of analysis/scripts/scripts_categories/B_transcript_context/03bis_find_pc_hervs.py
(sha256 f6eeb9ec...; identical to 03_find_lncrna_hervs.py apart from its paths), with the D028
correction (the D027 rule: HERV geometry is BED half-open).

  universe      every gene of the order-19 protein_coding_overlapping_HERVs.gene_level.tsv (one
                row per same-strand protein-coding host gene); the gene must exist in the GTF
  HERV bases    the gene's overlapping_coords, same chromosome and strand as the gene;
                D027/D028: these are BED 0-based half-open and are converted once to 1-based
                inclusive [start + 1, end] (the historical script read them as 1-based)
  union         HERV intervals merged per gene (overlapping and bookended blocks merge),
                so no base is counted twice
  gene body     the GTF 'gene' feature span, 1-based inclusive: gene_body_bp = end - start + 1
  exons         the union of every GTF 'exon' record of the gene: exonic_bp
  covered bp    |HERV union  n  gene span|  and  |HERV union  n  exon union|  (clipped)
  coverage      fraction = covered / denominator, round(., 4); percent = round(fraction * 100, 2)
  row order     descending (exonic_herv_percent, gene_body_herv_percent); ties keep the
                order-19 table order (stable sort)

Continuous coverage only: no categories are assigned (the frozen plan bins only the
long non-coding table, at effective order 37). No thresholds are read, so config/parameters.yaml
is not used. Kept separate from compute_lncrna_herv_coverage.py (order 35), as 03bis was from 03.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, Tuple

Interval = Tuple[int, int]          # 1-based inclusive

OUTPUT_FILE = "protein_coding_overlapping_herv_coverage.tsv"
FIELDNAMES = [
    "gene_stable_id", "gene_name", "gene_type",
    "chrom", "strand", "gene_start", "gene_end",
    "gene_body_bp", "exonic_bp", "n_overlapping_hervs",
    "herv_bp_in_gene_body", "herv_bp_in_exons",
    "gene_body_herv_fraction", "gene_body_herv_percent",
    "exonic_herv_fraction", "exonic_herv_percent",
    "overlapping_subfamilies",
]
COORD = re.compile(r"^([^:]+):(\d+)-(\d+)\(([+-])\)$")
GTF_ATTR = re.compile(r'(\S+)\s+"([^"]+)"')


# ---------------------------------------------------------------- input guards
def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify_against_manifest(path: str, manifest: str) -> str:
    rel = os.path.relpath(os.path.abspath(path), os.path.dirname(os.path.abspath(manifest)))
    want = None
    for line in open(manifest):
        h, f = line.rstrip("\n").split(None, 1)
        if f in (rel, "./" + rel):
            want = h
    got = sha256(path)
    if want is None or got != want:
        raise ValueError(f"{path} does not match its accepted manifest {manifest}")
    return got


# ---------------------------------------------------------------- intervals (1-based inclusive)
def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    """Union of 1-based inclusive intervals; overlapping and bookended intervals merge."""
    merged: List[List[int]] = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def interval_length(intervals: List[Interval]) -> int:
    return sum(e - s + 1 for s, e in intervals)


def intersect_merged(a: List[Interval], b: List[Interval]) -> List[Interval]:
    """Intersection of two merged interval lists (two-pointer sweep)."""
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s <= e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return merge_intervals(out)


def parse_bed_coord(token: str) -> Tuple[str, Interval, str]:
    """'chr1:259065-259773(+)', a BED 0-based half-open interval -> 1-based inclusive (D027/D028)."""
    m = COORD.match(token.strip())
    if not m:
        raise ValueError(f"Could not parse coordinate: {token!r}")
    chrom, start, end, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    if not start < end:
        raise ValueError(f"Empty or inverted BED interval: {token!r}")
    return chrom, (start + 1, end), strand


# ---------------------------------------------------------------- loaders
def load_pc_universe(path: str) -> List[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    ids = [r["gene_stable_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("The protein-coding overlap table lists a gene more than once")
    bad = sorted({r["gene_type"] for r in rows} - {"protein_coding"})
    if bad:
        raise ValueError(f"Non-protein-coding gene types in the protein-coding overlap table: {bad}")
    return rows


def load_gtf_genes(path: str, wanted: set) -> Tuple[Dict[str, dict], Dict[str, List[Interval]]]:
    """GTF 'gene' spans and merged 'exon' records for the wanted genes (1-based inclusive)."""
    genes: Dict[str, dict] = {}
    exons: Dict[str, List[Interval]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] not in ("gene", "exon"):
                continue
            m = re.search(r'gene_id "([^"]+)"', f[8])
            if not m or m.group(1) not in wanted:
                continue
            gene_id = m.group(1)
            if f[2] == "gene":
                attrs = dict(GTF_ATTR.findall(f[8]))
                if gene_id in genes:
                    raise ValueError(f"GTF lists gene {gene_id} twice")
                genes[gene_id] = {"chrom": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                                  "gene_type": attrs.get("gene_type")}
            else:
                exons[gene_id].append((int(f[3]), int(f[4])))
    return genes, {g: merge_intervals(v) for g, v in exons.items()}


# ---------------------------------------------------------------- coverage
def gene_coverage(row: dict, gene: dict, exon_union: List[Interval]) -> dict:
    herv = []
    for token in row["overlapping_coords"].split(";"):
        if not token.strip():
            continue
        chrom, iv, strand = parse_bed_coord(token)
        if chrom == gene["chrom"] and strand == gene["strand"]:      # same strand only
            herv.append(iv)
    herv_union = merge_intervals(herv)
    span = [(gene["start"], gene["end"])]
    gene_bp, exon_bp = interval_length(span), interval_length(exon_union)
    in_gene = interval_length(intersect_merged(herv_union, span))
    in_exons = interval_length(intersect_merged(herv_union, exon_union))
    body_f = in_gene / gene_bp if gene_bp > 0 else 0.0
    exon_f = in_exons / exon_bp if exon_bp > 0 else 0.0
    return {
        "gene_stable_id": row["gene_stable_id"], "gene_name": row["gene_name"], "gene_type": row["gene_type"],
        "chrom": gene["chrom"], "strand": gene["strand"], "gene_start": gene["start"], "gene_end": gene["end"],
        "gene_body_bp": gene_bp, "exonic_bp": exon_bp, "n_overlapping_hervs": row["n_overlapping_hervs"],
        "herv_bp_in_gene_body": in_gene, "herv_bp_in_exons": in_exons,
        "gene_body_herv_fraction": round(body_f, 4), "gene_body_herv_percent": round(body_f * 100, 2),
        "exonic_herv_fraction": round(exon_f, 4), "exonic_herv_percent": round(exon_f * 100, 2),
        "overlapping_subfamilies": row["overlapping_subfamilies"],
        "_n_herv_tokens_other_strand": sum(1 for t in row["overlapping_coords"].split(";") if t.strip()) - len(herv),
    }


def compute_coverage(universe: List[dict], genes: Dict[str, dict], exons: Dict[str, List[Interval]]) -> List[dict]:
    missing = [r["gene_stable_id"] for r in universe if r["gene_stable_id"] not in genes]
    if missing:
        raise ValueError(f"{len(missing)} universe genes absent from the GTF, e.g. {missing[:3]}")
    wrong = [g for g in genes if genes[g]["gene_type"] != "protein_coding"]
    if wrong:
        raise ValueError(f"GTF gene_type is not protein_coding for {wrong[:3]}")
    rows = [gene_coverage(r, genes[r["gene_stable_id"]], exons.get(r["gene_stable_id"], [])) for r in universe]
    return sorted(rows, key=lambda x: (x["exonic_herv_percent"], x["gene_body_herv_percent"]), reverse=True)


def write_outputs(rows: List[dict], outdir: str) -> None:
    parent = os.path.dirname(os.path.abspath(outdir)) or "."
    os.makedirs(parent, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".stage_pc_herv_coverage.", dir=parent)
    with open(os.path.join(stage, OUTPUT_FILE), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.chmod(os.path.join(stage, OUTPUT_FILE), 0o644)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.rename(stage, outdir)
    os.chmod(outdir, 0o755)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--overlap-table", required=True, help="order-19 protein_coding_overlapping_HERVs.gene_level.tsv")
    ap.add_argument("--overlap-manifest", required=True, help="03 results/gene_overlap/checksums.sha256")
    ap.add_argument("--gtf", required=True, help="GENCODE v48 primary-assembly GTF (SRC02)")
    ap.add_argument("--gtf-sha256", required=True, help="expected SRC02 sha256 (01_sources/SOURCE_STATUS.tsv)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    print(f"[ok] overlap table sha256 {verify_against_manifest(a.overlap_table, a.overlap_manifest)} (manifest)", flush=True)
    got = sha256(a.gtf)
    if got != a.gtf_sha256:
        raise ValueError(f"GTF sha256 {got} is not SRC02 {a.gtf_sha256}")
    print(f"[ok] GTF sha256 {got} (SRC02)", flush=True)

    universe = load_pc_universe(a.overlap_table)
    genes, exons = load_gtf_genes(a.gtf, {r["gene_stable_id"] for r in universe})
    rows = compute_coverage(universe, genes, exons)
    write_outputs(rows, a.outdir)
    other = sum(r["_n_herv_tokens_other_strand"] for r in rows)
    over = sum(r["herv_bp_in_gene_body"] > r["gene_body_bp"] for r in rows)
    print(f"[ok] {len(rows)} protein-coding genes; HERV tokens off the gene's chromosome/strand: {other}; "
          f"genes with covered bp > gene body: {over}")
    print(f"[OK] wrote {os.path.join(a.outdir, OUTPUT_FILE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
