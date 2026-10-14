#!/usr/bin/env python3
"""
Build the canonical six-dataset FIMO manifest (Stage 7).

This is the interface every downstream null-model stage reads: one row per
dataset, naming the exact FASTA, background, motif library and FIMO output that
produced it, with checksums.

The observed row is NOT rescanned. It points at the Stage 3 artifact and is
labelled `generated_in_stage = Stage 3`; the five shuffles are `Stage 7`.

Usage:
    build_fimo_dataset_manifest.py --module-dir 02_rbp --config config/parameters.yaml \
        --out results/null_model/fimo_dataset_manifest.tsv \
        --checksums-out validation/fimo_null_datarow_checksums.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

REPLICATES = ["01", "02", "03", "04", "05"]

MANIFEST_COLUMNS = [
    "dataset", "replicate", "seed", "fasta", "fasta_sha256",
    "background", "background_sha256", "motif_library", "motif_library_sha256",
    "fimo_output", "n_hits", "n_positive_loci", "scientific_rows_sha256",
    "generated_in_stage", "status",
]

CHECKSUM_COLUMNS = [
    "dataset", "replicate", "seed", "n_hits", "n_positive_loci",
    "scientific_rows_sha256", "full_file_md5", "full_file_sha256", "bytes", "status",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def scan_fimo(path: Path):
    """Scientific rows only: skip the header, skip '#' trailer lines, keep order."""
    h = hashlib.sha256()
    loci = set()
    n = 0
    with path.open("rb") as fh:
        for i, raw in enumerate(fh):
            if i == 0 or raw.startswith(b"#"):
                continue
            line = raw.rstrip(b"\n")
            if not line.strip():
                continue
            if line.count(b"\t") != 9:
                continue
            n += 1
            loci.add(line.split(b"\t")[2])
            h.update(raw)
    return n, len(loci), h.hexdigest()


def cfg_map(text: str, key: str) -> "dict[str, str]":
    """Read an indented `key:` mapping block.

    The block ends at the first following line indented no deeper than the key
    itself. A naive `(?=^\\S|\\Z)` terminator runs past a sibling key -- both
    `dinuc_shuffle_seeds` and `fimo_null_expected_hits` are nested under `rbp:`,
    so it swallowed the second mapping and returned hit counts as seeds.
    """
    lines = text.split("\n")
    for i, line in enumerate(lines):
        m = re.match(r"^(\s*)" + re.escape(key) + r":\s*$", line)
        if not m:
            continue
        indent = len(m.group(1))
        out: "dict[str, str]" = {}
        for nxt in lines[i + 1:]:
            if not nxt.strip() or nxt.lstrip().startswith("#"):
                continue
            if len(nxt) - len(nxt.lstrip()) <= indent:
                break
            e = re.match(r'^\s*"(\d+)":\s*(\d+)\s*$', nxt)
            if e:
                out[e.group(1)] = e.group(2)
        return out
    return {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--checksums-out", required=True)
    # Dataset layout. The defaults are the primary profile, so an existing
    # primary invocation is unchanged; the strict profile passes its own paths.
    # See config/analysis_profiles.yaml.
    ap.add_argument("--observed-fimo", default="results/fimo_primary/fimo.tsv")
    ap.add_argument("--shuffle-template",
                    default="results/null_model/fimo/shuffle{rep}/fimo.tsv")
    ap.add_argument("--analysis-label", default="primary")
    ap.add_argument("--generated-in-stage-observed", default="Stage 3")
    ap.add_argument("--generated-in-stage-shuffle", default="Stage 7")
    ap.add_argument("--expect-observed", type=int, default=716075)
    ap.add_argument("--expect-shuffles", default=None,
                    help="comma-separated NN:COUNT overriding the config block")
    args = ap.parse_args()

    mod = Path(args.module_dir)
    text = Path(args.config).read_text()
    seeds = cfg_map(text, "dinuc_shuffle_seeds")
    expected = cfg_map(text, "fimo_null_expected_hits")
    if args.expect_shuffles:
        expected = dict(item.split(":") for item in args.expect_shuffles.split(","))

    motifs = mod / "results/motif_library/attract_hsa_cisbp.meme"
    bg = mod / "results/background/herv_internal_order1_background.txt"
    observed_fasta = Path("01_sources/frozen/herv/HERV_internal_sequences_v5.fasta")
    obs_fasta_abs = mod.parent / observed_fasta

    for p in (motifs, bg):
        if not p.is_file():
            sys.exit(f"missing fixed input: {p}")

    motif_sha = sha256_file(motifs)
    bg_sha = sha256_file(bg)

    rows, csums = [], []

    # ---- observed: Stage 3, referenced, never regenerated ----
    obs_tsv = mod / args.observed_fimo
    if not obs_tsv.is_file():
        sys.exit(f"Stage 3 observed output missing: {obs_tsv}")
    n, loci, sha = scan_fimo(obs_tsv)
    rows.append({
        "dataset": "observed", "replicate": "NA", "seed": "NA",
        "fasta": str(observed_fasta),
        "fasta_sha256": sha256_file(obs_fasta_abs) if obs_fasta_abs.is_file() else "NA",
        "background": str(bg.relative_to(mod)), "background_sha256": bg_sha,
        "motif_library": str(motifs.relative_to(mod)), "motif_library_sha256": motif_sha,
        "fimo_output": str(obs_tsv.relative_to(mod)),
        "n_hits": n, "n_positive_loci": loci, "scientific_rows_sha256": sha,
        "generated_in_stage": args.generated_in_stage_observed,
        "status": "PASS" if n == args.expect_observed else "FAIL",
    })
    csums.append({
        "dataset": "observed", "replicate": "NA", "seed": "NA", "n_hits": n,
        "n_positive_loci": loci, "scientific_rows_sha256": sha,
        "full_file_md5": md5_file(obs_tsv), "full_file_sha256": sha256_file(obs_tsv),
        "bytes": obs_tsv.stat().st_size,
        "status": "PASS" if n == args.expect_observed else "FAIL",
    })

    # ---- the five shuffles: Stage 7 ----
    for rep in REPLICATES:
        fasta = mod / f"results/null_model/fasta/HERV_internal_sequences_dinuc_shuffle_{rep}.fasta"
        tsv = mod / args.shuffle_template.format(rep=rep)
        if not tsv.is_file():
            sys.exit(f"missing Stage 7 output: {tsv}")
        n, loci, sha = scan_fimo(tsv)
        want = int(expected.get(rep, -1))
        status = "PASS" if n == want else "FAIL"
        rows.append({
            "dataset": f"shuffle{rep}", "replicate": rep, "seed": seeds.get(rep, "NA"),
            "fasta": str(fasta.relative_to(mod)), "fasta_sha256": sha256_file(fasta),
            "background": str(bg.relative_to(mod)), "background_sha256": bg_sha,
            "motif_library": str(motifs.relative_to(mod)), "motif_library_sha256": motif_sha,
            "fimo_output": str(tsv.relative_to(mod)),
            "n_hits": n, "n_positive_loci": loci, "scientific_rows_sha256": sha,
            "generated_in_stage": args.generated_in_stage_shuffle, "status": status,
        })
        csums.append({
            "dataset": f"shuffle{rep}", "replicate": rep, "seed": seeds.get(rep, "NA"),
            "n_hits": n, "n_positive_loci": loci, "scientific_rows_sha256": sha,
            "full_file_md5": md5_file(tsv), "full_file_sha256": sha256_file(tsv),
            "bytes": tsv.stat().st_size, "status": status,
        })

    def write(path: Path, cols, data):
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".tmp")
        with tmp.open("w", encoding="ascii") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in data:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        tmp.chmod(0o644)
        tmp.replace(path)

    write(Path(args.out), MANIFEST_COLUMNS, rows)
    write(Path(args.checksums_out), CHECKSUM_COLUMNS, csums)

    print(f"manifest : {args.out}")
    print(f"checksums: {args.checksums_out}")
    for r in rows:
        print(f"  {r['dataset']:<10} hits={r['n_hits']:>7,}  loci={r['n_positive_loci']:>6,}  "
              f"{r['generated_in_stage']:<8} {r['status']}")

    bgs = {r["background_sha256"] for r in rows}
    print(f"\ndistinct background models across all six datasets: {len(bgs)}")
    if len(bgs) != 1:
        print("FAIL: more than one background model")
        return 1
    if any(r["status"] != "PASS" for r in rows):
        print("FAIL: at least one dataset missed its expected hit count")
        return 1
    print("manifest OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
