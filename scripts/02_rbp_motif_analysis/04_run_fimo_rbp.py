#!/usr/bin/env python3
"""
04_run_fimo_rbp.py

Run FIMO to scan strand-oriented transcript-like FASTA sequences (e.g., HERV internals)
for RBP PWMs (ATtRACT/CisBP-RNA exported to MEME format).

Key defaults for RNA-like scanning:
- --norc ON by default (single-strand scanning)
- PGC (--parse-genomic-coord) OFF by default because custom headers often confuse it
- Optionally clean sequence_name (pos_chr11 -> chr11) after run

Supports scanning:
- sense only (default)
- antisense only (requires --antisense-fasta, which should be reverse-complemented sequences)
- both (runs twice and writes into outdir/sense and outdir/antisense)

Background handling:
- sense uses --bgfile (or --bgfile-sense if provided)
- antisense uses --bgfile-antisense (required for mode=antisense|both)

Outputs per run:
- <outdir>/<sense|antisense>/fimo.tsv
- <outdir>/<sense|antisense>/fimo.raw.tsv   (backup if cleaning changes anything)
- <outdir>/<sense|antisense>/cmd.txt        (command provenance)
"""

import argparse
import os
import re
import sys
import subprocess
from pathlib import Path


def require_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        print(f"❌ {label} not found: {path}", file=sys.stderr)
        sys.exit(1)


def run(cmd) -> None:
    print(f"🔧 Command: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def clean_sequence_name(name: str) -> str:
    m = re.match(r"^(pos|neg|plus|minus)_(chr[0-9A-Za-z]+)$", name)
    if m:
        return m.group(2)
    return name


def clean_fimo_tsv_sequence_names(tsv_path: Path) -> bool:
    if not tsv_path.exists():
        raise RuntimeError(f"Missing expected output: {tsv_path}")

    raw_path = tsv_path.with_name("fimo.raw.tsv")
    tmp_path = tsv_path.with_suffix(".tmp")
    changed = False

    with tsv_path.open() as inp, tmp_path.open("w") as out:
        header = inp.readline()
        if not header:
            raise RuntimeError(f"Empty TSV: {tsv_path}")
        out.write(header)

        for line in inp:
            line = line.rstrip("\n")
            if not line:
                out.write("\n")
                continue
            fields = line.split("\t")
            if len(fields) >= 3:
                orig = fields[2]
                new = clean_sequence_name(orig)
                if new != orig:
                    fields[2] = new
                    changed = True
                out.write("\t".join(fields) + "\n")
            else:
                out.write(line + "\n")

    if changed:
        if raw_path.exists():
            raw_path.unlink()
        tsv_path.replace(raw_path)
        tmp_path.replace(tsv_path)
    else:
        tmp_path.unlink(missing_ok=True)

    return changed


def validate_fimo_tsv(tsv_path: Path, enforce_plus_only: bool = False) -> None:
    if not tsv_path.exists():
        raise RuntimeError(f"Missing expected output: {tsv_path}")

    with tsv_path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        expected = [
            "motif_id", "motif_alt_id", "sequence_name", "start", "stop",
            "strand", "score", "p-value", "q-value", "matched_sequence"
        ]
        if header[:len(expected)] != expected:
            print(f"⚠️  Unexpected TSV header in {tsv_path}:\n  {header}", file=sys.stderr)

        bad_nf = 0
        seen_strands = set()
        for i, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 10:
                bad_nf += 1
                if bad_nf <= 5:
                    print(f"⚠️  Line {i} has {len(fields)} columns (expected 10): {line}", file=sys.stderr)
                continue
            seen_strands.add(fields[5])

        if bad_nf > 0:
            raise RuntimeError(
                f"FIMO TSV appears malformed: {bad_nf} lines do not have 10 columns. "
                f"Example issues printed above."
            )

        if enforce_plus_only and (seen_strands - {"+"}):
            raise RuntimeError(
                f"Unexpected strand values in FIMO TSV despite --norc: {sorted(seen_strands)}"
            )


def build_fimo_cmd(args, fasta_path: str, outdir: Path, bgfile: str) -> list[str]:
    cmd = [
        args.fimo_path,
        "--thresh", str(args.thresh),
        "--max-stored-scores", str(args.max_scores),
        "--bgfile", bgfile,
        "--oc", str(outdir),
    ]

    if args.qv_thresh:
        cmd.append("--qv-thresh")

    if not args.pgc:
        cmd.append("--no-pgc")

    if args.norc:
        cmd.append("--norc")

    cmd += [args.motifs, fasta_path]
    return cmd


def run_one(args, label: str, fasta_path: str, base_outdir: Path, bgfile: str) -> None:
    outdir = base_outdir / label
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== ▶ Running mode: {label} ===", file=sys.stderr)
    print(f"📄 FASTA: {fasta_path}", file=sys.stderr)
    print(f"🌫️  BGFILE: {bgfile}", file=sys.stderr)

    cmd = build_fimo_cmd(args, fasta_path=fasta_path, outdir=outdir, bgfile=bgfile)

    with open(outdir / "cmd.txt", "w") as fh:
        fh.write(" ".join(cmd) + "\n")

    run(cmd)

    tsv_path = outdir / "fimo.tsv"

    if args.clean_seqname:
        changed = clean_fimo_tsv_sequence_names(tsv_path)
        if changed:
            print("🧼 Cleaned sequence_name in fimo.tsv (backup saved as fimo.raw.tsv).", file=sys.stderr)

    if args.validate:
        validate_fimo_tsv(tsv_path, enforce_plus_only=args.enforce_plus_only and args.norc)

    print(f"✅ Done. Results: {tsv_path}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run FIMO for RBP motif scanning on transcript-oriented FASTA."
    )
    ap.add_argument("--motifs", "-m", required=True, help="MEME motif file (e.g., attract_hsa_cisbp.meme)")
    ap.add_argument("--fasta", "-f", required=True, help="Sense FASTA to scan (e.g., ERV_internal.fasta)")
    ap.add_argument("--outdir", "-o", required=True, help="Base output directory (created if needed)")

    ap.add_argument("--mode", choices=["sense", "antisense", "both"], default="sense",
                    help="Which RNA orientation to scan. Default: sense.")
    ap.add_argument("--antisense-fasta", default=None,
                    help="FASTA for antisense scanning (reverse-complemented sequences). Required if mode=antisense/both.")

    ap.add_argument("--thresh", "-t", type=float, default=1e-4,
                    help="Threshold (p-value by default; q-value if --qv-thresh). Default: 1e-4")
    ap.add_argument("--qv-thresh", action="store_true",
                    help="Use q-value thresholding instead of p-value.")
    ap.add_argument("--max-scores", type=int, default=1_000_000,
                    help="Max motif matches to store. Default: 1e6")

    # Backgrounds
    ap.add_argument("--bgfile", required=True,
                    help="Background file for sense scan (default).")
    ap.add_argument("--bgfile-sense", default=None,
                    help="Optional explicit background for sense scan (overrides --bgfile).")
    ap.add_argument("--bgfile-antisense", default=None,
                    help="Background file for antisense scan (required if mode=antisense|both).")

    ap.add_argument("--fimo-path", default="fimo", help="Path to FIMO executable.")

    ap.add_argument("--norc", dest="norc", action="store_true", default=True,
                    help="Do not score reverse complement (RNA-like). Default: ON.")
    ap.add_argument("--allow-rc", dest="norc", action="store_false",
                    help="Allow reverse-complement scanning (DNA-style).")

    ap.add_argument("--pgc", action="store_true", default=False,
                    help="Enable parse-genomic-coord (PGC). OFF by default.")

    ap.add_argument("--clean-seqname", dest="clean_seqname", action="store_true", default=True,
                    help="Clean sequence_name (pos_chr11 -> chr11) in fimo.tsv. Default: ON.")
    ap.add_argument("--no-clean-seqname", dest="clean_seqname", action="store_false",
                    help="Disable cleaning of sequence_name.")

    ap.add_argument("--validate", dest="validate", action="store_true", default=True,
                    help="Validate fimo.tsv integrity after run. Default: ON.")
    ap.add_argument("--no-validate", dest="validate", action="store_false",
                    help="Disable validation of fimo.tsv.")
    ap.add_argument("--enforce-plus-only", dest="enforce_plus_only", action="store_true", default=True,
                    help="Require only '+' strand hits (expected with --norc). Default: ON.")
    ap.add_argument("--no-enforce-plus-only", dest="enforce_plus_only", action="store_false",
                    help="Do not enforce '+'-only hits.")

    args = ap.parse_args()

    require_file(args.motifs, "Motif file")
    require_file(args.fasta, "Sense FASTA file")
    require_file(args.bgfile, "Sense background file (--bgfile)")

    # Resolve sense bg
    sense_bg = args.bgfile_sense if args.bgfile_sense else args.bgfile
    require_file(sense_bg, "Sense background file (resolved)")

    # Antisense requirements
    if args.mode in ("antisense", "both"):
        if not args.antisense_fasta:
            print("❌ --antisense-fasta is required when --mode is antisense or both.", file=sys.stderr)
            sys.exit(1)
        require_file(args.antisense_fasta, "Antisense FASTA file")

        if not args.bgfile_antisense:
            print("❌ --bgfile-antisense is required when --mode is antisense or both.", file=sys.stderr)
            sys.exit(1)
        require_file(args.bgfile_antisense, "Antisense background file (--bgfile-antisense)")

    base_outdir = Path(args.outdir)
    base_outdir.mkdir(parents=True, exist_ok=True)

    # Info banner
    if args.qv_thresh:
        print(f"🔬 Using q-value threshold: {args.thresh}", file=sys.stderr)
    else:
        print(f"🔬 Using p-value threshold: {args.thresh}", file=sys.stderr)

    if not args.pgc:
        print("🧩 PGC disabled: --no-pgc (recommended for custom FASTA headers).", file=sys.stderr)
    else:
        print("🧩 PGC enabled.", file=sys.stderr)

    if args.norc:
        print("🧬 RNA-like mode: --norc enabled (single-strand scanning).", file=sys.stderr)
    else:
        print("🧬 DNA-like mode: reverse-complement scanning allowed.", file=sys.stderr)

    try:
        if args.mode == "sense":
            run_one(args, "sense", args.fasta, base_outdir, bgfile=sense_bg)
        elif args.mode == "antisense":
            run_one(args, "antisense", args.antisense_fasta, base_outdir, bgfile=args.bgfile_antisense)
        else:  # both
            run_one(args, "sense", args.fasta, base_outdir, bgfile=sense_bg)
            run_one(args, "antisense", args.antisense_fasta, base_outdir, bgfile=args.bgfile_antisense)

    except subprocess.CalledProcessError:
        print("❌ FIMO failed.", file=sys.stderr)
        sys.exit(2)
    except RuntimeError as e:
        print(f"❌ Post-run check failed: {e}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
