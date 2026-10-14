#!/usr/bin/env python3
"""
Canonical QC of the dinucleotide-preserving shuffled FASTA replicates.

Step 20 of the frozen execution plan; clean-release successor of the historical
`scripts/rbp/06b_qc_dinuc_shuffles.py`, which recomputed QC from the shuffled
FASTAs on disk rather than trusting what the generator reported.

This is the standalone QC artifact. It re-reads the original and every replicate
from disk and recomputes every metric from scratch. It shares no code with the
Stage 5 producer, and it does not read Stage 5's validation output.

WHAT IS CHECKED, per locus, for every locus in every replicate
--------------------------------------------------------------
  identifier      preserved, modulo the documented `__dinuc_shufNN` tag
  description     preserved verbatim
  length          preserved
  first base      preserved
  last base       preserved
  mononucleotide  counts equal per symbol
  dinucleotide    counts equal per ordered pair

DINUCLEOTIDE DEFINITION
-----------------------
Overlapping adjacent ORDERED pairs:

    seq[i:i+2]   for i = 0 .. len(seq) - 2

Every position except the last contributes one pair, so a sequence of length L
has exactly L-1 pairs. Pairs containing `N` are counted like any other pair --
`N` is an ordinary symbol here, not a masked position. A non-overlapping
(step-2) interpretation is NOT used and would not preserve the property the null
model depends on.

UNCHANGED SEQUENCES
-------------------
A locus is "unchanged" when its shuffled sequence is character-for-character
equal to the original. This is expected at a low rate and is NOT a failure: a
short homopolymer or simple repeat may admit only one Eulerian trail, in which
case returning it unchanged is the only result consistent with preserving its
dinucleotide counts. It is reported, not flagged.

A no-op implementation would pass every preservation check above, so the report
also records how many loci actually changed.

Usage:
    qc_dinuc_shuffles.py --original IN.fasta --fasta-dir DIR --out-tsv REPORT.tsv
                         --seeds 01:26073001,02:26073002,...
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from collections import Counter
from pathlib import Path

# Canonical symbol order, matching config `dinuc_shuffle_allowed_bases`.
# Reported in this order rather than sorted(), which would give "ACGNT".
CANONICAL_ALPHABET = "ACGTN"
EXPECTED_ALPHABET = set(CANONICAL_ALPHABET)


def die(message: str) -> None:
    sys.stderr.write(f"[error] {message}\n")
    sys.exit(1)


def read_fasta(path: Path):
    """Yield (identifier, description, uppercase sequence). Fails on malformed input."""
    identifier = None
    description = ""
    parts: list[str] = []
    seen_any = False

    with path.open("r", encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                seen_any = True
                if identifier is not None:
                    yield identifier, description, "".join(parts).upper()
                fields = line[1:].split(maxsplit=1)
                if not fields:
                    raise ValueError(f"{path.name}: empty header at line {line_number}")
                identifier = fields[0]
                description = fields[1] if len(fields) > 1 else ""
                parts = []
            else:
                if identifier is None:
                    raise ValueError(
                        f"{path.name}: sequence before the first header at line {line_number}"
                    )
                parts.append(line.strip())

    if not seen_any:
        raise ValueError(f"{path.name}: no FASTA records found")
    if identifier is not None:
        yield identifier, description, "".join(parts).upper()


def dinucleotides(sequence: str) -> "Counter[str]":
    """Overlapping adjacent ordered pairs; see the module docstring."""
    return Counter(sequence[i:i + 2] for i in range(len(sequence) - 1))


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_seeds(text: str) -> "dict[str, int]":
    mapping = {}
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        rep, _, seed = item.partition(":")
        if not seed.isdigit():
            die(f"malformed --seeds entry: {item!r} (expected NN:SEED)")
        mapping[rep.strip()] = int(seed)
    if not mapping:
        die("--seeds parsed to an empty mapping")
    return mapping


COLUMNS = [
    "replicate", "seed", "n_records", "total_bases",
    "n_changed", "n_unchanged", "pct_changed", "pct_unchanged",
    "identifier_violations", "description_violations", "length_violations",
    "first_base_violations", "last_base_violations",
    "mononucleotide_violations", "dinucleotide_violations",
    "alphabet", "n_count", "unexpected_symbols",
    "n_dinucleotide_types", "bytes", "n_lines", "md5", "sha256", "status",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original", required=True)
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--seeds", required=True,
                    help="comma-separated NN:SEED, e.g. 01:26073001,02:26073002")
    ap.add_argument("--pattern", default="HERV_internal_sequences_dinuc_shuffle_{rep}.fasta")
    ap.add_argument("--tag-template", default="__dinuc_shuf{rep}")
    ap.add_argument("--expect-records", type=int, default=None)
    ap.add_argument("--expect-bases", type=int, default=None)
    args = ap.parse_args()

    seeds = parse_seeds(args.seeds)
    fasta_dir = Path(args.fasta_dir)
    original_path = Path(args.original)

    if not original_path.is_file():
        die(f"original FASTA not found: {original_path}")

    # Refuse to report on an incomplete replicate set: a QC artifact covering four
    # of five replicates would look like a pass.
    missing = [rep for rep in sorted(seeds)
               if not (fasta_dir / args.pattern.format(rep=rep)).is_file()]
    if missing:
        die(f"missing replicate FASTA(s): {', '.join(missing)}")

    print(f"reading original: {original_path}")
    original = {}
    for identifier, description, sequence in read_fasta(original_path):
        if identifier in original:
            die(f"duplicate identifier in the original FASTA: {identifier}")
        original[identifier] = (description, sequence)

    src_alphabet: "Counter[str]" = Counter()
    for _, sequence in original.values():
        src_alphabet.update(sequence)
    src_bases = sum(len(s) for _, s in original.values())

    print(f"  records {len(original):,}  bases {src_bases:,}  "
          f"alphabet {dict(sorted(src_alphabet.items()))}")

    unexpected_src = sorted(set(src_alphabet) - EXPECTED_ALPHABET)
    if unexpected_src:
        die(f"original FASTA contains unexpected symbols: {unexpected_src}")

    rows = []
    overall_ok = True

    for rep in sorted(seeds):
        path = fasta_dir / args.pattern.format(rep=rep)
        tag = args.tag_template.format(rep=rep)
        print(f"\nreplicate {rep}  seed {seeds[rep]}  {path.name}")

        v = Counter()
        n_records = 0
        n_bases = 0
        n_unchanged = 0
        alphabet: "Counter[str]" = Counter()
        dinuc_types: set[str] = set()
        seen: set[str] = set()

        for identifier, description, shuffled in read_fasta(path):
            if not identifier.endswith(tag):
                v["identifier"] += 1
                base_id = identifier
            else:
                base_id = identifier[: -len(tag)]

            if base_id in seen:
                die(f"replicate {rep}: duplicate identifier {identifier}")
            seen.add(base_id)

            entry = original.get(base_id)
            if entry is None:
                v["identifier"] += 1
                n_records += 1
                continue
            orig_desc, orig_seq = entry

            if description != orig_desc:
                v["description"] += 1
            if len(shuffled) != len(orig_seq):
                v["length"] += 1
            if shuffled[0:1] != orig_seq[0:1]:
                v["first_base"] += 1
            if shuffled[-1:] != orig_seq[-1:]:
                v["last_base"] += 1
            if Counter(shuffled) != Counter(orig_seq):
                v["mononucleotide"] += 1

            shuffled_dinucs = dinucleotides(shuffled)
            if shuffled_dinucs != dinucleotides(orig_seq):
                v["dinucleotide"] += 1
            dinuc_types.update(shuffled_dinucs)

            if shuffled == orig_seq:
                n_unchanged += 1

            alphabet.update(shuffled)
            n_bases += len(shuffled)
            n_records += 1

        # Every original locus must appear exactly once.
        for identifier in original:
            if identifier not in seen:
                v["identifier"] += 1

        unexpected = sorted(set(alphabet) - EXPECTED_ALPHABET)
        if unexpected:
            v["unexpected_symbol"] += len(unexpected)
        if alphabet != src_alphabet:
            v["mononucleotide"] += 1     # file-level symbol loss or gain

        n_changed = n_records - n_unchanged
        if n_records and n_changed == 0:
            v["no_randomization"] += 1

        if args.expect_records is not None and n_records != args.expect_records:
            v["expected_records"] += 1
        if args.expect_bases is not None and n_bases != args.expect_bases:
            v["expected_bases"] += 1

        status = "PASS" if not v else "FAIL"
        if status == "FAIL":
            overall_ok = False
            print(f"  FAIL  violations: {dict(v)}")
        else:
            print(f"  PASS  {n_records:,} records  {n_bases:,} bases  "
                  f"{n_unchanged} unchanged  {len(dinuc_types)} dinucleotide types")

        data = path.read_bytes()
        rows.append({
            "replicate": rep,
            "seed": seeds[rep],
            "n_records": n_records,
            "total_bases": n_bases,
            "n_changed": n_changed,
            "n_unchanged": n_unchanged,
            "pct_changed": f"{100 * n_changed / n_records:.4f}" if n_records else "NA",
            "pct_unchanged": f"{100 * n_unchanged / n_records:.4f}" if n_records else "NA",
            "identifier_violations": v.get("identifier", 0),
            "description_violations": v.get("description", 0),
            "length_violations": v.get("length", 0),
            "first_base_violations": v.get("first_base", 0),
            "last_base_violations": v.get("last_base", 0),
            "mononucleotide_violations": v.get("mononucleotide", 0),
            "dinucleotide_violations": v.get("dinucleotide", 0),
            "alphabet": "".join(c for c in CANONICAL_ALPHABET if c in alphabet),
            "n_count": alphabet.get("N", 0),
            "unexpected_symbols": ",".join(unexpected) if unexpected else "none",
            "n_dinucleotide_types": len(dinuc_types),
            "bytes": len(data),
            "n_lines": data.count(b"\n"),
            "md5": hashlib.md5(data).hexdigest(),
            "sha256": hashlib.sha256(data).hexdigest(),
            "status": status,
        })

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with tmp.open("w", encoding="ascii") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    tmp.chmod(0o644)
    tmp.replace(out)

    print(f"\nQC report: {out}")
    print(f"unchanged vector: {','.join(str(r['n_unchanged']) for r in rows)}")

    if not overall_ok:
        print("\nQC FAILED")
        return 1
    print("\nQC PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
