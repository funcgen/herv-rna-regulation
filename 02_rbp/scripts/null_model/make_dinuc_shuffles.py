#!/usr/bin/env python3
"""
Generate ONE dinucleotide-preserving shuffled FASTA replicate.

Clean-release reimplementation of the historical canonical producer
`scripts/rbp/06_make_one_dinuc_shuffle.py` (forensic decisions D002 and D005).
The randomization is reproduced exactly; only path handling, configuration,
validation and error reporting are improved.

ALGORITHM (unchanged from the historical implementation)
--------------------------------------------------------
Each sequence is a directed multigraph whose vertices are nucleotide symbols and
whose edges are the observed adjacent dinucleotides. A randomized Eulerian trail
(Hierholzer) uses every edge exactly once, so the shuffled sequence has exactly
the same dinucleotide counts -- and therefore the same mononucleotide counts and
the same length -- as the original. Starting the trail at the original first
nucleotide preserves the first base; the trail's structure preserves the last.

FOUR PROPERTIES THAT MUST NOT CHANGE, because they determine the exact output
----------------------------------------------------------------------------
1. ONE RNG FOR THE WHOLE FILE. `random.Random(seed)` is created once and shared
   by every record, in input order. Record N's shuffle therefore depends on every
   preceding record. Re-seeding per record, or reordering records, would preserve
   every dinucleotide property and still produce different sequences.
2. ADJACENCY INSERTION ORDER. `rng.shuffle` is applied to each adjacency list in
   `dict.values()` order, which is insertion order -- the order in which each
   left-hand symbol is first seen in the sequence. Sorting the vertices would
   change the stream of random numbers consumed.
3. POP FROM THE END. Hierholzer takes the next edge with `list.pop()`, i.e. the
   last element of the shuffled adjacency list.
4. UPPERCASING ON READ. Sequences are uppercased before shuffling, so a
   lowercase-masked input contributes to the same vertex as its uppercase form.

SYMBOLS
-------
`ALLOWED_BASES = ACGTN`. N is an ordinary vertex and participates in the shuffle;
it is NOT fixed in place. This is the behaviour that actually ran, and it differs
from the abandoned MEME `fasta-shuffle-letters` plan, which specified `-fix N`
(decisions D003 and D005).

Usage:
    make_dinuc_shuffles.py --input IN.fasta --output OUT.fasta --seed N --tag TAG
"""

from __future__ import annotations

import argparse
import os
import random
import stat
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator, TextIO

# Historical: ALLOWED_BASES = set("ACGTN"). N participates in the shuffle.
ALLOWED_BASES = set("ACGTN")

# Release artifacts must be group/world readable.
OUTPUT_MODE = 0o644


def die(message: str) -> "None":
    sys.stderr.write(f"[error] {message}\n")
    sys.exit(1)


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header text, uppercase sequence). Historical parser, unchanged."""
    header: "str | None" = None
    parts: list[str] = []

    with path.open("r", encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header = line[1:]
                parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"sequence before the first FASTA header at line {line_number}"
                    )
                parts.append(line)

    if header is not None:
        yield header, "".join(parts).upper()


def append_tag(header: str, tag: str) -> str:
    """Append the replicate tag to the identifier, preserving any description."""
    fields = header.split(maxsplit=1)
    tagged = fields[0] + tag
    return tagged if len(fields) == 1 else tagged + " " + fields[1]


def dinucleotide_counts(sequence: str) -> "Counter[str]":
    """Overlapping adjacent pairs: sequence[i:i+2] for i in 0 .. len-2."""
    return Counter(sequence[i:i + 2] for i in range(len(sequence) - 1))


def dinucleotide_shuffle(sequence: str, rng: random.Random) -> str:
    """Randomized Eulerian trail (Hierholzer). Historical logic, unchanged."""
    if len(sequence) < 2:
        return sequence

    adjacency: "dict[str, list[str]]" = defaultdict(list)
    for left, right in zip(sequence[:-1], sequence[1:]):
        adjacency[left].append(right)

    # Insertion order, NOT sorted order -- see property 2 in the module docstring.
    for destinations in adjacency.values():
        rng.shuffle(destinations)

    stack = [sequence[0]]
    reverse_path: list[str] = []

    while stack:
        current = stack[-1]
        if adjacency[current]:
            stack.append(adjacency[current].pop())   # property 3: from the end
        else:
            reverse_path.append(stack.pop())

    shuffled = "".join(reversed(reverse_path))

    if len(shuffled) != len(sequence):
        raise RuntimeError(
            f"Eulerian trail has incorrect length: {len(sequence)} -> {len(shuffled)}"
        )

    return shuffled


def write_wrapped(handle: TextIO, sequence: str, line_width: int) -> None:
    for start in range(0, len(sequence), line_width):
        handle.write(sequence[start:start + line_width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate one dinucleotide-preserving shuffled FASTA replicate."
    )
    parser.add_argument("--input", required=True, help="canonical HERV internal FASTA")
    parser.add_argument("--output", required=True, help="canonical shuffled FASTA")
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--tag", required=True, help="appended to each identifier")
    parser.add_argument("--line-width", type=int, default=100)
    parser.add_argument("--expect-records", type=int, default=None)
    parser.add_argument("--expect-bases", type=int, default=None)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    # Staging file next to the destination, so the publish below is a same-filesystem
    # rename and therefore atomic. Not mkstemp: that creates 0600.
    temporary_path = Path(str(output_path) + ".tmp")

    if not input_path.is_file():
        die(f"input FASTA not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if temporary_path.exists():
        temporary_path.unlink()

    # One RNG for the whole file -- see property 1 in the module docstring.
    rng = random.Random(args.seed)

    n_records = 0
    n_bases = 0
    n_unchanged = 0
    symbols: "Counter[str]" = Counter()

    try:
        with temporary_path.open("w", encoding="ascii") as out:
            for header, original in read_fasta(input_path):
                invalid = set(original) - ALLOWED_BASES
                if invalid:
                    raise ValueError(
                        f"invalid characters in {header}: {sorted(invalid)}"
                    )

                symbols.update(original)
                shuffled = dinucleotide_shuffle(original, rng)

                # Per-record verification before anything is written, exactly as
                # the historical producer did. A failure aborts the whole replicate.
                if len(shuffled) != len(original):
                    raise RuntimeError(
                        f"length changed for {header}: "
                        f"{len(original)} -> {len(shuffled)}"
                    )
                if shuffled[0:1] != original[0:1]:
                    raise RuntimeError(f"first nucleotide changed for {header}")
                if shuffled[-1:] != original[-1:]:
                    raise RuntimeError(f"last nucleotide changed for {header}")
                if Counter(shuffled) != Counter(original):
                    raise RuntimeError(f"mononucleotide counts changed for {header}")
                if dinucleotide_counts(shuffled) != dinucleotide_counts(original):
                    raise RuntimeError(f"dinucleotide counts changed for {header}")

                if shuffled == original:
                    n_unchanged += 1

                out.write(">" + append_tag(header, args.tag) + "\n")
                write_wrapped(out, shuffled, args.line_width)

                n_records += 1
                n_bases += len(original)

                if n_records % 5000 == 0:
                    print(f"Processed {n_records:,} records; {n_bases:,} bases", flush=True)

        # Publication gates. The historical producer checked neither; a truncated
        # input would have yielded a short, plausible-looking replicate.
        if args.expect_records is not None and n_records != args.expect_records:
            raise RuntimeError(
                f"expected {args.expect_records:,} records, produced {n_records:,}"
            )
        if args.expect_bases is not None and n_bases != args.expect_bases:
            raise RuntimeError(
                f"expected {args.expect_bases:,} bases, produced {n_bases:,}"
            )

        os.chmod(temporary_path, OUTPUT_MODE)
        os.replace(temporary_path, output_path)

    except Exception:
        if temporary_path.exists():
            temporary_path.unlink()
        raise

    mode = stat.S_IMODE(output_path.stat().st_mode)

    print()
    print(f"Output: {output_path}")
    print(f"Seed: {args.seed}")
    print(f"Tag: {args.tag}")
    print(f"Records: {n_records:,}")
    print(f"Bases: {n_bases:,}")
    print(f"Symbols: {dict(sorted(symbols.items()))}")
    print(f"Unchanged sequences: {n_unchanged:,} ({100 * n_unchanged / n_records:.3f}%)")
    print(f"Permissions: {oct(mode)}")
    print("Per-record dinucleotide validation: PASS")


if __name__ == "__main__":
    main()
