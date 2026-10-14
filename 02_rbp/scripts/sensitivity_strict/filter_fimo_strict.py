#!/usr/bin/env python3
"""
Strict-threshold FIMO filter (Stage 11).

Streams one primary FIMO TSV and emits the strict-threshold file: the same file
with data rows removed unless their p-value satisfies `float(field 8) <= 1e-5`.

    keep a data row  <=>  float(field 8) <= 1e-5

FIMO IS NOT RERUN. Decision D016 established that the historical strict datasets
are a pure row filter of the primary p <= 1e-4 outputs, applied to the raw
fimo.tsv before any parsing. This script reimplements that and nothing else.

FOUR PROPERTIES THAT DECIDE BYTE IDENTITY
-----------------------------------------
1. `<=`, NOT `<`. Rows sitting exactly on 1e-05 exist in every dataset (34 in
   observed). Using `<` drops them, changes the observed file from 1,672,920 to
   1,669,142 bytes and yields the wrong sha256. The operators are not
   interchangeable here.
2. RETAINED ROWS ARE COPIED VERBATIM. The line is never rebuilt from parsed
   fields. Parsing happens only to read field 8 for the decision; the bytes
   written are the bytes read. That preserves score and q-value notation, tabs,
   and every one of the ten columns exactly.
3. HEADER AND TRAILER ARE EMITTED UNCONDITIONALLY. FIMO writes a 4-line trailer:
   one blank line then three comment lines (version, format URL, original
   command). Any empty line or line starting with '#' passes through untouched.
   The trailer therefore still cites the UPSTREAM run's --thresh 0.0001 and its
   --oc path. That is correct and deliberate: the filtering step is recorded in
   separate provenance rather than by editing the trailer.
4. ROW ORDER IS UNTOUCHED. The output is an ordered subsequence of the input.

Binary mode throughout, so no encoding or newline normalisation can occur.

Usage:
    filter_fimo_strict.py --input IN.tsv --output OUT.tsv [--expect-rows N]
                          [--expect-strict N] [--expect-boundary N]
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

# The one place the threshold is defined.
STRICT_P_THRESHOLD = 1e-5
P_VALUE_FIELD = 8          # 1-based; index 7 in a zero-based split
EXPECTED_COLUMNS = 10
EXPECTED_TRAILER_LINES = 4
OUTPUT_MODE = 0o644


def die(message: str) -> "None":
    sys.stderr.write(f"[error] {message}\n")
    sys.exit(1)


def filter_stream(source, sink):
    """Filter one FIMO TSV. Returns counts; never rewrites a retained row."""
    stats = {
        "header": 0, "trailer": 0, "data": 0, "kept": 0,
        "dropped": 0, "boundary": 0, "passthrough_malformed": 0,
    }
    seen_header = False

    for raw in source:
        stripped = raw.rstrip(b"\r\n")

        # 1. Header: the first line, emitted unconditionally.
        if not seen_header:
            seen_header = True
            stats["header"] += 1
            sink.write(raw)
            continue

        # 2. Trailer / comments: empty lines and '#' lines, unconditionally.
        if not stripped or stripped.startswith(b"#"):
            stats["trailer"] += 1
            sink.write(raw)
            continue

        # 3. Data rows.
        fields = stripped.split(b"\t")
        if len(fields) < P_VALUE_FIELD:
            # Defensive: emit verbatim rather than silently dropping. No such
            # line occurs in these six files.
            stats["passthrough_malformed"] += 1
            sink.write(raw)
            continue

        try:
            p_value = float(fields[P_VALUE_FIELD - 1])
        except ValueError:
            stats["passthrough_malformed"] += 1
            sink.write(raw)
            continue

        stats["data"] += 1
        if p_value <= STRICT_P_THRESHOLD:
            stats["kept"] += 1
            if p_value == STRICT_P_THRESHOLD:
                stats["boundary"] += 1
            sink.write(raw)          # verbatim, not reconstructed
        else:
            stats["dropped"] += 1

    return stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--expect-rows", type=int, default=None,
                    help="expected primary data-row count")
    ap.add_argument("--expect-strict", type=int, default=None,
                    help="expected strict data-row count")
    ap.add_argument("--expect-boundary", type=int, default=None,
                    help="expected rows with p-value exactly 1e-5")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    if not in_path.is_file():
        die(f"input not found: {in_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Atomic publication: stage beside the destination (same filesystem, so the
    # replace is atomic), publish only after every gate passes.
    fd, tmp_name = tempfile.mkstemp(dir=str(out_path.parent),
                                    prefix=".strict.", suffix=".tmp")
    tmp = Path(tmp_name)
    try:
        with open(fd, "wb") as sink, in_path.open("rb") as source:
            stats = filter_stream(source, sink)

        if stats["header"] != 1:
            raise RuntimeError(f"expected 1 header line, saw {stats['header']}")
        if stats["trailer"] != EXPECTED_TRAILER_LINES:
            raise RuntimeError(
                f"expected {EXPECTED_TRAILER_LINES} trailer lines "
                f"(one blank plus three comments), saw {stats['trailer']}")
        if args.expect_rows is not None and stats["data"] != args.expect_rows:
            raise RuntimeError(
                f"input has {stats['data']} data rows, expected {args.expect_rows}")
        if args.expect_strict is not None and stats["kept"] != args.expect_strict:
            raise RuntimeError(
                f"kept {stats['kept']} rows, expected {args.expect_strict}")
        if args.expect_boundary is not None and stats["boundary"] != args.expect_boundary:
            raise RuntimeError(
                f"{stats['boundary']} rows sit exactly on the threshold, "
                f"expected {args.expect_boundary}. A '<' comparison would give 0.")
        if stats["passthrough_malformed"]:
            sys.stderr.write(
                f"[warn] {stats['passthrough_malformed']} malformed lines "
                f"emitted verbatim\n")

        os.chmod(tmp, OUTPUT_MODE)
        tmp.replace(out_path)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise

    print(f"{in_path.name}: {stats['data']} data rows -> {stats['kept']} kept "
          f"({stats['boundary']} exactly at {STRICT_P_THRESHOLD:g}), "
          f"{stats['dropped']} dropped; header 1, trailer {stats['trailer']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
