#!/usr/bin/env python3
"""Read one value from config/analysis_profiles.yaml.

Standalone so shell callers never interpolate into a regex, and so a block ends
at the first line indented no deeper than its key -- the sibling-key bug found
in Stage 7's reader.
"""
import argparse
import re
import sys


def block_lines(text, key, min_indent=0):
    lines = text.split("\n")
    for i, line in enumerate(lines):
        m = re.match(r"^(\s*)" + re.escape(key) + r":\s*$", line)
        if not m or len(m.group(1)) < min_indent:
            continue
        indent = len(m.group(1))
        out = []
        for nxt in lines[i + 1:]:
            if not nxt.strip() or nxt.lstrip().startswith("#"):
                continue
            if len(nxt) - len(nxt.lstrip()) <= indent:
                break
            out.append(nxt)
        return out
    return []


def scalar(lines, key):
    for l in lines:
        m = re.match(r'^\s*' + re.escape(key) + r':\s*"?([^"#\n]+?)"?\s*(?:#.*)?$', l)
        if m:
            return m.group(1).strip()
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--profile")
    ap.add_argument("--key")
    ap.add_argument("--shared")
    args = ap.parse_args()

    text = open(args.file, encoding="utf-8").read()

    if args.shared:
        v = scalar(block_lines(text, "shared"), args.shared)
        if v is None:
            sys.stderr.write(f"shared key not found: {args.shared}\n")
            return 1
        sys.stdout.write(v)
        return 0

    prof = block_lines(text, args.profile, min_indent=1)
    if not prof:
        sys.stderr.write(f"profile not found: {args.profile}\n")
        return 1

    # Derived keys assembled from the expected_hits sub-block.
    if args.key in ("expect_observed", "expect_shuffles"):
        sub = "\n".join(prof)
        hits = block_lines(sub, "expected_hits")
        obs = scalar(hits, "observed")
        reps = dict(re.findall(r'^\s*"(\d+)":\s*(\d+)\s*$', "\n".join(hits), re.M))
        if args.key == "expect_observed":
            if obs is None:
                sys.stderr.write("expected_hits.observed not found\n")
                return 1
            sys.stdout.write(obs)
        else:
            if len(reps) != 5:
                sys.stderr.write(f"expected 5 shuffle counts, found {len(reps)}\n")
                return 1
            sys.stdout.write(",".join(f"{k}:{v}" for k, v in sorted(reps.items())))
        return 0

    v = scalar(prof, args.key)
    if v is None:
        sys.stderr.write(f"key not found in profile {args.profile}: {args.key}\n")
        return 1
    sys.stdout.write(v)
    return 0


if __name__ == "__main__":
    sys.exit(main())
