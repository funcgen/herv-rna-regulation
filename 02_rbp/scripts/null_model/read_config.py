#!/usr/bin/env python3
"""Read one scientific parameter out of config/parameters.yaml.

Standalone so that shell callers never interpolate anything into a regex.
Prints the value and exits non-zero if the key is absent, so a missing
parameter fails the caller loudly instead of silently becoming an empty string.
"""
import argparse
import re
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True)
    ap.add_argument("--key")
    ap.add_argument("--seed-for", dest="seed_for")
    ap.add_argument("--hits-for", dest="hits_for")
    args = ap.parse_args()

    text = open(args.file, encoding="utf-8").read()

    if args.key:
        m = re.search(
            r"^\s*" + re.escape(args.key) + r":\s*\"?([^\"#\n]+?)\"?\s*(?:#.*)?$",
            text, re.M)
        if not m:
            sys.stderr.write(f"parameter not found: {args.key}\n")
            return 1
        sys.stdout.write(m.group(1).strip())
        return 0

    if args.seed_for:
        block = re.search(r"^\s*dinuc_shuffle_seeds:\s*\n(.*?)(?=^\S|\Z)",
                          text, re.M | re.S)
        if not block:
            sys.stderr.write("no dinuc_shuffle_seeds block\n")
            return 1
        m = re.search(r"^\s*\"" + re.escape(args.seed_for) + r"\":\s*(\d+)\s*$",
                      block.group(1), re.M)
        if not m:
            sys.stderr.write(f"no seed configured for replicate {args.seed_for}\n")
            return 1
        sys.stdout.write(m.group(1))
        return 0

    if args.hits_for:
        block = re.search(r"^\s*fimo_null_expected_hits:\s*\n(.*?)(?=^\S|\Z)",
                          text, re.M | re.S)
        if not block:
            sys.stderr.write("no fimo_null_expected_hits block\n")
            return 1
        m = re.search(r"^\s*\"" + re.escape(args.hits_for) + r"\":\s*(\d+)\s*$",
                      block.group(1), re.M)
        if not m:
            sys.stderr.write(f"no expected hit count for replicate {args.hits_for}\n")
            return 1
        sys.stdout.write(m.group(1))
        return 0

    ap.error("pass --key, --seed-for or --hits-for")
    return 2


if __name__ == "__main__":
    sys.exit(main())
