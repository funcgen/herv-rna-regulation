#!/usr/bin/env python3
"""
Independent validator for 03_transcript_context step 3 (effective order 10).

Step 3 changes ORDER ONLY. For every asset this checks, without sorting anything
itself and without calling bedtools:

  record set   output multiset of lines == input multiset (exact bytes), so
               row count, columns, coordinates, names, scores, strands and
               duplicates are all preserved; column count 6 everywhere
  order        chromosome rank (chrom-sizes file order) never decreases;
               start never decreases within a chromosome -- the complete
               historical key (bedtools sort -g sorts on chrom, then start)
  sequences    input and output sequence sets identical; every sequence is
               in the chrom-sizes file
  ties         reports (chrom,start) tie groups, identical-coordinate records
               with different names, and tie groups spanning both strands
  D022         for the corrected step-2 assets: 443,234 / 364,340 rows, one
               gene per interval, strand + or -, no mixed-gene / mixed-strand
               intervals; set preserved exactly from the accepted input

With --historical-root, each output is also compared with the historical sorted
asset: BYTE_IDENTICAL / SCIENTIFICALLY_IDENTICAL_SERIALIZATION_DIFFERS /
SCIENTIFICALLY_EQUIVALENT_ORDERING_SERIALIZATION_DIFFERS / SCIENTIFICALLY_DIFFERENT,
and for the two D022 assets INTENTIONALLY_DIFFERENT_D022_CORRECTION.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
from collections import Counter

ASSETS = [  # class, name, d022, anchor rows
    ("genes", "gencode.v48.genes", False, 78894),
    ("transcripts", "gencode.v48.transcripts", False, 387954),
    ("features", "gencode.v48.exons", False, 2165096),
    ("features", "gencode.v48.cds", False, 902508),
    ("features", "gencode.v48.utr", False, 391840),
    ("features", "gencode.v48.gene_exons_merged", True, 443234),
    ("features", "gencode.v48.gene_introns", True, 364340),
]
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{('  ' + detail) if detail else ''}")


def lines(path):
    with gzip.open(path, "rb") as fh:
        return fh.read().split(b"\n")[:-1]


def sha(data_lines):
    h = hashlib.sha256()
    for l in data_lines:
        h.update(l + b"\n")
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True, help="canonical results/gencode_features")
    ap.add_argument("--sorted-root", help="where <class>/sorted/ lives (default: results root)")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--historical-root", help="historical context/assets (read-only) for comparison")
    ap.add_argument("--report", help="write a per-asset TSV report here")
    a = ap.parse_args()
    sroot = a.sorted_root or a.results_root

    rank = {}
    for l in open(a.genome, "rb"):
        rank.setdefault(l.split(b"\t")[0], len(rank))

    report = []
    for cls, name, d022, anchor in ASSETS:
        src = f"{a.results_root}/{cls}/{name}.bed.gz"
        out = f"{sroot}/{cls}/sorted/{name}.sorted.bed.gz"
        print(f"\n{name}")
        A, B = lines(src), lines(out)
        check(f"{name}: row count preserved", len(A) == len(B), f"{len(A):,} -> {len(B):,}")
        check(f"{name}: row count == accepted anchor", len(B) == anchor, f"{anchor:,}")
        ca, cb = Counter(A), Counter(B)
        check(f"{name}: record multiset identical (coords, names, scores, strands)", ca == cb,
              f"{sum((ca - cb).values())} missing, {sum((cb - ca).values())} added")
        dup_a = sum(v - 1 for v in ca.values() if v > 1)
        dup_b = sum(v - 1 for v in cb.values() if v > 1)
        check(f"{name}: duplicate rows preserved", dup_a == dup_b, f"{dup_a} in, {dup_b} out")
        rows = [l.split(b"\t") for l in B]
        check(f"{name}: 6 columns on every row", all(len(r) == 6 for r in rows))
        unknown = {r[0] for r in rows if r[0] not in rank}
        check(f"{name}: every sequence in the chrom-sizes file", not unknown, ",".join(x.decode() for x in unknown))
        keys = [(rank.get(r[0], 10**9), int(r[1])) for r in rows]
        bad = sum(1 for i in range(1, len(keys)) if keys[i] < keys[i - 1])
        check(f"{name}: genome order (chrom rank, then start) non-decreasing", bad == 0, f"{bad} inversions")
        seq_in = {l.split(b"\t", 1)[0] for l in A}
        seq_out = [r[0] for r in rows]
        seq_blocks = [seq_out[0]] + [seq_out[i] for i in range(1, len(seq_out)) if seq_out[i] != seq_out[i - 1]]
        check(f"{name}: sequence set unchanged", seq_in == set(seq_out), f"{len(seq_in)} sequences")
        check(f"{name}: each sequence one contiguous block", len(seq_blocks) == len(set(seq_blocks)))
        # ties
        tie_cs, same_coord_diff_name, both_strands = Counter(), 0, 0
        grp = {}
        for r in rows:
            grp.setdefault((r[0], r[1]), []).append(r)
        n_tie = sum(1 for v in grp.values() if len(v) > 1)
        both_strands = sum(1 for v in grp.values() if len({x[5] for x in v}) > 1)
        cse = {}
        for r in rows:
            cse.setdefault((r[0], r[1], r[2]), set()).add(r[3])
        same_coord_diff_name = sum(1 for v in cse.values() if len(v) > 1)
        main_chroms = [c.decode() for c in (b"chrX", b"chrY", b"chrM") if c in seq_in]
        alt = sum(1 for c in seq_in if not c.startswith(b"chr"))
        rec = dict(asset=name, cls=cls, d022=d022, rows_in=len(A), rows_out=len(B),
                   content_sha256_in=sha(A), content_sha256_out=sha(B),
                   sequences=len(seq_in), first_sequence=seq_blocks[0].decode(), last_sequence=seq_blocks[-1].decode(),
                   chrX_chrY_chrM=",".join(main_chroms), alt_contigs=alt,
                   tie_groups_chrom_start=n_tie, same_coords_different_names=same_coord_diff_name,
                   tie_groups_both_strands=both_strands)
        if d022:
            mixed_gene = sum(1 for r in rows if b"," in r[3])
            mixed_strand = sum(1 for r in rows if r[5] not in (b"+", b"-"))
            want = 4 if name.endswith("introns") else 3
            badname = sum(1 for r in rows if len(r[3].split(b"|")) != want)
            check(f"{name}: D022 mixed-gene intervals == 0", mixed_gene == 0, str(mixed_gene))
            check(f"{name}: D022 mixed-strand intervals == 0", mixed_strand == 0, str(mixed_strand))
            check(f"{name}: D022 one gene identity per interval", badname == 0, str(badname))
        if a.historical_root:
            hist = f"{a.historical_root}/{cls}/sorted/{name}.sorted.bed.gz"
            H = lines(hist)
            same_gz = open(hist, "rb").read() == open(out, "rb").read()
            if d022:
                verdict = "INTENTIONALLY_DIFFERENT_D022_CORRECTION"
                ok = sha(H) != sha(B)  # the buggy historical content must NOT be reproduced
                check(f"{name}: historical {len(H):,} rows vs canonical {len(B):,}; {verdict}", ok)
            elif same_gz:
                verdict = "BYTE_IDENTICAL"
            elif H == B:
                verdict = "SCIENTIFICALLY_IDENTICAL_SERIALIZATION_DIFFERS"
            elif Counter(H) == cb:
                verdict = "SCIENTIFICALLY_EQUIVALENT_ORDERING_SERIALIZATION_DIFFERS"
            else:
                verdict = "SCIENTIFICALLY_DIFFERENT"
            if not d022:
                check(f"{name}: historical comparison", verdict != "SCIENTIFICALLY_DIFFERENT", verdict)
            rec.update(historical_rows=len(H), historical_content_sha256=sha(H), historical_verdict=verdict)
        report.append(rec)

    if a.report:
        cols = list(dict.fromkeys(k for r in report for k in r))
        with open(a.report, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in report:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    nfail = sum(1 for _, ok, _ in results if not ok)
    print(f"\n{'ALL STEP 3 VALIDATIONS PASSED' if nfail == 0 else f'{nfail} CHECK(S) FAILED'} ({len(results)} checks)")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
