#!/usr/bin/env python3
"""
03_transcript_context step 2 (effective order 9): per-gene merged exons and
per-gene introns.

Canonical definition (D022, scientific implementation correction):

    for each gene_id independently:
        E(g) = union of all annotated exons of g, across all its transcripts,
               merging overlapping OR directly touching intervals
        I(g) = gene body of g  minus  E(g)

Grouping is by gene identity, explicitly. Exons of other genes, on either
strand, never merge with or mask the current gene -- including same-strand
overlapping genes, which a strand-only grouping (bedtools merge -s) would
still fuse.

This replaces the historical genome-wide `bedtools merge` / `bedtools subtract`,
which merged exons across genes and strands. That behaviour is preserved only as
a forensic reproduction under validation/historical_reproduction/step2/.

Inputs  (accepted step-1 BED6, 0-based half-open):
  genes: name = gene_id|gene_name|gene_type
  exons: name = transcript_id|gene_id|gene_name|exon|gene_type|transcript_type
Outputs (BED6, gzip, genome order):
  gencode.v48.gene_exons_merged.bed.gz   name = gene_id|gene_name|gene_type
  gencode.v48.gene_introns.bed.gz        name = gene_id|gene_name|gene_type|intron
  score 0; strand = the gene's strand; exactly one gene per interval.

Genome order: chromosome rank from the chrom-sizes file, then start, end,
gene_id. The gzip header carries mtime 0 and no file name, so reruns are
byte-identical.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tempfile


def die(msg):
    sys.exit(f"[error] {msg}")


def read_bed6(path):
    with gzip.open(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 6:
                die(f"{path}:{ln}: expected 6 columns, found {len(f)}")
            yield ln, f[0], int(f[1]), int(f[2]), f[3], f[5]


def union_touching(intervals):
    """Merge overlapping or book-ended [start, end) intervals of one gene."""
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def complement(body_start, body_end, merged):
    """Gene body minus that gene's own merged exons (merged: sorted, disjoint)."""
    out, cursor = [], body_start
    for s, e in merged:
        if s > cursor:
            out.append((cursor, s))
        cursor = max(cursor, e)
    if cursor < body_end:
        out.append((cursor, body_end))
    return out


def write_bed_gz(path, rows):
    """Deterministic gzip: mtime 0, no embedded file name."""
    h = hashlib.sha256()
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        for r in rows:
            line = ("\t".join(map(str, r)) + "\n").encode()
            gz.write(line)
            h.update(line)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--genes", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--genome", required=True, help="chrom sizes (defines chromosome order)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    rank = {}
    with open(a.genome) as fh:
        for line in fh:
            c = line.split("\t")[0]
            rank.setdefault(c, len(rank))

    genes = {}
    for ln, c, s, e, name, st in read_bed6(a.genes):
        p = name.split("|")
        if len(p) != 3:
            die(f"gene line {ln}: name {name!r} is not gene_id|gene_name|gene_type")
        if st not in ("+", "-"):
            die(f"gene line {ln}: strand {st!r}")
        if c not in rank:
            die(f"gene line {ln}: chromosome {c} absent from {a.genome}")
        if not s < e:
            die(f"gene line {ln}: empty interval")
        if p[0] in genes:
            die(f"duplicate gene_id {p[0]}")
        genes[p[0]] = (c, s, e, name, st)

    exons = {g: [] for g in genes}
    n_exon = 0
    for ln, c, s, e, name, st in read_bed6(a.exons):
        n_exon += 1
        p = name.split("|")
        if len(p) != 6 or p[3] != "exon":
            die(f"exon line {ln}: unexpected name {name!r}")
        gid = p[1]
        if gid not in genes:
            die(f"exon line {ln}: gene {gid} absent from the gene BED")
        gc, gs, ge, _, gst = genes[gid]
        if c != gc or st != gst or s < gs or e > ge or not s < e:
            die(f"exon line {ln}: not inside its own gene {gid} on the same chromosome/strand")
        exons[gid].append((s, e))

    merged_rows, intron_rows = [], []
    for gid, (c, s, e, name, st) in genes.items():
        if not exons[gid]:
            die(f"gene {gid} has no exons")
        m = union_touching(exons[gid])
        for x, y in m:
            merged_rows.append((c, x, y, name, 0, st))
        for x, y in complement(s, e, m):
            intron_rows.append((c, x, y, name + "|intron", 0, st))

    key = lambda r: (rank[r[0]], r[1], r[2], r[3].split("|")[0])
    merged_rows.sort(key=key)
    intron_rows.sort(key=key)
    if not merged_rows or not intron_rows:
        die("empty output")

    os.makedirs(a.outdir, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".stage_gene_introns.", dir=a.outdir)
    try:
        out = {}
        for fname, rows in (("gencode.v48.gene_exons_merged.bed.gz", merged_rows),
                            ("gencode.v48.gene_introns.bed.gz", intron_rows)):
            tmp = os.path.join(stage, fname)
            out[fname] = (write_bed_gz(tmp, rows), len(rows))
            os.chmod(tmp, 0o644)
        for fname in out:
            os.replace(os.path.join(stage, fname), os.path.join(a.outdir, fname))
    finally:
        for f in os.listdir(stage):
            os.remove(os.path.join(stage, f))
        os.rmdir(stage)

    print(f"[ok] genes {len(genes)}  source exons {n_exon}")
    for fname, (digest, n) in out.items():
        print(f"[ok] {fname}  rows {n}  decompressed sha256 {digest}")


if __name__ == "__main__":
    main()
