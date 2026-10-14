#!/usr/bin/env python3
"""
03_transcript_context: transcript-level TSS BED (D023 dependency repair, a
prerequisite of effective order 13).

Reimplements scripts/dsrna/bidirectional/01_extract_gencode_tss_bed.py, which
wrote the historical input/annotations/tss_gencode.transcripts.bed from the GTF
'transcript' features. The accepted step-1 transcript BED carries the same
records, in the same GTF order, already converted to BED coordinates
(start0 = start1 - 1, end = end1), so the GTF is not re-parsed.

Rule (one row per transcript; nothing filtered, collapsed or deduplicated):
  strand +  TSS = transcript start  -> BED [start0, start0 + 1)
  strand -  TSS = transcript end    -> BED [end - 1, end)
  name  = transcript_id|gene_id|gene_name   (Ensembl versions kept)
  score = 0, strand = the transcript's strand
Equivalently, in GTF terms: + -> [start1-1, start1), - -> [end1-1, end1),
exactly the historical bed6_line_tss().

The historical script wrote a missing gene_name as an empty string, whereas step 1
writes "NA". GENCODE v48 has no transcript without gene_name, so the two agree; the
producer refuses to run if an "NA" gene_name ever appears rather than guess.

Output: gencode.v48.tss.transcripts.bed.gz in input (GTF) order, deterministic
gzip (mtime 0, no name), staged then renamed. Genome-ordered copy: see
build_tss.sh, which applies the accepted step-3 sort command.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tempfile


def main():
    ap = argparse.ArgumentParser(description="transcript-level TSS BED from the accepted step-1 transcript BED")
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if "historical_reproduction" in a.transcripts:
        sys.exit("[error] refusing a historical-reproduction input")

    h = hashlib.sha256()
    n = {"+": 0, "-": 0}
    outdir = os.path.dirname(os.path.abspath(a.out))
    os.makedirs(outdir, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".stage_tss.", dir=outdir)
    os.close(fd)
    try:
        with gzip.open(a.transcripts, "rt") as fin, open(tmp, "wb") as raw, \
                gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            for ln, line in enumerate(fin, 1):
                f = line.rstrip("\n").split("\t")
                if len(f) != 6:
                    sys.exit(f"[error] transcript line {ln}: {len(f)} columns")
                chrom, s, e, name, _, strand = f
                s, e = int(s), int(e)
                p = name.split("|")
                if len(p) != 4:
                    sys.exit(f"[error] transcript line {ln}: name {name!r}")
                if p[2] == "NA":
                    sys.exit(f"[error] transcript line {ln}: gene_name NA is ambiguous (historically written as empty)")
                if strand == "+":
                    ts, te = s, s + 1
                elif strand == "-":
                    ts, te = e - 1, e
                else:
                    sys.exit(f"[error] transcript line {ln}: strand {strand!r}")
                row = f"{chrom}\t{ts}\t{te}\t{p[0]}|{p[1]}|{p[2]}\t0\t{strand}\n".encode()
                gz.write(row)
                h.update(row)
                n[strand] += 1
        os.chmod(tmp, 0o644)
        os.replace(tmp, a.out)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    print(f"[ok] {a.out}: {n['+'] + n['-']} rows ({n['+']} +, {n['-']} -), decompressed sha256 {h.hexdigest()}")


if __name__ == "__main__":
    main()
