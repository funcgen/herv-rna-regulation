#!/usr/bin/env python3
"""
07_sparcs step 3 (frozen global order 44), executed as effective order 16 under
erratum D017: the minimal LTR feature projection.

Replaces, for the integrated HERV-locus table, the unpublished 173 MB
LTR_fully_annotated.tsv (SRC11, DO_NOT_MIGRATE, decision D014). The historical
integration (extract_ltr_metrics in 01_build_integrated_internal_table_simplify.py)
reads exactly four of its 35 columns; this projection rebuilds only those four,
from three published, checksummed files:

  column          definition                                   source
  sequence_name   LTR identifier (BED name)                    SRC07 HERV_ltr.bed, column 4
  LTR_length      end - start (BED, 0-based half-open)         SRC07 HERV_ltr.bed
  n_motifs        number of rows of the published LTR FIMO     SRC09 merged_fimo_parsed_v4.tsv
                  table whose sequence_name is this LTR;       (no threshold applied: the
                  0 when it has none                            published table is used as is)
  status          'status' of the published U3-R-U5            SRC08 HERV_LTR_U3_R_U5_catalogue.tsv
                  reconstruction for this LTR                   (column 'status', joined on 'name')

One row per SRC07 LTR, in SRC07 order. Every LTR must be in SRC08 exactly once;
every FIMO sequence_name must be an SRC07 LTR; duplicate keys are errors. No other
field is carried: integration derives the has_ltr / confidence columns itself, and
SRC10 (the internal-to-LTR map) is read by integration, not here.

The FIMO member (6.27 GB) is streamed once; only a per-LTR counter is kept.

Output: gzip TSV (mtime 0, no name), header line, staged then renamed.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tempfile

EXPECTED_MD5 = {  # accepted 01_sources identities
    "ltr_bed": "ef318953ea35e2821431f943c5e20007",
    "catalogue": "8c685d0dbb30372cbab5e10167f308fb",
    "fimo": "0a86e40e8a6648faad2b214079b1d255",
}
COLUMNS = ["sequence_name", "LTR_length", "n_motifs", "status"]


def md5(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def die(msg):
    sys.exit(f"[error] {msg}")


def main():
    ap = argparse.ArgumentParser(description="minimal LTR feature projection (D017)")
    ap.add_argument("--ltr-bed", required=True, help="SRC07 HERV_ltr.bed")
    ap.add_argument("--catalogue", required=True, help="SRC08 HERV_LTR_U3_R_U5_catalogue.tsv")
    ap.add_argument("--fimo", required=True, help="SRC09 member merged_fimo_parsed_v4.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    for key, path in (("ltr_bed", a.ltr_bed), ("catalogue", a.catalogue), ("fimo", a.fimo)):
        if "LTR_fully_annotated" in os.path.basename(path):
            die("LTR_fully_annotated.tsv (SRC11) is excluded from the release (D014) and must not be an input")
        got = md5(path)
        if got != EXPECTED_MD5[key]:
            die(f"{key} md5 {got} != accepted {EXPECTED_MD5[key]}")
        print(f"[ok] {key}: {path} md5 {got}")

    # SRC07: universe, order, length
    order, length = [], {}
    with open(a.ltr_bed) as fh:
        for ln, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                die(f"LTR BED line {ln}: {len(f)} columns")
            name, s, e = f[3], int(f[1]), int(f[2])
            if name in length:
                die(f"duplicate LTR name {name}")
            if e <= s:
                die(f"LTR {name}: non-positive length")
            order.append(name)
            length[name] = e - s

    # SRC08: status by name
    status = {}
    with open(a.catalogue) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i_name, i_status = hdr.index("name"), hdr.index("status")
        for ln, line in enumerate(fh, 2):
            f = line.rstrip("\n").split("\t")
            if f[i_name] in status:
                die(f"catalogue: duplicate name {f[i_name]}")
            status[f[i_name]] = f[i_status]
    missing = [n for n in order if n not in status]
    extra = set(status) - set(length)
    if missing or extra:
        die(f"catalogue/BED key mismatch: {len(missing)} BED LTRs not in catalogue, {len(extra)} catalogue names not in BED")

    # SRC09: stream once, count rows per sequence_name
    counts = dict.fromkeys(order, 0)
    n_rows = n_headers = 0
    with open(a.fimo, "rb") as fh:
        hdr = fh.readline().rstrip(b"\n").split(b"\t")
        i_seq = hdr.index(b"sequence_name")
        for line in fh:
            key = line.split(b"\t", i_seq + 1)[i_seq].decode()
            if key == "sequence_name":
                n_headers += 1
                continue
            try:
                counts[key] += 1
            except KeyError:
                die(f"FIMO sequence_name {key!r} is not an SRC07 LTR")
            n_rows += 1
    if n_headers:
        die(f"{n_headers} repeated header lines inside the FIMO table")
    print(f"[ok] FIMO: {n_rows} data rows, {sum(1 for v in counts.values() if v)} LTRs with >= 1 hit")

    outdir = os.path.dirname(os.path.abspath(a.out))
    os.makedirs(outdir, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".stage_ltr_projection.", dir=outdir)
    os.close(fd)
    h = hashlib.sha256()
    try:
        with open(tmp, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            for row in ["\t".join(COLUMNS)] + [f"{n}\t{length[n]}\t{counts[n]}\t{status[n]}" for n in order]:
                b = (row + "\n").encode()
                gz.write(b); h.update(b)
        os.chmod(tmp, 0o644)
        os.replace(tmp, a.out)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    print(f"[ok] {a.out}: {len(order)} LTRs x {len(COLUMNS)} columns; sum(n_motifs) = {n_rows}; "
          f"decompressed sha256 {h.hexdigest()}")


if __name__ == "__main__":
    main()
