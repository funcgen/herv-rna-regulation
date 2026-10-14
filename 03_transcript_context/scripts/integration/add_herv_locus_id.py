#!/usr/bin/env python3
"""
03_transcript_context step 10 (frozen global order 17, effective order 18):
attach HERV_id to the integrated locus table, producing the canonical shared
table transcript_context.with_herv_id.tsv (D010).

Refactor of scripts/integration/02_add_hervid.py (sha256 5a2a7654...). The
historical join rule is kept exactly:
  * mapping source: GTF records with geneRegion "internal"; HERV_id = gene_id
    (the GTF is SRC13 hervs_internal_and_ltrs.gtf, byte-identical to the
    historical hervs_full_ltr_int_v2.gtf; D019)
  * key: chrom, start, end, strand, with table start (0-based) + 1 == GTF start
    (1-based) and table end == GTF end
  * attributes parsed with the historical regex  (\\S+)\\s+"([^"]+)"
  * a coordinate key appearing twice among internal records is an error
  * HERV_id becomes the FIRST column; row order is the table's
What changed:
  * text-level, not pandas: each output row is HERV_id, a TAB, and the input
    line verbatim, so the 41 existing columns cannot be altered by parsing or
    re-serialization (removing column 1 gives back the input byte for byte)
  * an unmapped locus, or a HERV_id used twice, is an error (historically
    unmapped rows were only warned about and written with an empty HERV_id)
  * the input table and GTF must be the accepted artifacts (digests below)
"""
import argparse
import hashlib
import os
import re
import sys
import tempfile

N_LOCI = 53176
TABLE_SHA256 = "b38df4e9c3b84e58f393f58c4440545fb96415cb74e118fa878fdf3fc994e74a"  # effective order 17
GTF_MD5 = "3570d4709785a2890dd9bd86b8783d91"  # SRC13, 01_sources/derived/hervs_internal_and_ltrs.gtf
_attr_re = re.compile(r'(\S+)\s+"([^"]+)"')  # historical parse_gtf_attributes()


def die(msg):
    sys.exit(f"[error] {msg}")


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def internal_map(gtf):
    """(chrom, start1, end1, strand) -> gene_id over geneRegion 'internal' records."""
    m, n = {}, 0
    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ad = {k: v for k, v in _attr_re.findall(p[8])}
            if ad.get("geneRegion") != "internal" or ad.get("gene_id") is None:
                continue
            key = (p[0], int(p[3]), int(p[4]), p[6])
            if key in m:
                die(f"duplicate coordinate key among internal GTF records: {key}")
            m[key] = ad["gene_id"]
            n += 1
    return m, n


def main():
    ap = argparse.ArgumentParser(description="attach HERV_id (effective order 18)")
    ap.add_argument("--table", required=True, help="results/integration/transcript_context.tsv (order 17)")
    ap.add_argument("--gtf", required=True, help="SRC13 01_sources/derived/hervs_internal_and_ltrs.gtf")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if digest(a.table, "sha256") != TABLE_SHA256:
        die("--table is not the accepted effective-order-17 transcript_context.tsv")
    if digest(a.gtf, "md5") != GTF_MD5:
        die("--gtf is not the accepted SRC13 hervs_internal_and_ltrs.gtf")
    print(f"[ok] table sha256 {TABLE_SHA256}; GTF md5 {GTF_MD5}", file=sys.stderr)

    gmap, n_internal = internal_map(a.gtf)
    print(f"[ok] GTF: {n_internal} internal records, {len(gmap)} unique coordinate keys", file=sys.stderr)

    with open(a.table, newline="") as fh:
        lines = fh.readlines()
    header = lines[0].rstrip("\n").split("\t")
    if "HERV_id" in header:
        die("input already carries HERV_id")
    ci = {c: header.index(c) for c in ("chrom", "start", "end", "strand", "locid")}
    rows = lines[1:]
    if len(rows) != N_LOCI:
        die(f"{len(rows)} data rows, expected {N_LOCI}")

    out_lines, ids, locids, missing = ["HERV_id\t" + lines[0]], set(), set(), []
    for ln in rows:
        f = ln.rstrip("\n").split("\t")
        key = (f[ci["chrom"]], int(f[ci["start"]]) + 1, int(f[ci["end"]]), f[ci["strand"]])
        hid = gmap.get(key)
        if hid is None:
            missing.append(f[ci["locid"]])
            continue
        if hid in ids:
            die(f"HERV_id {hid} would be attached to more than one locus")
        ids.add(hid)
        locids.add(f[ci["locid"]])
        out_lines.append(hid + "\t" + ln)
    if missing:
        die(f"{len(missing)} loci have no internal GTF record, e.g. {missing[:5]}")
    if len(locids) != N_LOCI or len(out_lines) != N_LOCI + 1:
        die("one-row-per-locus invariant broken")
    print(f"[ok] mapped {len(ids)}/{N_LOCI} loci; {n_internal - len(ids)} internal GTF records unused",
          file=sys.stderr)

    outdir = os.path.dirname(os.path.abspath(a.out))
    os.makedirs(outdir, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".stage_with_herv_id.", dir=outdir)
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            fh.writelines(out_lines)
        os.chmod(tmp, 0o644)
        os.replace(tmp, a.out)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    print(f"[OK] Wrote: {a.out} ({N_LOCI} rows, {len(header) + 1} cols)", file=sys.stderr)


if __name__ == "__main__":
    main()
