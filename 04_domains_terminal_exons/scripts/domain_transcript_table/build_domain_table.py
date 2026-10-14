#!/usr/bin/env python3
"""
04_domains_terminal_exons step 1 (frozen global order 29, effective order 30):
the domain-centred transcript table.

Refactor of scripts/integration/05_build_domain_centered_terminal_exon_table.py
(sha256 7b683741...). One run per parameter profile (config/parameters.yaml, domains):

  conserved     coverage >= domains.min_hmm_profile_coverage (0.40)   -> VAL31, VAL32, VAL33
  all_coverage  coverage >= domains.all_coverage_threshold (0.0)      -> conservation
                comparisons (effective orders 32-33; VAL35 188 + 237 = 425 domains)
  both          exonic overlap fraction >= domains.min_exonic_overlap_fraction (0.80);
                terminal exon under the stringent rule when the terminal-exon overlap
                fraction >= domains.min_terminal_exon_overlap_fraction (0.80)

Semantics (historical, unchanged)
  domain        one exact SRC05 BED record; domain_uid = "DOM_" + sha1("chrom:start-end:
                strand|bed_name")[:16]; metadata (DomainType, Coverage, Subfamily, Locus,
                Name, Conservedness) from SRC06, joined on Query + Domain, one to one;
                SRC06 Coverage is authoritative and must equal the BED-name coverage
  association   one row per unique domain x transcript with >= 1 bp of SAME-STRAND exon
                overlap; exonic_overlap_bp = length of the UNION of the overlap intervals
                over all exons of that transcript (never a raw sum); fraction = that union
                / domain length (BED end - start); terminal overlap likewise over the
                transcript's terminal exon
  terminal exon the exon furthest 3' in transcript orientation
  domain table  one row per qualifying domain; terminal_all / isoform_dependent /
                nonterminal_all over its qualifying transcripts, and the same over
                multi-exon transcripts only

What changed: the shared 03_transcript_context layer is consumed instead of recomputed
(MODULE_DEPENDENCIES: "domain-vs-feature overlaps; terminal-exon definitions ... must not
be recomputed inside 04"): the order-14 same-strand domain x transcript-exon -wo table
supplies the overlaps, the step-4 last_exon_per_tx supplies the terminal exon, and the
step-3 exon BED supplies per-transcript exon counts and annotations. Thresholds come
from parameters.yaml, not from code. Every input is checked against its accepted digest.

Outputs, in <outdir>/<profile>/ (historical file names):
  domain_centered_qualifying_domain_transcript_associations.tsv
  domain_centered_unique_domains.tsv
  domain_centered_counts.tsv
  domain_centered_qc.tsv
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tempfile
from collections import defaultdict

import re

import pandas as pd

BED_COLUMNS = ["chrom", "start", "end", "bed_name", "score", "strand"]
PINNED = {  # accepted identities
    "domains_bed": ("md5", "f02b8bb89f19a0a2a3ac45efcd900898"),   # SRC05
    "domains_tsv": ("md5", "162cba73027d1e4500dc1dbb2671c84a"),   # SRC06
}
WO_HEADER = "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\toverlap_bp"
OUT_FILES = ("domain_centered_qualifying_domain_transcript_associations.tsv", "domain_centered_unique_domains.tsv",
             "domain_centered_counts.tsv", "domain_centered_qc.tsv")


DOMAIN_KEYS = ("min_hmm_profile_coverage", "min_exonic_overlap_fraction", "all_coverage_threshold",
               "min_terminal_exon_overlap_fraction")


def read_domain_params(path):
    """The `domains:` block of config/parameters.yaml, without a YAML dependency (the CNAG
    environment has none; same convention as 02_rbp read_config.py). A missing key is fatal."""
    text = open(path, encoding="utf-8").read()
    m = re.search(r"^domains:\s*\n(.*?)(?=^\S|\Z)", text, re.M | re.S)
    if not m:
        die("no domains: block in parameters.yaml")
    vals = {}
    for key in DOMAIN_KEYS:
        k = re.search(r"^\s+" + re.escape(key) + r":\s*([0-9.eE+-]+)\s*(?:#.*)?$", m.group(1), re.M)
        if not k:
            die(f"parameter domains.{key} not found")
        vals[key] = float(k.group(1))
    return vals


def die(msg):
    sys.exit(f"[error] {msg}")


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def verify_in_manifest(path):
    """03_transcript_context assets: verify against the checksums.sha256 of their results directory."""
    p = os.path.abspath(path)
    d = os.path.dirname(p)
    while not os.path.exists(os.path.join(d, "checksums.sha256")):
        if d == os.path.dirname(d):
            die(f"no checksums.sha256 above {path}")
        d = os.path.dirname(d)
    rel = "./" + os.path.relpath(p, d)
    want = None
    for line in open(os.path.join(d, "checksums.sha256")):
        h, f = line.rstrip("\n").split(None, 1)
        if f in (rel, rel[2:]):
            want = h
    if want is None:
        die(f"{rel} not listed in {d}/checksums.sha256")
    if digest(p, "sha256") != want:
        die(f"{path} does not match its accepted digest")
    print(f"[ok] {rel} sha256 {want} (manifest)", file=sys.stderr)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


# ---------------------------------------------------------------- domains (historical)
def domain_uid(chrom, start, end, strand, bed_name):
    return "DOM_" + hashlib.sha1(f"{chrom}:{start}-{end}:{strand}|{bed_name}".encode("utf-8")).hexdigest()[:16]


def load_domains(bed_path, tsv_path):
    d = pd.read_csv(bed_path, sep="\t", header=None, names=BED_COLUMNS, dtype={"chrom": str, "bed_name": str, "strand": str})
    parsed = d["bed_name"].str.rsplit("|", n=4, expand=True)
    if parsed.shape[1] != 5:
        die("domain BED names could not be parsed with rsplit('|', 4)")
    parsed.columns = ["Query", "Domain", "Gene", "Conservedness", "Coverage_bed"]
    d = pd.concat([d, parsed], axis=1)
    d["start"] = pd.to_numeric(d["start"], errors="raise").astype(int)
    d["end"] = pd.to_numeric(d["end"], errors="raise").astype(int)
    d["domain_length_bp"] = d["end"] - d["start"]
    d["Coverage_bed"] = pd.to_numeric(d["Coverage_bed"], errors="coerce")
    if (d["domain_length_bp"] <= 0).any():
        die("zero- or negative-length domain intervals")
    d["domain_uid"] = [domain_uid(c, s, e, st, n) for c, s, e, st, n in
                       zip(d["chrom"], d["start"], d["end"], d["strand"], d["bed_name"])]
    if d["domain_uid"].duplicated().any():
        die("non-unique exact domain BED records")
    n_bed = len(d)
    meta = pd.read_csv(tsv_path, sep="\t", dtype=str)
    if meta.duplicated(["Query", "Domain"]).any():
        die("Query + Domain is not unique in the domain metadata table")
    meta["Coverage_metadata"] = pd.to_numeric(meta["Coverage"], errors="coerce")
    keep = [x for x in ["Query", "Domain", "DomainType", "Coverage_metadata", "Conservedness", "Subfamily", "Locus",
                        "Name", "Gene"] if x in meta.columns]
    d = d.merge(meta[keep], on=["Query", "Domain"], how="left", suffixes=("_bed", "_metadata"), validate="one_to_one")
    if d["Coverage_metadata"].isna().any():
        die(f"{int(d['Coverage_metadata'].isna().sum())} BED domain rows have no metadata")
    d["Coverage"] = d["Coverage_metadata"]
    if ((d["Coverage_bed"] - d["Coverage"]).abs() > 1e-9).any():
        die("coverage differs between the domain BED and the metadata table")
    return d, n_bed


# ---------------------------------------------------------------- shared 03 assets
def load_transcripts(exons_path, last_exons_path):
    tx = {}
    with open_text(exons_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            p = f[3].split("|")
            if len(p) < 6:
                die(f"malformed exon name {f[3]}")
            rec = tx.get(p[0])
            if rec is None:
                tx[p[0]] = {"gene_id": p[1], "gene_name": p[2], "gene_type": p[4], "tx_type": p[5], "n_exons": 1,
                            "strand": f[5]}
            else:
                rec["n_exons"] += 1
    with open_text(last_exons_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            t = f[3].split("|", 1)[0]
            if t not in tx or "last" in tx[t]:
                die(f"last-exon BED: transcript {t} unknown or repeated")
            tx[t]["last"] = (int(f[1]), int(f[2]))
    missing = [t for t, r in tx.items() if "last" not in r]
    if missing:
        die(f"{len(missing)} transcripts have no terminal exon in the step-4 asset")
    return tx


def association_intervals(wo_path, domains_by_key, tx):
    """Same-strand domain x exon overlaps from the order-14 -wo table, grouped by domain x transcript."""
    iv = defaultdict(lambda: {"all": [], "terminal": []})
    n_rows = 0
    with open(wo_path) as fh:
        if fh.readline().rstrip("\n") != WO_HEADER:
            die("domain x exon overlap table header is not the order-14 -wo header")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a_s, a_e, b_s, b_e, ov = int(f[1]), int(f[2]), int(f[7]), int(f[8]), int(f[12])
            if f[5] != f[11] or f[0] != f[6]:
                die("overlap table row is not same-strand")
            s, e = max(a_s, b_s), min(a_e, b_e)
            if e - s != ov or ov <= 0:
                die(f"overlap_bp {ov} disagrees with the intervals ({s}, {e})")
            n_rows += 1
            dom = domains_by_key.get((f[0], a_s, a_e, f[5], f[3]))
            if dom is None:
                continue  # below the coverage threshold of this profile
            t = f[9].split("|", 1)[0]
            rec = tx[t]
            last = rec["last"]
            terminal = (b_s, b_e) == last
            iv[(dom, t)]["all"].append((s, e))
            if terminal:
                iv[(dom, t)]["terminal"].append((s, e))
    return iv, n_rows


def union_length(intervals):
    total, cur = 0, None
    for s, e in sorted(intervals):
        if cur is None or s > cur[1]:
            if cur:
                total += cur[1] - cur[0]
            cur = [s, e]
        else:
            cur[1] = max(cur[1], e)
    return total + (cur[1] - cur[0] if cur else 0)


# ---------------------------------------------------------------- tables (historical)
def build_associations(domains, iv, tx, terminal_rule):
    lookup = domains.set_index("domain_uid").to_dict("index")
    rows = []
    for (uid, t), sets in iv.items():
        dom, rec = lookup[uid], tx[t]
        rows.append({
            "domain_uid": uid, "tx_id": t, "gene_id": rec["gene_id"], "gene_name": rec["gene_name"],
            "gene_type": rec["gene_type"], "tx_type": rec["tx_type"], "n_exons": int(rec["n_exons"]),
            "is_single_exon_transcript": bool(rec["n_exons"] == 1),
            "exonic_overlap_bp": union_length(sets["all"]), "terminal_exon_overlap_bp": union_length(sets["terminal"]),
            "overlaps_terminal_exon": bool(union_length(sets["terminal"]) > 0),
            "chrom": dom["chrom"], "start": int(dom["start"]), "end": int(dom["end"]), "strand": dom["strand"],
            "domain_length_bp": int(dom["domain_length_bp"]), "bed_name": dom["bed_name"], "Query": dom["Query"],
            "Domain": dom["Domain"], "DomainType": dom["DomainType"], "Coverage": float(dom["Coverage"]),
            "Conservedness": dom.get("Conservedness_metadata", dom.get("Conservedness_bed", ".")),
            "Subfamily": dom["Subfamily"], "Locus": dom["Locus"], "Name": dom.get("Name", "."),
        })
    a = pd.DataFrame(rows)
    a["exonic_overlap_fraction"] = a["exonic_overlap_bp"] / a["domain_length_bp"]
    a["terminal_exon_overlap_fraction"] = a["terminal_exon_overlap_bp"] / a["domain_length_bp"]
    a["terminal_under_stringent_rule"] = a["terminal_exon_overlap_fraction"] >= terminal_rule
    return a


def classify(flags, kind_prefix, suffix=""):
    if flags.all():
        return f"terminal_all{suffix}"
    if flags.any():
        return f"isoform_dependent{suffix}"
    return f"nonterminal_all{suffix}"


def summarise_domains(q):
    rows = []
    for uid, g in q.groupby("domain_uid", sort=False):
        first = g.iloc[0]
        tf = g["terminal_under_stringent_rule"].astype(bool)
        multi = g[~g["is_single_exon_transcript"]]
        rows.append({
            "domain_uid": uid, "chrom": first["chrom"], "start": int(first["start"]), "end": int(first["end"]),
            "strand": first["strand"], "domain_length_bp": int(first["domain_length_bp"]), "Query": first["Query"],
            "Domain": first["Domain"], "DomainType": first["DomainType"], "Coverage": float(first["Coverage"]),
            "Conservedness": first["Conservedness"], "Subfamily": first["Subfamily"], "Locus": first["Locus"],
            "n_qualifying_transcripts": int(g["tx_id"].nunique()), "n_qualifying_genes": int(g["gene_id"].nunique()),
            "n_terminal_transcripts": int(tf.sum()), "n_nonterminal_transcripts": int((~tf).sum()),
            "n_single_exon_transcripts": int(g["is_single_exon_transcript"].sum()),
            "n_multi_exon_transcripts": int((~g["is_single_exon_transcript"]).sum()),
            "domain_terminal_classification": classify(tf, ""),
            "multi_exon_classification": "single_exon_only" if multi.empty
            else classify(multi["terminal_under_stringent_rule"].astype(bool), "", "_multi_exon"),
            "terminal_in_any_qualifying_transcript": bool(tf.any()),
            "terminal_in_all_qualifying_transcripts": bool(tf.all()),
        })
    return pd.DataFrame(rows)


def build_counts(ds, q):
    nan = float("nan")
    total = len(ds)
    c = ds["domain_terminal_classification"].value_counts().to_dict()
    m = ds[ds["multi_exon_classification"] != "single_exon_only"]
    mt = len(m)
    mc = m["multi_exon_classification"].value_counts().to_dict()
    pct = lambda n, d: 100 * n / d if d else nan
    any_multi = int((m["multi_exon_classification"] != "nonterminal_all_multi_exon").sum())
    rows = [
        ("unique_qualifying_domains", total, 100.0 if total else nan),
        ("terminal_all", c.get("terminal_all", 0), pct(c.get("terminal_all", 0), total)),
        ("isoform_dependent", c.get("isoform_dependent", 0), pct(c.get("isoform_dependent", 0), total)),
        ("nonterminal_all", c.get("nonterminal_all", 0), pct(c.get("nonterminal_all", 0), total)),
        ("terminal_in_any_transcript", int(ds["terminal_in_any_qualifying_transcript"].sum()),
         100 * ds["terminal_in_any_qualifying_transcript"].mean() if total else nan),
        ("domains_with_multi_exon_transcript", mt, pct(mt, total)),
        ("terminal_all_multi_exon", mc.get("terminal_all_multi_exon", 0), pct(mc.get("terminal_all_multi_exon", 0), mt)),
        ("isoform_dependent_multi_exon", mc.get("isoform_dependent_multi_exon", 0), pct(mc.get("isoform_dependent_multi_exon", 0), mt)),
        ("nonterminal_all_multi_exon", mc.get("nonterminal_all_multi_exon", 0), pct(mc.get("nonterminal_all_multi_exon", 0), mt)),
        ("terminal_in_any_multi_exon_transcript", any_multi,
         100 * (m["multi_exon_classification"] != "nonterminal_all_multi_exon").mean() if mt else nan),
        ("qualifying_domain_transcript_associations", len(q), 100.0),
        ("unique_herv_loci", ds["Locus"].nunique(), nan),
        ("unique_genes", q["gene_id"].nunique(), nan),
    ]
    return pd.DataFrame(rows, columns=["metric", "n", "percent"])


def run(domains_bed, domains_tsv, wo_path, exons_path, last_exons_path, min_cov, min_ov, term_rule):
    """The whole computation for one profile (no input guards; main() applies them)."""
    domains, n_bed = load_domains(domains_bed, domains_tsv)
    domains = domains[domains["Coverage"] >= min_cov].copy()
    by_key = {(c, s, e, st, n): u for c, s, e, st, n, u in
              zip(domains["chrom"], domains["start"], domains["end"], domains["strand"], domains["bed_name"], domains["domain_uid"])}
    tx = load_transcripts(exons_path, last_exons_path)
    iv, n_wo = association_intervals(wo_path, by_key, tx)
    assoc_all = build_associations(domains, iv, tx, term_rule)
    if assoc_all.empty:
        die("no same-strand domain-exon overlaps")
    q = assoc_all[assoc_all["exonic_overlap_fraction"] >= min_ov].copy()
    ds = summarise_domains(q)
    counts = build_counts(ds, q)
    qc = pd.DataFrame({
        "metric": ["input_domain_bed_rows", "domains_after_coverage_filter", "all_same_strand_domain_transcript_overlaps",
                   "qualifying_associations", "qualifying_unique_domains",
                   "terminal_associations_with_terminal_fraction_below_0.8"],
        "value": [n_bed, len(domains), len(assoc_all), len(q), len(ds),
                  int(((q["overlaps_terminal_exon"]) & (~q["terminal_under_stringent_rule"])).sum())]})
    return q, ds, counts, qc, (n_wo, n_bed, len(domains), assoc_all)


def main():
    ap = argparse.ArgumentParser(description="domain-centred transcript table (effective order 30)")
    ap.add_argument("--domains-bed", required=True, help="SRC05 01_sources/frozen/herv/HERV_domains_v5.bed")
    ap.add_argument("--domains-tsv", required=True, help="SRC06 01_sources/frozen/herv/HERV_loci_annotated_domains_v5.tsv")
    ap.add_argument("--domain-exon-overlaps", required=True, help="03 overlap/domains_vs_gencode_tx_exons.same_strand.wo.tsv")
    ap.add_argument("--tx-exons-bed", required=True, help="03 features/sorted/gencode.v48.exons.sorted.bed.gz")
    ap.add_argument("--last-exons-bed", required=True, help="03 features/sorted/gencode.v48.last_exon_per_tx.sorted.bed.gz")
    ap.add_argument("--params", required=True, help="config/parameters.yaml")
    ap.add_argument("--profile", required=True, choices=["conserved", "all_coverage"])
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    P = read_domain_params(a.params)
    min_cov = P["min_hmm_profile_coverage"] if a.profile == "conserved" else P["all_coverage_threshold"]
    min_ov, term_rule = P["min_exonic_overlap_fraction"], P["min_terminal_exon_overlap_fraction"]
    print(f"[params] profile={a.profile} coverage>={min_cov} exonic_overlap>={min_ov} terminal_rule>={term_rule}", file=sys.stderr)

    for key in PINNED:
        algo, want = PINNED[key]
        if digest(getattr(a, key), algo) != want:
            die(f"{key} is not the accepted source ({algo} {want})")
    for p in (a.domain_exon_overlaps, a.tx_exons_bed, a.last_exons_bed):
        verify_in_manifest(p)
    if "gene_exons_merged" in a.tx_exons_bed or "gene_introns" in a.tx_exons_bed:
        die("transcript-resolved exons required; gene-level merged assets (D022) are not an input")

    q, ds, counts, qc, info = run(a.domains_bed, a.domains_tsv, a.domain_exon_overlaps, a.tx_exons_bed,
                                  a.last_exons_bed, min_cov, min_ov, term_rule)
    n_wo, n_bed, n_dom, assoc_all = info

    out = os.path.join(a.outdir, a.profile)
    os.makedirs(a.outdir, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=f".stage_{a.profile}.", dir=a.outdir)
    q.sort_values(["chrom", "start", "end", "DomainType", "tx_id"], kind="mergesort").to_csv(os.path.join(stage, OUT_FILES[0]), sep="\t", index=False)
    ds.sort_values(["chrom", "start", "end", "DomainType"], kind="mergesort").to_csv(os.path.join(stage, OUT_FILES[1]), sep="\t", index=False)
    counts.to_csv(os.path.join(stage, OUT_FILES[2]), sep="\t", index=False)
    qc.to_csv(os.path.join(stage, OUT_FILES[3]), sep="\t", index=False)
    for f in OUT_FILES:
        os.chmod(os.path.join(stage, f), 0o644)
    if os.path.exists(out):
        for f in os.listdir(out):
            os.remove(os.path.join(out, f))
        os.rmdir(out)
    os.rename(stage, out)
    os.chmod(out, 0o755)
    print(f"[qc] -wo rows {n_wo}; domains {n_bed} -> {n_dom}; domain x transcript overlaps {len(assoc_all)}; "
          f"qualifying {len(q)} associations, {len(ds)} domains", file=sys.stderr)
    print(counts.to_string(index=False), file=sys.stderr)
    print(f"[OK] wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
