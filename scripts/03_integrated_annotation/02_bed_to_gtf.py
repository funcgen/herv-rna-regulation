#!/usr/bin/env python3
import argparse, csv, re
from bisect import bisect_left
from typing import Dict, List, Optional, Tuple

# -----------------------
# Parsing helpers
# -----------------------

def parse_full_name(name: str) -> dict:
    parts = name.strip().split("|")
    meta = {
        "internal_subfam": parts[0],
        "locid": None, "has5": "0", "has3": "0",
        "ltr5": ".", "ltr3": "."
    }
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            meta[k] = v
    return meta

def infer_internal_family(internal_subfam: str) -> Optional[str]:
    s = internal_subfam.upper()
    if s.startswith("HERVK"): return "ERVK"
    if s.startswith("HERVH"): return "ERV1"
    if s.startswith("HERVW"): return "ERV1"
    if s.startswith("HERVE"): return "ERV1"
    if s.startswith("ERVL") or s.startswith("MER") or "MALR" in s: return "ERVL"
    if "ERVK" in s: return "ERVK"
    if "ERVL" in s: return "ERVL"
    if "ERV1" in s: return "ERV1"
    return None

def infer_family_from_classfamily(cf: str) -> Optional[str]:
    if not cf: return None
    toks = [t.strip() for t in cf.split("/")]
    for t in toks:
        if t.startswith("ERV"): return t
    for t in toks:
        if "ERV" in t: return t
    return None

def load_ltr_family_map(path: Optional[str]) -> Dict[str, str]:
    m = {}
    if not path: return m
    with open(path, "r") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        has_cf = "ClassFamily" in r.fieldnames
        for row in r:
            rep = (row.get("RepeatName") or "").strip()
            if not rep: continue
            if has_cf:
                cf = (row.get("ClassFamily") or "").strip()
            else:
                cf = (row.get("Class","").strip() + "/" + row.get("Family","").strip()).strip("/")
            fam = infer_family_from_classfamily(cf)
            if fam: m[rep] = fam
    return m

def parse_internal_name(name: str) -> str:
    m = re.search(r'([^|:]+-int)', name)
    return m.group(1) if m else name.split("|")[0]

def parse_ltr_subfamily(name: str) -> str:
    m = re.match(r'^(.+?)_merged', name)
    return m.group(1) if m else name.split("|")[0]

# -----------------------
# Interval indexing
# -----------------------

class BedRec:
    __slots__ = ("chrom","start","end","name","score","strand","extra")
    def __init__(self, chrom:str, start:int, end:int, name:str, score:str, strand:str, extra=None):
        self.chrom, self.start, self.end, self.name, self.score, self.strand = chrom, start, end, name, score, strand
        self.extra = extra or {}

def load_bed(path: str, kind="generic") -> List["BedRec"]:
    out = []
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            chrom, s, e, name, score, strand = parts[:6]
            rec = BedRec(chrom, int(s), int(e), name, score, strand)
            if kind == "full":
                rec.extra = parse_full_name(name)
            elif kind == "internal":
                rec.extra = {"internal_subfam": parse_internal_name(name)}
            elif kind == "ltr":
                rec.extra = {"ltr_subfam": parse_ltr_subfamily(name)}
            out.append(rec)
    return out

class ByChromStrand:
    def __init__(self, recs: List[BedRec]):
        self.idx: Dict[Tuple[str,str], Dict[str, List]] = {}
        for r in recs:
            key = (r.chrom, r.strand)
            bucket = self.idx.setdefault(key, {"starts": [], "recs": []})
            bucket["recs"].append(r)
        for key in self.idx:
            self.idx[key]["recs"].sort(key=lambda x: x.start)
            self.idx[key]["starts"] = [r.start for r in self.idx[key]["recs"]]

    def query_overlap(self, chrom: str, strand: str, qstart: int, qend: int) -> List[BedRec]:
        key = (chrom, strand)
        if key not in self.idx: return []
        starts = self.idx[key]["starts"]; recs = self.idx[key]["recs"]
        i = bisect_left(starts, qstart)
        i = max(0, i-16)
        out = []
        for r in recs[i:]:
            if r.start > qend: break
            if r.end > qstart and r.start < qend:
                out.append(r)
        return out

# -----------------------
# GTF helpers
# -----------------------

def write_gtf_exon_line(chrom, start0, end0, strand, source, score, attrs: dict) -> str:
    """
    Build a single GTF exon line string (no newline at end).
    BED 0-based half-open -> GTF 1-based inclusive.
    """
    start1 = start0 + 1
    end1   = end0
    attr_str = " ".join(f'{k} "{v}";' for k, v in attrs.items() if v is not None)
    return "\t".join([chrom, source, "exon", str(start1), str(end1), str(score), strand, ".", attr_str])

def chrom_numeric_key(chrom: str) -> int:
    """
    Map chromosome string to a numeric rank:
      chr1..chr22 -> 1..22
      chrX        -> 23
      chrY        -> 24
      chrM/chrMT  -> 25
      others      -> 1000+ (after canonical)
    Works also if 'chr' prefix is missing.
    """
    s = chrom
    if s.lower().startswith("chr"):
        s = s[3:]
    s_upper = s.upper()
    if s.isdigit():
        return int(s)
    if s_upper == "X":
        return 23
    if s_upper == "Y":
        return 24
    if s_upper in ("M", "MT"):
        return 25
    # non-canonical contigs (GL000..., KI..., etc.) go after
    return 1000

def gtf_sort_key(line: str):
    """
    Sorting key for GTF lines:
      1) numeric chromosome order (via chrom_numeric_key)
      2) start position
      3) chromosome name as tiebreaker
    """
    parts = line.split("\t")
    chrom = parts[0]
    start = int(parts[3])
    return (chrom_numeric_key(chrom), chrom, start)

# -----------------------
# Main
# -----------------------

def main(full_bed: str, internal_bed: str, ltr_bed: str, out_gtf: str,
         ltr_map: Optional[str], source: str):
    full = load_bed(full_bed, kind="full")
    internals = load_bed(internal_bed, kind="internal")
    ltrs = load_bed(ltr_bed, kind="ltr")

    idx_internal = ByChromStrand(internals)
    idx_ltr = ByChromStrand(ltrs)

    ltr_fam_map = load_ltr_family_map(ltr_map)

    # Collect all exon lines here
    rows: List[str] = []

    for F in full:
        meta = F.extra  # internal_subfam, locid, has5, has3, ltr5, ltr3
        locid = meta.get("locid") or f"{meta['internal_subfam']}_{F.chrom}:{F.start+1}-{F.end}{F.strand}"
        internal_subfam = meta["internal_subfam"]
        internal_fam = infer_internal_family(internal_subfam) or "ERV"
        score = 0

        # Common IDs (minimal set Stellarscope uses)
        common = {
            "gene_id": locid,
            "transcript_id": locid,
            "locus": locid,
            "locid": locid,
            "repClass": "LTR",
        }

        # ---------- FULL exon (counts the whole provirus)
        attrs_full = {
            **common,
            "repName": internal_subfam,
            "geneRegion": "full",
            "repFamily": internal_fam,
        }
        rows.append(
            write_gtf_exon_line(F.chrom, F.start, F.end, F.strand, source, score, attrs_full)
        )

        # ---------- INTERNAL exon (component)
        cand_int = idx_internal.query_overlap(F.chrom, F.strand, F.start, F.end)
        best_int = None; best_ov = -1
        for r in cand_int:
            if r.extra.get("internal_subfam") != internal_subfam:
                continue
            ov = min(F.end, r.end) - max(F.start, r.start)
            if ov > best_ov:
                best_ov = ov; best_int = r
        if best_int:
            attrs_internal = {
                **common,
                "repName": internal_subfam,
                "geneRegion": "internal",
                "repFamily": internal_fam,
            }
            rows.append(
                write_gtf_exon_line(best_int.chrom, best_int.start, best_int.end,
                                    best_int.strand, source, score, attrs_internal)
            )

        # ---------- LTR exons (components)
        want5 = meta.get("ltr5", ".")
        want3 = meta.get("ltr3", ".")
        cand_ltr = idx_ltr.query_overlap(F.chrom, F.strand, F.start, F.end)

        # group candidates by subfamily
        ltr_matches: Dict[str, List[BedRec]] = {}
        for r in cand_ltr:
            subfam = r.extra.get("ltr_subfam")
            if subfam:
                ltr_matches.setdefault(subfam, []).append(r)

        def pick_left_boundary(cands: List[BedRec]) -> Optional[BedRec]:
            """Closest to left boundary (F.start)."""
            if not cands:
                return None
            return min(cands, key=lambda r: abs(r.end - F.start))

        def pick_right_boundary(cands: List[BedRec]) -> Optional[BedRec]:
            """Closest to right boundary (F.end)."""
            if not cands:
                return None
            return min(cands, key=lambda r: abs(r.start - F.end))

        # Decide which boundary corresponds to 5'/3' based on strand
        if F.strand == "+":
            pick5 = pick_left_boundary
            pick3 = pick_right_boundary
        else:
            pick5 = pick_right_boundary
            pick3 = pick_left_boundary

        r5 = None
        r3 = None

        # Pick 5'
        if want5 and want5 != "." and want5 in ltr_matches:
            r5 = pick5(ltr_matches[want5])

        # Pick 3'
        if want3 and want3 != "." and want3 in ltr_matches:
            r3 = pick3(ltr_matches[want3])

        # ✅ If both sides collapse to the same exact interval, keep only one
        if r5 and r3 and (r5.chrom, r5.start, r5.end, r5.strand) == (r3.chrom, r3.start, r3.end, r3.strand):
            # if both intended but only one overlap exists, drop r3 (or drop r5, either is fine)
            r3 = None

        # Emit 5'
        if r5:
            repFam5 = ltr_fam_map.get(want5, internal_fam)
            attrs_ltr5 = {
                **common,
                "repName": want5,
                "geneRegion": "ltr",   # or "ltr5" if you want to keep that info
                "repFamily": repFam5,
            }
            rows.append(write_gtf_exon_line(r5.chrom, r5.start, r5.end,
                                            r5.strand, source, score, attrs_ltr5))

        # Emit 3'
        if r3:
            repFam3 = ltr_fam_map.get(want3, internal_fam)
            attrs_ltr3 = {
                **common,
                "repName": want3,
                "geneRegion": "ltr",   # or "ltr3"
                "repFamily": repFam3,
            }
            rows.append(write_gtf_exon_line(r3.chrom, r3.start, r3.end,
                                            r3.strand, source, score, attrs_ltr3))


    # ---- Sort all exon lines by numeric chromosome -> start
    rows.sort(key=gtf_sort_key)

    # ---- Write header + sorted lines
    with open(out_gtf, "w") as out:
        out.write("##gtf-version 2.2\n")
        out.write(f"# Generated with build_stellarscope_gtf_full_plus_components.py; source={source}\n")
        for line in rows:
            out.write(line + "\n")

    print(f"✅ Wrote GTF to {out_gtf}")
    print("Each locus (locid) has: 1 full-span exon (+ internal exon + up to 2 LTR exons when available).")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build Stellarscope GTF with full + component exons per provirus.")
    ap.add_argument("--full-bed", required=True, help="BED of merged full provirus spans (with locid, has5/has3, ltr5/ltr3 in name)")
    ap.add_argument("--internal-bed", required=True, help="BED of merged internal parts")
    ap.add_argument("--ltr-bed", required=True, help="BED of merged LTRs")
    ap.add_argument("--out-gtf", required=True, help="Output GTF path")
    ap.add_argument("--ltr-map", default=None, help="Optional ltr_subfamilies.tsv (RepeatName, ClassFamily) to set repFamily")
    ap.add_argument("--source", default="rmsk", help="GTF source field (default: rmsk)")
    args = ap.parse_args()
    main(args.full_bed, args.internal_bed, args.ltr_bed, args.out_gtf, args.ltr_map, args.source)
