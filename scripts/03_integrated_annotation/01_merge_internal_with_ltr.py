#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict
from bisect import bisect_left, bisect_right

import pandas as pd
from tqdm import tqdm

# ==============================
# Helpers
# ==============================

def extract_internal_subfamily(name: str) -> str:
    """From 'MLT1K-int_pos_chr1_...' extract 'MLT1K-int' base token (keeps '-int')."""
    m = re.search(r'([^|:]+-int)', str(name))
    return m.group(1) if m else ""

def extract_ltr_subfamily(name: str) -> str:
    """From 'LTR5_Hs_merged_pos_chr1_...' extract 'LTR5_Hs'."""
    m = re.match(r'^(.+?)_merged', str(name))   # non-greedy up to '_merged'
    return m.group(1) if m else name.split("|")[0]

def build_classfam_maps(df_map: pd.DataFrame):
    """
    Build:
      - internal_to_family: internal base (w/o '-int') -> set of allowed class/family strings
      - ltr_to_family:      ltr subfamily              -> single class/family string (string as-is)
    """
    internal_to_family = defaultdict(set)
    ltr_to_family = {}

    for _, row in df_map.iterrows():
        rep = str(row["RepeatName"])
        cf = str(row["ClassFamily"])
        if rep.endswith("-int"):
            internal_to_family[rep[:-4]].add(cf)  # remove "-int"
        else:
            # if duplicate mappings exist, last one wins (OK, they should be consistent)
            ltr_to_family[rep] = cf
    return internal_to_family, ltr_to_family

def build_group_index(df):
    """Index a BED-like DataFrame per (chrom, strand) with sorted arrays for fast window scans."""
    groups = {}
    for (chrom, strand), sub in df.groupby(["chrom", "strand"], sort=False):
        sub = sub.sort_values(["start", "end"], kind="mergesort").reset_index(drop=True)
        groups[(chrom, strand)] = {
            "df": sub,
            "starts": sub["start"].to_numpy(),
            "ends": sub["end"].to_numpy(),
        }
    return groups

def choose_upstream(up_hits: pd.DataFrame, istart: int, strategy: str):
    """Return (new_start, chosen_row or None) for upstream side."""
    if up_hits.empty:
        return None, None
    H = up_hits.copy()
    if strategy == "nearest":
        H["dist"] = istart - H["end"]
        H["len_neg"] = -(H["end"] - H["start"])
        H = H.sort_values(["dist", "len_neg"], kind="mergesort")
        r = H.iloc[0]
        return int(r["start"]), r
    else:  # furthest within-gap
        idx = H["start"].idxmin()
        r = H.loc[idx]
        return int(r["start"]), r

def choose_downstream(dn_hits: pd.DataFrame, iend: int, strategy: str):
    """Return (new_end, chosen_row or None) for downstream side."""
    if dn_hits.empty:
        return None, None
    H = dn_hits.copy()
    if strategy == "nearest":
        H["dist"] = H["start"] - iend
        H["len_neg"] = -(H["end"] - H["start"])
        H = H.sort_values(["dist", "len_neg"], kind="mergesort")
        r = H.iloc[0]
        return int(r["end"]), r
    else:  # furthest within-gap
        idx = H["end"].idxmax()
        r = H.loc[idx]
        return int(r["end"]), r

# ==============================
# Core
# ==============================

def main(
    internal_bed_path: str,
    ltr_bed_path: str,
    ltr_dict_path: str,
    output_bed_path: str,
    ltr_gap: int,
    strategy: str,
    map_out: str = None,
):
    # -------- Load internal BED
    internal_df = pd.read_csv(
        internal_bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
        dtype={"chrom": str, "start": int, "end": int, "name": str, "score": str, "strand": str},
        engine="c",
    )
    # keep internal base with '-int' (e.g., 'HERVK13-int')
    internal_df["internal_base"] = internal_df["name"].map(extract_internal_subfamily)
    # also store subfamily without '-int'
    internal_df["subfamily"] = internal_df["internal_base"].str.replace(r"-int$", "", regex=True)

    # -------- Load LTR dictionary
    ltr_dict_df = pd.read_csv(ltr_dict_path, sep="\t")
    if not {"RepeatName", "ClassFamily"}.issubset(ltr_dict_df.columns):
        raise ValueError("`ltr_subfamilies.tsv` must contain columns: RepeatName, ClassFamily")
    internal_to_family, ltr_to_family = build_classfam_maps(ltr_dict_df)

    # -------- Load merged LTR BED
    ltr_df = pd.read_csv(
        ltr_bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
        dtype={"chrom": str, "start": int, "end": int, "name": str, "score": str, "strand": str},
        engine="c",
    )
    ltr_df["ltr_subfamily"] = ltr_df["name"].map(extract_ltr_subfamily)
    ltr_df["classfam"] = ltr_df["ltr_subfamily"].map(lambda s: ltr_to_family.get(s, ""))

    # Optionally restrict to mapped LTRs
    ltr_df = ltr_df[ltr_df["classfam"].astype(bool)].copy()

    # -------- Index LTRs per (chrom,strand)
    ltr_groups = build_group_index(ltr_df)

    # -------- Propose extensions and capture chosen LTRs
    proposed = []
    it = tqdm(internal_df.itertuples(index=False), total=len(internal_df), desc="Proposing extensions")
    for row in it:
        chrom, istart, iend, name, score, strand, internal_base, subf = (
            row.chrom, int(row.start), int(row.end), row.name, row.score, row.strand, row.internal_base, row.subfamily
        )

        allowed_cfs = internal_to_family.get(subf, set())
        ltr5_name = "."
        ltr3_name = "."
        ltr5_dist = ""
        ltr3_dist = ""

        new_start, new_end = istart, iend

        if allowed_cfs and (chrom, strand) in ltr_groups:
            G = ltr_groups[(chrom, strand)]
            dfL = G["df"]
            starts = G["starts"]

            # ---- Upstream (LTR5 relative position)
            left_bound = istart - (ltr_gap + (iend - istart) + 1000)
            i0 = max(0, bisect_left(starts, left_bound))
            i1 = bisect_right(starts, istart)
            upstream = dfL.iloc[i0:i1]
            if len(upstream):
                up_hits = upstream[
                    (upstream["classfam"].isin(allowed_cfs))
                    & ((istart - upstream["end"]).between(0, ltr_gap))
                ]
                if len(up_hits):
                    cand_start, r5 = choose_upstream(up_hits, istart, strategy)
                    if cand_start is not None:
                        new_start = min(new_start, cand_start)
                        ltr5_name = r5["ltr_subfamily"]
                        ltr5_dist = str(istart - int(r5["end"]))

            # ---- Downstream (LTR3)
            j0 = bisect_left(starts, iend)
            j1 = bisect_right(starts, iend + ltr_gap)
            downstream = dfL.iloc[j0:j1]
            if len(downstream):
                dn_hits = downstream[
                    (downstream["classfam"].isin(allowed_cfs))
                    & ((downstream["start"] - iend).between(0, ltr_gap))
                ]
                if len(dn_hits):
                    cand_end, r3 = choose_downstream(dn_hits, iend, strategy)
                    if cand_end is not None:
                        new_end = max(new_end, cand_end)
                        ltr3_name = r3["ltr_subfamily"]
                        ltr3_dist = str(int(r3["start"]) - iend)

        proposed.append(
            {
                "chrom": chrom,
                "start": int(new_start),
                "end": int(new_end),
                "original_start": istart,
                "original_end": iend,
                "strand": strand,
                "internal_base": internal_base,  # e.g., 'HERVK13-int'
                "subfamily": subf,               # e.g., 'HERVK13'
                "score": score,
                "ltr5": ltr5_name,
                "ltr3": ltr3_name,
                "dist5": ltr5_dist,
                "dist3": ltr3_dist,
            }
        )

    proposed_df = pd.DataFrame(proposed)

    # -------- Neighbor-only overlap resolution (per chrom,strand)
    final_rows = []
    audit_rows = []
    gb = proposed_df.groupby(["chrom", "strand"], sort=False)

    # Per-subfamily stable counters
    counters = defaultdict(int)

    it2 = tqdm(gb, total=gb.ngroups, desc="Validating overlaps & naming")
    for (chrom, strand), grp in it2:
        grp = grp.sort_values(["start", "end"], kind="mergesort").reset_index(drop=True)
        for i, r in grp.iterrows():
            start = int(r["start"])
            end = int(r["end"])

            # Left neighbor
            if i > 0:
                left = grp.iloc[i - 1]
                if start < int(left["end"]):
                    start = int(r["original_start"])  # revert left extension only

            # Right neighbor
            if i + 1 < len(grp):
                right = grp.iloc[i + 1]
                if end > int(right["start"]):
                    end = int(r["original_end"])  # revert right extension only

            subfam_int = str(r["internal_base"])      # with '-int'
            subfam = str(r["subfamily"])              # without '-int'
            counters[subfam_int] += 1
            locid = f"{subfam_int}_{counters[subfam_int]:05d}"  # e.g., HERVK13-int_00027

            has5 = "1" if r["ltr5"] != "." else "0"
            has3 = "1" if r["ltr3"] != "." else "0"

            bed_name = (
                f"{subfam_int}|locid={locid}|has5={has5}|has3={has3}|"
                f"ltr5={r['ltr5']}|ltr3={r['ltr3']}"
            )

            # BED row
            final_rows.append([r["chrom"], start, end, bed_name, r["score"], r["strand"]])

            # Mapping/audit row
            audit_rows.append({
                "locid": locid,
                "chrom": r["chrom"],
                "start": start,
                "end": end,
                "strand": r["strand"],
                "internal_base": subfam_int,
                "subfamily": subfam,
                "original_start": int(r["original_start"]),
                "original_end": int(r["original_end"]),
                "has5": has5,
                "has3": has3,
                "ltr5": r["ltr5"],
                "ltr3": r["ltr3"],
                "dist5": r["dist5"],
                "dist3": r["dist3"],
            })

    final_df = pd.DataFrame(final_rows, columns=["chrom", "start", "end", "name", "score", "strand"])
    final_df = final_df.sort_values(["chrom", "start", "end"], kind="mergesort")
    final_df.to_csv(output_bed_path, sep="\t", header=False, index=False)

    if map_out:
        pd.DataFrame(audit_rows).sort_values(["chrom","start","end"]).to_csv(map_out, sep="\t", index=False)

    print(
        f"✅ Saved: {output_bed_path}\n"
        f"   Total loci: {len(final_df)} | gap={ltr_gap} | strategy={strategy} (nearest/furthest)\n"
        f"{'   Mapping TSV: ' + map_out if map_out else ''}"
    )

# ==============================
# CLI
# ==============================

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge internal ERV regions with nearby LTRs (strand & class/family aware) and produce stable locus IDs."
    )
    p.add_argument("--internal", required=True, help="Internal BED (merged fragment set)")
    p.add_argument("--ltr", required=True, help="Merged LTR BED")
    p.add_argument("--map", required=True, help="LTR dictionary TSV with columns: RepeatName, ClassFamily")
    p.add_argument("--out", required=True, help="Output BED")
    p.add_argument("--map-out", default=None, help="Optional: write mapping TSV (per-locus audit)")
    p.add_argument("--gap", type=int, default=150, help="Max gap (bp) to consider an LTR as adjacent (default: 150)")
    p.add_argument(
        "--strategy",
        choices=["nearest", "furthest"],
        default="nearest",
        help="Choose nearest or furthest valid LTR per side (default: nearest)",
    )
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(
        internal_bed_path=args.internal,
        ltr_bed_path=args.ltr,
        ltr_dict_path=args.map,
        output_bed_path=args.out,
        ltr_gap=args.gap,
        strategy=args.strategy,
        map_out=args.map_out,
    )
