#!/usr/bin/env python3
import argparse
import re
import sys
from collections import defaultdict
import pandas as pd

HDR_RE = re.compile(r"^>(\S+)\s+(\d+)\s*$")

def die(msg: str, code: int = 1):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def read_attract(attract_tsv: str):
    df = pd.read_csv(attract_tsv, sep="\t", dtype=str)

    required = {"Gene_name", "Matrix_id"}
    missing = required - set(df.columns)
    if missing:
        die(f"ATtRACT TSV missing required columns: {', '.join(sorted(missing))}")

    df = df[df["Matrix_id"].notna()].copy()
    df["Matrix_id"] = df["Matrix_id"].astype(str)
    df = df[(df["Matrix_id"] != "") & (df["Matrix_id"] != "N/A")]

    wanted = set(df["Matrix_id"].tolist())

    # Deterministic representative gene name per matrix_id
    mid_to_gene = (
        df.groupby("Matrix_id")["Gene_name"]
          .apply(lambda x: sorted(set(map(str, x)))[0])
          .to_dict()
    )

    # Also keep all genes per matrix (sometimes shared)
    mid_to_allgenes = defaultdict(set)
    for _, r in df.iterrows():
        mid_to_allgenes[r["Matrix_id"]].add(str(r["Gene_name"]))

    return wanted, mid_to_gene, mid_to_allgenes

def parse_pwm(pwm_txt: str, wanted_ids: set, log_every: int = 500):
    pwm = {}  # mid -> list of rows [A,C,G,U]
    cur_id = None
    cur_rows = []
    seen_headers = 0
    kept = 0

    def save_current():
        nonlocal kept
        if cur_id is not None and cur_id in wanted_ids:
            pwm[cur_id] = cur_rows
            kept += 1

    with open(pwm_txt) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            m = HDR_RE.match(line)
            if m:
                # reached a new PWM header; save previous
                save_current()
                cur_id = m.group(1)
                cur_rows = []
                seen_headers += 1

                if seen_headers % log_every == 0:
                    print(f"[parse_pwm] seen={seen_headers} kept={kept} (wanted={len(wanted_ids)})",
                          file=sys.stderr)
                continue

            parts = line.split()
            if len(parts) != 4 or cur_id is None:
                continue
            try:
                row = [float(x) for x in parts]
            except ValueError:
                continue
            cur_rows.append(row)

    save_current()
    print(f"[parse_pwm] done: seen={seen_headers} kept={kept}", file=sys.stderr)
    return pwm


def write_meme(out_meme: str, pwm: dict, mid_to_gene: dict, mid_to_allgenes: dict, label: str):
    # Deterministic ordering: by (gene, matrix_id)
    ordered = sorted(pwm.keys(), key=lambda mid: (mid_to_gene.get(mid, "NA"), str(mid)))

    with open(out_meme, "w") as out:
        out.write("MEME version 4\n\n")
        out.write(f"# Built from ATtRACT; label={label}\n")
        out.write("ALPHABET= ACGT\n")
        out.write("# NOTE: ATtRACT PWMs are A C G U; the 4th column is written as T for DNA scanning.\n")
        out.write("# Recommended usage: scan strand-oriented genomic FASTA (bedtools getfasta -s) and keep FIMO strands=+.\n")
        out.write("strands: +\n\n")
        out.write("Background letter frequencies:\n")
        out.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

        for mid in ordered:
            rows = pwm[mid]
            if not rows:
                continue

            gene = mid_to_gene.get(mid, "NA")
            allgenes = sorted(mid_to_allgenes.get(mid, []))
            # keep a short deterministic alt name; full list in comment
            motif_id = f"ATtRACT_{mid}"
            motif_name = f"{gene}|{mid}"

            out.write(f"MOTIF {motif_id} {motif_name}\n")
            out.write(f"# matrix_id={mid}; genes={','.join(allgenes) if allgenes else 'NA'}\n")

            w = len(rows)
            out.write(f"letter-probability matrix: alength= 4 w= {w} nsites= 20 E= 0\n")

            # ATtRACT rows are A C G U; we treat U as T for DNA scanning
            for a, c, g, u in rows:
                out.write(f"{a:.6f}\t{c:.6f}\t{g:.6f}\t{u:.6f}\n")
            out.write("\n")

def main():
    ap = argparse.ArgumentParser(description="Build a MEME motif file from ATtRACT_db.tsv + pwm.txt")
    ap.add_argument("--attract-tsv", required=True, help="Filtered ATtRACT TSV containing Gene_name and Matrix_id")
    ap.add_argument("--pwm-txt", required=True, help="ATtRACT pwm.txt file")
    ap.add_argument("--out-meme", required=True, help="Output MEME file")
    ap.add_argument("--label", default="ATtRACT", help="Label to write in MEME header")
    args = ap.parse_args()

    wanted, mid_to_gene, mid_to_allgenes = read_attract(args.attract_tsv)
    pwm = parse_pwm(args.pwm_txt, wanted)

    if len(pwm) == 0:
        die("No PWMs were parsed. Check that Matrix_id values match pwm.txt headers (e.g., '>M001_0.6 7').")

    write_meme(args.out_meme, pwm, mid_to_gene, mid_to_allgenes, args.label)

    print(f"[OK] Wrote {args.out_meme} with {len(pwm)} motifs")

if __name__ == "__main__":
    main()
