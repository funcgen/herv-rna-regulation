#!/usr/bin/env python3
"""
Build the canonical RBP motif library in MEME format from the frozen ATtRACT source.

Selection logic (unchanged from the historical implementation):

  1. Restrict the ATtRACT database table to one organism   (Organism == Homo_sapiens)
  2. Restrict to one source database                       (Database == C, CisBP-RNA)
  3. Require a usable matrix identifier                    (Matrix_id not in {"", "N/A"})
  4. Keep the PWM of every retained Matrix_id from pwm.txt
  5. Emit MEME v4 with a DNA alphabet, writing the ATtRACT U column as T

Motifs are ordered deterministically by (representative gene name, matrix id), and each
matrix id is given a single representative gene: the lexicographically first of the gene
names ATtRACT associates with it. Both rules are inherited from the historical code and
are what make the output reproducible.

Scientific behaviour is identical to the historical scripts. Only paths, configuration,
error handling and reporting have changed.

Reads nothing but its declared inputs. Writes nothing but its declared outputs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# ATtRACT pwm.txt record header, e.g. ">M001_0.6 7"  (matrix id, then width)
PWM_HEADER_RE = re.compile(r"^>(\S+)\s+(\d+)\s*$")

# MEME probabilities are written with this precision. Changing it changes the
# output bytes, so it is pinned here rather than left to the formatter.
PROBABILITY_FORMAT = "{:.6f}"


def die(message: str) -> None:
    print(f"[ERROR] {message}", file=sys.stderr)
    sys.exit(1)


def read_attract_table(path: Path, args: argparse.Namespace):
    """Filter the ATtRACT table and return the retained matrix ids and gene mappings."""
    if not path.is_file():
        die(f"ATtRACT database table not found: {path}")

    null_matrix_ids = set(args.matrix_id_null_values)
    wanted: set[str] = set()
    genes_by_matrix: dict[str, set[str]] = defaultdict(set)
    n_rows = 0

    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            die(f"ATtRACT table has no header row: {path}")

        for column in (
            args.organism_column,
            args.database_column,
            args.matrix_id_column,
            args.gene_name_column,
        ):
            if column not in reader.fieldnames:
                die(
                    f"ATtRACT table is missing the required column '{column}'. "
                    f"Columns present: {', '.join(reader.fieldnames)}"
                )

        for row in reader:
            n_rows += 1
            if row.get(args.organism_column) != args.organism_value:
                continue
            if row.get(args.database_column) != args.database_value:
                continue

            matrix_id = (row.get(args.matrix_id_column) or "").strip()
            if matrix_id in null_matrix_ids:
                continue

            gene = str(row.get(args.gene_name_column) or "").strip()
            wanted.add(matrix_id)
            genes_by_matrix[matrix_id].add(gene)

    if n_rows == 0:
        die(f"ATtRACT table contains no data rows: {path}")
    if not wanted:
        die(
            "No matrix ids survived filtering. Check the organism and database "
            f"filters ({args.organism_value!r}, {args.database_value!r})."
        )

    # Representative gene per matrix: lexicographically first, for determinism.
    representative_gene = {
        matrix_id: sorted(genes)[0] for matrix_id, genes in genes_by_matrix.items()
    }
    return wanted, representative_gene, genes_by_matrix, n_rows


def parse_pwm_file(path: Path, wanted: set[str]):
    """Return {matrix_id: [[A, C, G, U], ...]} for the wanted matrix ids only."""
    if not path.is_file():
        die(f"ATtRACT PWM file not found: {path}")

    matrices: dict[str, list[list[float]]] = {}
    current_id: str | None = None
    current_rows: list[list[float]] = []
    declared_width: int | None = None
    n_headers = 0

    def store_current() -> None:
        if current_id is None or current_id not in wanted:
            return
        if declared_width is not None and len(current_rows) != declared_width:
            die(
                f"PWM {current_id}: header declares width {declared_width} "
                f"but {len(current_rows)} rows were read"
            )
        matrices[current_id] = current_rows

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            header = PWM_HEADER_RE.match(line)
            if header:
                store_current()
                current_id = header.group(1)
                declared_width = int(header.group(2))
                current_rows = []
                n_headers += 1
                continue

            if current_id is None:
                die(f"{path}:{line_number}: data row before any PWM header")

            fields = line.split()
            if len(fields) != 4:
                die(
                    f"{path}:{line_number}: expected 4 columns (A C G U), "
                    f"found {len(fields)}"
                )
            try:
                current_rows.append([float(value) for value in fields])
            except ValueError:
                die(f"{path}:{line_number}: non-numeric PWM value: {line!r}")

    store_current()

    if not matrices:
        die(
            "No PWMs matched the filtered matrix ids. Check that Matrix_id values "
            "correspond to pwm.txt headers such as '>M001_0.6 7'."
        )

    missing = sorted(wanted - set(matrices))
    return matrices, n_headers, missing


def write_meme(path: Path, matrices, representative_gene, genes_by_matrix, label: str) -> None:
    """Write the MEME v4 library. Ordering is (representative gene, matrix id)."""
    ordered_ids = sorted(
        matrices,
        key=lambda matrix_id: (representative_gene.get(matrix_id, "NA"), str(matrix_id)),
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        out.write("MEME version 4\n\n")
        out.write(f"# Built from ATtRACT; label={label}\n")
        out.write("ALPHABET= ACGT\n")
        out.write(
            "# NOTE: ATtRACT PWMs are A C G U; the 4th column is written as T "
            "for DNA scanning.\n"
        )
        out.write(
            "# Recommended usage: scan strand-oriented genomic FASTA "
            "(bedtools getfasta -s) and keep FIMO strands=+.\n"
        )
        out.write("strands: +\n\n")
        out.write("Background letter frequencies:\n")
        out.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

        for matrix_id in ordered_ids:
            rows = matrices[matrix_id]
            if not rows:
                continue

            gene = representative_gene.get(matrix_id, "NA")
            all_genes = sorted(genes_by_matrix.get(matrix_id, []))
            genes_field = ",".join(all_genes) if all_genes else "NA"

            out.write(f"MOTIF ATtRACT_{matrix_id} {gene}|{matrix_id}\n")
            out.write(f"# matrix_id={matrix_id}; genes={genes_field}\n")
            out.write(
                "letter-probability matrix: alength= 4 "
                f"w= {len(rows)} nsites= 20 E= 0\n"
            )
            for a, c, g, u in rows:
                out.write(
                    "\t".join(
                        PROBABILITY_FORMAT.format(value) for value in (a, c, g, u)
                    )
                    + "\n"
                )
            out.write("\n")


def write_manifest(path: Path, matrices, representative_gene, genes_by_matrix, args) -> None:
    """One row per matrix id, recording what entered the library and why."""
    ordered_ids = sorted(
        matrices,
        key=lambda matrix_id: (representative_gene.get(matrix_id, "NA"), str(matrix_id)),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "motif_id",
                "rbp",
                "source_database",
                "species",
                "matrix_length",
                "source_record_id",
                "included",
                "notes",
            ]
        )
        for matrix_id in ordered_ids:
            all_genes = sorted(genes_by_matrix.get(matrix_id, []))
            note = ""
            if len(all_genes) > 1:
                note = f"matrix shared by {len(all_genes)} genes: {','.join(all_genes)}"
            writer.writerow(
                [
                    f"ATtRACT_{matrix_id}",
                    representative_gene.get(matrix_id, "NA"),
                    args.database_value,
                    args.organism_value,
                    len(matrices[matrix_id]),
                    matrix_id,
                    "yes",
                    note,
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the canonical RBP MEME motif library from frozen ATtRACT inputs."
    )
    parser.add_argument("--attract-db", required=True, type=Path)
    parser.add_argument("--pwm", required=True, type=Path)
    parser.add_argument("--out-meme", required=True, type=Path)
    parser.add_argument("--out-manifest", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--organism-column", default="Organism")
    parser.add_argument("--organism-value", default="Homo_sapiens")
    parser.add_argument("--database-column", default="Database")
    parser.add_argument("--database-value", default="C")
    parser.add_argument("--matrix-id-column", default="Matrix_id")
    parser.add_argument("--gene-name-column", default="Gene_name")
    parser.add_argument("--matrix-id-null-values", nargs="*", default=["N/A", ""])
    parser.add_argument(
        "--expect-unique-rbps",
        type=int,
        required=True,
        help="Acceptance target (VAL04). The build fails if this is not reproduced.",
    )
    parser.add_argument("--expect-motifs", type=int, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    wanted, representative_gene, genes_by_matrix, n_rows = read_attract_table(
        args.attract_db, args
    )
    matrices, n_headers, missing = parse_pwm_file(args.pwm, wanted)

    n_motifs = len(matrices)
    unique_rbps = {representative_gene[matrix_id] for matrix_id in matrices}
    n_rbps = len(unique_rbps)

    print(f"[info] ATtRACT rows read           : {n_rows}")
    print(f"[info] matrix ids after filtering  : {len(wanted)}")
    print(f"[info] PWM records scanned         : {n_headers}")
    print(f"[info] motifs written              : {n_motifs}")
    print(f"[info] unique RBPs                 : {n_rbps}")
    if missing:
        print(
            f"[warn] {len(missing)} filtered matrix ids had no PWM record: "
            f"{', '.join(missing[:10])}"
            + (" ..." if len(missing) > 10 else ""),
            file=sys.stderr,
        )

    # Fail rather than silently emit a library with the wrong biological universe.
    if n_rbps != args.expect_unique_rbps:
        die(
            f"unique RBP count is {n_rbps}, expected {args.expect_unique_rbps}. "
            "The tested RBP universe does not match the canonical library; "
            "refusing to write a library downstream analyses would misinterpret."
        )
    if n_motifs != args.expect_motifs:
        die(f"motif count is {n_motifs}, expected {args.expect_motifs}.")

    write_meme(args.out_meme, matrices, representative_gene, genes_by_matrix, args.label)
    write_manifest(
        args.out_manifest, matrices, representative_gene, genes_by_matrix, args
    )

    print(f"[ok] wrote {args.out_meme} ({n_motifs} motifs, {n_rbps} RBPs)")
    print(f"[ok] wrote {args.out_manifest}")


if __name__ == "__main__":
    main()
