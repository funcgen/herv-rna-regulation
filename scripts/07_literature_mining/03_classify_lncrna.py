import pandas as pd
import re
from pathlib import Path

# =========================================================
# INPUT
# =========================================================

FILES = [
    "../lncRNAs/interaction-human-tsv.tsv",
    "../lncRNAs/lncRNA-human-tsv.tsv",
    "../lncRNAs/not-disease-human-tsv.tsv",
    "../lncRNAs/Disease-human-tsv.tsv"
]

GTF_FILE = "../data/gencode.v48.primary_assembly.annotation.gtf"

OUT_PREFIX = "../lncRNAs/atlas/lncrna_master_resource"


# =========================================================
# HELPERS
# =========================================================

def norm_col(c):
    c = str(c).strip().lower()
    c = re.sub(r'[\s\-/]+', '_', c)
    c = re.sub(r'[^a-z0-9_]', '', c)
    return c


def find_first_matching_column(cols, patterns):
    for p in patterns:
        for c in cols:
            if re.search(p, c):
                return c
    return None


def extract_coordinates(text):
    """
    Parse coordinates from strings like:
      chr1:12345-67890(+)
      chr1:12345-67890
      1:12345-67890
      chr1 12345 67890 +
    """
    if pd.isna(text):
        return None, None, None, None

    s = str(text).strip()

    m = re.search(r'(chr[\w]+)[:\s]+(\d+)[-\s]+(\d+)\s*\(?([+-])?\)?', s, re.I)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)

    m = re.search(r'(^|\s)([\w]+)[:\s]+(\d+)[-\s]+(\d+)\s*\(?([+-])?\)?', s, re.I)
    if m:
        chrom = m.group(2)
        if not str(chrom).lower().startswith("chr"):
            chrom = f"chr{chrom}"
        return chrom, int(m.group(3)), int(m.group(4)), m.group(5)

    m = re.search(r'(chr[\w]+)\s+(\d+)\s+(\d+)\s+([+-])', s, re.I)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)

    return None, None, None, None


def looks_like_ensembl_id(x):
    if pd.isna(x):
        return False
    s = str(x).strip()
    s = re.sub(r'\.\d+$', '', s)
    return s.startswith("ENSG") or s.startswith("ENST") or s.startswith("ENSP")


def parse_gtf_attributes(attr_string):
    attrs = {}
    if pd.isna(attr_string):
        return attrs

    for m in re.finditer(r'(\S+)\s+"([^"]+)"', str(attr_string)):
        key, value = m.group(1), m.group(2)
        attrs[key] = value
    return attrs


def bool_to_str(x):
    if pd.isna(x):
        return ""
    return "TRUE" if bool(x) else "FALSE"


def standardize_flag_columns(df, flag_cols):
    for col in flag_cols:
        if col in df.columns:
            df[col] = df[col].apply(bool_to_str)
    return df


def clean_ensembl(x):
    """
    Normalize Ensembl-like IDs by removing version suffixes.
    Works for ENSG, ENST, ENSP, etc.
    Example:
      ENSG00000230544.1 -> ENSG00000230544
      ENST00000422494.1 -> ENST00000422494
    """
    if pd.isna(x):
        return None
    s = str(x).strip()
    s = re.sub(r'\.\d+$', '', s)
    return s if s else None


def clean_text(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    s = re.sub(r'\s+', ' ', s)
    return s if s else None


def normalize_strand(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if s in ["+", "-"]:
        return s
    if s.lower() in ["forward", "plus", "positive"]:
        return "+"
    if s.lower() in ["reverse", "minus", "negative"]:
        return "-"
    return s


def is_missing(x):
    return pd.isna(x) or str(x).strip() == ""


def infer_source_kind(file):
    file = str(file).lower()
    if "interaction" in file:
        return "interaction"
    if "not-disease" in file:
        return "not_disease"
    if "disease" in file:
        return "disease"
    if "lncrna-human" in file:
        return "annotation"
    return "other"


def build_function_description(row, source_kind):
    desc = clean_text(row.get("description"))
    dysfunction = clean_text(row.get("dysfunction_type"))
    disease = clean_text(row.get("disease"))
    sample = clean_text(row.get("sample"))
    expr = clean_text(row.get("expression_pattern"))
    interaction_target = clean_text(row.get("interaction_target"))
    interaction_type = clean_text(row.get("type_of_interaction"))
    interaction_level = clean_text(row.get("level_of_interaction"))

    generic_terms = {"expression", "interaction", "regulation", "mutation", "locus", "epigenetics"}

    if source_kind == "interaction":
        parts = []
        if desc:
            parts.append(desc)
        else:
            if interaction_type and interaction_target:
                parts.append(f"{interaction_type} with {interaction_target}")
            elif interaction_type:
                parts.append(f"{interaction_type}-based interaction")

        if interaction_target and desc and interaction_target.lower() not in desc.lower():
            parts.append(f"target: {interaction_target}")

        if interaction_level:
            parts.append(f"interaction level: {interaction_level}")

        return " | ".join(parts) if parts else None

    if source_kind == "not_disease":
        parts = []
        if desc and desc.lower() not in generic_terms:
            parts.append(desc)
        elif dysfunction and dysfunction.lower() not in generic_terms:
            parts.append(dysfunction)

        context_bits = []
        if expr:
            context_bits.append(f"expression: {expr}")
        if sample:
            context_bits.append(f"sample: {sample}")
        if disease:
            context_bits.append(f"context: {disease}")

        if context_bits:
            parts.append(" / ".join(context_bits))

        return " | ".join(parts) if parts else None

    if source_kind == "disease":
        parts = []
        if desc and desc.lower() not in generic_terms:
            parts.append(desc)

        context_bits = []
        if disease:
            context_bits.append(f"disease: {disease}")
        if expr:
            context_bits.append(f"expression: {expr}")
        if sample:
            context_bits.append(f"sample: {sample}")
        if dysfunction:
            context_bits.append(f"evidence: {dysfunction}")

        if context_bits:
            parts.append(" / ".join(context_bits))

        return " | ".join(parts) if parts else None

    if source_kind == "annotation":
        return desc

    return desc


def extract_evidence_type(row, source_kind):
    dysfunction = clean_text(row.get("dysfunction_type"))
    interaction_type = clean_text(row.get("type_of_interaction"))

    if source_kind == "interaction":
        return interaction_type.lower() if interaction_type else "interaction"
    if source_kind == "disease":
        return dysfunction.lower() if dysfunction else "disease_association"
    if source_kind == "not_disease":
        return dysfunction.lower() if dysfunction else "functional_annotation"
    if source_kind == "annotation":
        return "annotation"

    return "unknown"


def collapse_unique(values):
    vals = []
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s and s not in vals:
            vals.append(s)
    return "; ".join(vals)


def read_tsv_flexibly(path):
    try:
        return pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.read_csv(path, sep="\t", dtype=str, engine="python", on_bad_lines="skip")


def load_gencode_lookups(gtf_file):
    """
    Parse GENCODE GTF and build lookup tables for:
      - lncRNA genes
      - all genes
      - all transcripts

    This allows:
      1. preferred matching to lncRNA genes
      2. fallback matching to any gene/transcript
    """
    lnc_types = {
        "lncrna",
        "lincrna",
        "antisense",
        "sense_intronic",
        "sense_overlapping",
        "processed_transcript",
        "3prime_overlapping_ncrna",
        "non_coding",
        "macro_lncrna",
        "bidirectional_promoter_lncrna"
    }

    gene_rows = []
    tx_rows = []

    with open(gtf_file, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            attr = parse_gtf_attributes(attrs)

            gene_id = clean_ensembl(attr.get("gene_id"))
            gene_name = clean_text(attr.get("gene_name"))
            gene_type = clean_text(attr.get("gene_type"))

            transcript_id = clean_ensembl(attr.get("transcript_id"))
            transcript_name = clean_text(attr.get("transcript_name"))
            transcript_type = clean_text(attr.get("transcript_type"))

            start = pd.to_numeric(start, errors="coerce")
            end = pd.to_numeric(end, errors="coerce")
            strand = normalize_strand(strand)

            if feature == "gene":
                gene_rows.append({
                    "ensembl_id": gene_id,
                    "lncrna_name": gene_name,
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_type_gencode": gene_type,
                    "is_lnc_like": (gene_type is not None and gene_type.lower() in lnc_types)
                })

            elif feature == "transcript":
                tx_rows.append({
                    "transcript_id": transcript_id,
                    "transcript_name": transcript_name,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_type_gencode": gene_type,
                    "transcript_type_gencode": transcript_type,
                    "is_lnc_like": (
                        (gene_type is not None and gene_type.lower() in lnc_types) or
                        (transcript_type is not None and transcript_type.lower() in lnc_types)
                    )
                })

    genes_df = pd.DataFrame(gene_rows).drop_duplicates()
    tx_df = pd.DataFrame(tx_rows).drop_duplicates()

    lnc_gene_df = genes_df[genes_df["is_lnc_like"] == True].copy()

    def make_lookup(df, key_col):
        if df.empty or key_col not in df.columns:
            return None
        tmp = df[df[key_col].notna() & (df[key_col].astype(str).str.strip() != "")]
        if tmp.empty:
            return None
        return tmp.drop_duplicates(key_col).set_index(key_col)

    lookups = {
        "lnc_gene_by_ensg": make_lookup(lnc_gene_df, "ensembl_id"),
        "lnc_gene_by_name": make_lookup(lnc_gene_df, "lncrna_name"),

        "all_gene_by_ensg": make_lookup(genes_df, "ensembl_id"),
        "all_gene_by_name": make_lookup(genes_df, "lncrna_name"),

        "tx_by_enst": make_lookup(tx_df, "transcript_id"),
        "tx_by_txname": make_lookup(tx_df, "transcript_name"),
    }

    return genes_df, tx_df, lookups


# =========================================================
# LOAD TABLES
# =========================================================

all_rows = []
annotation_lookup_rows = []

for file in FILES:
    path = Path(file)
    if not path.exists():
        print(f"WARNING: file not found: {file}")
        continue

    df = read_tsv_flexibly(path)

    if df.empty:
        print(f"WARNING: empty file: {file}")
        continue

    df.columns = [norm_col(c) for c in df.columns]
    cols = list(df.columns)

    print(f"\nProcessing: {file}")
    print(f"Columns: {cols}")

    name_col = find_first_matching_column(cols, [
        r'lncrna', r'^name$', r'gene_name', r'gene$', r'symbol', r'lnc_name'
    ])
    ensembl_col = find_first_matching_column(cols, [
        r'ensembl', r'gene_id', r'ensg'
    ])
    coord_col = find_first_matching_column(cols, [
        r'coordinate', r'genomic_location', r'location', r'locus', r'position'
    ])
    chrom_col = find_first_matching_column(cols, [
        r'^chr$', r'^chrom$', r'chromosome'
    ])
    start_col = find_first_matching_column(cols, [
        r'^start$', r'^start_site$', r'genomic_start'
    ])
    end_col = find_first_matching_column(cols, [
        r'^end$', r'^end_site$', r'genomic_end'
    ])
    strand_col = find_first_matching_column(cols, [
        r'^strand$', r'^chain$'
    ])

    source_kind = infer_source_kind(file)

    for _, row in df.iterrows():
        row_dict = {c: row[c] for c in df.columns}

        lncrna_name = clean_text(row_dict.get(name_col)) if name_col else None
        ensembl_id = clean_ensembl(row_dict.get(ensembl_col)) if ensembl_col else None

        if is_missing(ensembl_id) and looks_like_ensembl_id(lncrna_name):
            ensembl_id = clean_ensembl(lncrna_name)

        chromosome = clean_text(row_dict.get(chrom_col)) if chrom_col else None
        start = row_dict.get(start_col) if start_col else None
        end = row_dict.get(end_col) if end_col else None
        strand = normalize_strand(row_dict.get(strand_col)) if strand_col else None
        coordinates = clean_text(row_dict.get(coord_col)) if coord_col else None

        if is_missing(chromosome) or is_missing(start) or is_missing(end):
            c_chr, c_start, c_end, c_strand = extract_coordinates(coordinates)
            chromosome = chromosome if not is_missing(chromosome) else c_chr
            start = start if not is_missing(start) else c_start
            end = end if not is_missing(end) else c_end
            strand = strand if not is_missing(strand) else c_strand

        if not is_missing(chromosome):
            chromosome = str(chromosome).strip()
            if not chromosome.lower().startswith("chr"):
                chromosome = f"chr{chromosome}"

        start = pd.to_numeric(start, errors="coerce")
        end = pd.to_numeric(end, errors="coerce")

        if source_kind == "annotation":
            annotation_lookup_rows.append({
                "lncrna_name": lncrna_name,
                "ensembl_id": ensembl_id,
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "strand": strand
            })

        function_description_raw = clean_text(row_dict.get("description"))
        function_description_clean = build_function_description(row_dict, source_kind)
        evidence_type = extract_evidence_type(row_dict, source_kind)

        interaction_target = clean_text(row_dict.get("interaction_target"))
        interaction_type = clean_text(row_dict.get("type_of_interaction"))
        level_of_interaction = clean_text(row_dict.get("level_of_interaction"))
        disease_context = clean_text(row_dict.get("disease"))
        sample_context = clean_text(row_dict.get("sample"))
        expression_pattern = clean_text(row_dict.get("expression_pattern"))

        if all(is_missing(x) for x in [lncrna_name, ensembl_id, coordinates]):
            continue

        all_rows.append({
            "lncrna_name": lncrna_name,
            "ensembl_id": ensembl_id,
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "strand": strand,
            "coordinates": coordinates,

            "record_type": source_kind,
            "evidence_type": evidence_type,

            "function_description_raw": function_description_raw,
            "function_description_clean": function_description_clean,

            "disease_context": disease_context,
            "sample_context": sample_context,
            "expression_pattern": expression_pattern,

            "interaction_target": interaction_target,
            "interaction_type": interaction_type,
            "level_of_interaction": level_of_interaction,

            "source_file": file
        })


# =========================================================
# BUILD MASTER
# =========================================================

master = pd.DataFrame(all_rows)

if master.empty:
    raise ValueError("No valid rows were extracted from the TSV files.")

for c in [
    "lncrna_name", "ensembl_id", "chromosome", "strand", "coordinates",
    "record_type", "evidence_type",
    "function_description_raw", "function_description_clean",
    "disease_context", "sample_context", "expression_pattern",
    "interaction_target", "interaction_type", "level_of_interaction",
    "source_file", "coordinate_ref_source",
    "gencode_match_type", "gencode_biotype"
]:
    if c in master.columns:
        master[c] = master[c].astype("string").str.strip()


# =========================================================
# FILL / REFINE COORDINATES
# Priority:
#   1) GENCODE
#   2) annotation table (lncRNA-human-tsv.tsv)
# =========================================================

annotation_lookup = pd.DataFrame(annotation_lookup_rows)

by_annot_ensembl = None
by_annot_name = None

if not annotation_lookup.empty:
    annotation_lookup = annotation_lookup.drop_duplicates()

    by_annot_ensembl = (
        annotation_lookup[
            annotation_lookup["ensembl_id"].notna() &
            (annotation_lookup["ensembl_id"].astype(str).str.strip() != "")
        ]
        .drop_duplicates("ensembl_id")
        .set_index("ensembl_id")
    )

    by_annot_name = (
        annotation_lookup[
            annotation_lookup["lncrna_name"].notna() &
            (annotation_lookup["lncrna_name"].astype(str).str.strip() != "")
        ]
        .drop_duplicates("lncrna_name")
        .set_index("lncrna_name")
    )

genes_df, tx_df, gencode_lookups = load_gencode_lookups(GTF_FILE)

print("\nLoaded GENCODE genes:", 0 if genes_df is None else len(genes_df))
print("Loaded GENCODE transcripts:", 0 if tx_df is None else len(tx_df))


def fill_coords_with_priority(row):
    """
    Coordinate priority:
      1. GENCODE lncRNA gene by ENSG
      2. GENCODE lncRNA gene by gene name
      3. GENCODE transcript by ENST
      4. GENCODE any gene by ENSG
      5. GENCODE any gene by gene name
      6. annotation table by Ensembl
      7. annotation table by gene name
      8. keep original
    """
    raw_id = clean_ensembl(row["ensembl_id"])
    name = clean_text(row["lncrna_name"])

    if is_missing(raw_id) and looks_like_ensembl_id(name):
        raw_id = clean_ensembl(name)

    hit = None
    coord_source = None
    gencode_biotype = None
    gencode_match_type = None

    is_ensg = pd.notna(raw_id) and str(raw_id).startswith("ENSG")
    is_enst = pd.notna(raw_id) and str(raw_id).startswith("ENST")

    if is_ensg and gencode_lookups["lnc_gene_by_ensg"] is not None and raw_id in gencode_lookups["lnc_gene_by_ensg"].index:
        hit = gencode_lookups["lnc_gene_by_ensg"].loc[raw_id]
        coord_source = "GENCODE_lnc_gene_ENSG"
        gencode_biotype = hit.get("gene_type_gencode")
        gencode_match_type = "gene"

    elif pd.notna(name) and gencode_lookups["lnc_gene_by_name"] is not None and name in gencode_lookups["lnc_gene_by_name"].index:
        hit = gencode_lookups["lnc_gene_by_name"].loc[name]
        coord_source = "GENCODE_lnc_gene_name"
        gencode_biotype = hit.get("gene_type_gencode")
        gencode_match_type = "gene"

    elif is_enst and gencode_lookups["tx_by_enst"] is not None and raw_id in gencode_lookups["tx_by_enst"].index:
        hit = gencode_lookups["tx_by_enst"].loc[raw_id]
        coord_source = "GENCODE_transcript_ENST"
        gencode_biotype = hit.get("transcript_type_gencode") or hit.get("gene_type_gencode")
        gencode_match_type = "transcript"

        if is_missing(row["lncrna_name"]) and pd.notna(hit.get("gene_name")):
            row["lncrna_name"] = hit.get("gene_name")
        if is_missing(row["ensembl_id"]) and pd.notna(hit.get("gene_id")):
            row["ensembl_id"] = hit.get("gene_id")

    elif is_ensg and gencode_lookups["all_gene_by_ensg"] is not None and raw_id in gencode_lookups["all_gene_by_ensg"].index:
        hit = gencode_lookups["all_gene_by_ensg"].loc[raw_id]
        coord_source = "GENCODE_any_gene_ENSG"
        gencode_biotype = hit.get("gene_type_gencode")
        gencode_match_type = "gene"

    elif pd.notna(name) and gencode_lookups["all_gene_by_name"] is not None and name in gencode_lookups["all_gene_by_name"].index:
        hit = gencode_lookups["all_gene_by_name"].loc[name]
        coord_source = "GENCODE_any_gene_name"
        gencode_biotype = hit.get("gene_type_gencode")
        gencode_match_type = "gene"

    elif by_annot_ensembl is not None and pd.notna(raw_id) and raw_id in by_annot_ensembl.index:
        hit = by_annot_ensembl.loc[raw_id]
        coord_source = "annotation_ensembl"
        gencode_match_type = "annotation"

    elif by_annot_name is not None and pd.notna(name) and name in by_annot_name.index:
        hit = by_annot_name.loc[name]
        coord_source = "annotation_name"
        gencode_match_type = "annotation"

    if hit is not None:
        row["chromosome"] = hit.get("chromosome")
        row["start"] = hit.get("start")
        row["end"] = hit.get("end")
        row["strand"] = hit.get("strand")
        row["coordinate_ref_source"] = coord_source

        if is_missing(row.get("ensembl_id")):
            if pd.notna(hit.get("ensembl_id")):
                row["ensembl_id"] = hit.get("ensembl_id")
            elif pd.notna(hit.get("gene_id")):
                row["ensembl_id"] = hit.get("gene_id")

        if is_missing(row.get("lncrna_name")):
            if pd.notna(hit.get("lncrna_name")):
                row["lncrna_name"] = hit.get("lncrna_name")
            elif pd.notna(hit.get("gene_name")):
                row["lncrna_name"] = hit.get("gene_name")
    else:
        row["coordinate_ref_source"] = "original"

    row["gencode_match_type"] = gencode_match_type
    row["gencode_biotype"] = gencode_biotype

    return row


master = master.apply(fill_coords_with_priority, axis=1)

master["ensembl_id"] = master["ensembl_id"].apply(clean_ensembl)
master["lncrna_name"] = master["lncrna_name"].apply(clean_text)

master["merge_key"] = master["ensembl_id"].fillna("")
master.loc[master["merge_key"] == "", "merge_key"] = master["lncrna_name"].fillna("")
master.loc[master["merge_key"] == "", "merge_key"] = master["coordinates"].fillna("UNKNOWN")


# =========================================================
# COLLAPSE DUPLICATES
# =========================================================

collapsed = (
    master.groupby("merge_key", dropna=False)
    .agg({
        "lncrna_name": collapse_unique,
        "ensembl_id": collapse_unique,
        "chromosome": collapse_unique,
        "start": "first",
        "end": "first",
        "strand": collapse_unique,
        "coordinates": collapse_unique,

        "record_type": collapse_unique,
        "evidence_type": collapse_unique,

        "function_description_raw": collapse_unique,
        "function_description_clean": collapse_unique,

        "disease_context": collapse_unique,
        "sample_context": collapse_unique,
        "expression_pattern": collapse_unique,

        "interaction_target": collapse_unique,
        "interaction_type": collapse_unique,
        "level_of_interaction": collapse_unique,

        "source_file": collapse_unique,
        "coordinate_ref_source": collapse_unique,

        "gencode_match_type": collapse_unique,
        "gencode_biotype": collapse_unique,
    })
    .reset_index()
)

collapsed["has_ensembl"] = collapsed["ensembl_id"].notna() & (collapsed["ensembl_id"] != "")

collapsed["has_coordinates"] = (
    collapsed["chromosome"].notna() &
    (collapsed["chromosome"] != "") &
    collapsed["start"].notna() &
    collapsed["end"].notna()
)

collapsed["has_function"] = (
    collapsed["function_description_clean"].notna() &
    (collapsed["function_description_clean"] != "")
)

n_with_coordinates = collapsed["has_coordinates"].sum()

collapsed["coordinates_final"] = collapsed.apply(
    lambda r: (
        f"{r['chromosome']}:{int(r['start'])}-{int(r['end'])}({r['strand']})"
        if pd.notna(r["chromosome"]) and str(r["chromosome"]).strip() != ""
        and pd.notna(r["start"]) and pd.notna(r["end"])
        and pd.notna(r["strand"]) and str(r["strand"]).strip() != ""
        else (
            f"{r['chromosome']}:{int(r['start'])}-{int(r['end'])}"
            if pd.notna(r["chromosome"]) and str(r["chromosome"]).strip() != ""
            and pd.notna(r["start"]) and pd.notna(r["end"])
            else r["coordinates"]
        )
    ),
    axis=1
)

collapsed = collapsed.sort_values(
    by=["has_ensembl", "has_coordinates", "has_function", "lncrna_name"],
    ascending=[False, False, False, True]
)

flag_cols = ["has_ensembl", "has_coordinates", "has_function"]
collapsed = standardize_flag_columns(collapsed, flag_cols)


# =========================================================
# EXPORT
# =========================================================

master.to_csv(f"{OUT_PREFIX}.raw_harmonized.tsv", sep="\t", index=False)
collapsed.to_csv(f"{OUT_PREFIX}.tsv", sep="\t", index=False)

with pd.ExcelWriter(f"{OUT_PREFIX}.xlsx", engine="openpyxl") as writer:
    collapsed.to_excel(writer, sheet_name="master_resource", index=False)
    master.to_excel(writer, sheet_name="raw_harmonized_rows", index=False)

print("\nDone.")
print(f"Raw rows: {len(master)}")
print(f"Collapsed rows: {len(collapsed)}")
print("Files written:")
print(f" - {OUT_PREFIX}.raw_harmonized.tsv")
print(f" - {OUT_PREFIX}.tsv")
print(f" - {OUT_PREFIX}.xlsx")

print("\nRows with coordinates:", n_with_coordinates, "/", len(collapsed))