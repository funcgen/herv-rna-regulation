import pandas as pd
import spacy
import re
from tqdm import tqdm

print("Loading NLP model...")
nlp = spacy.load("en_core_web_sm")

# ------------------------------------------------------------
# Keyword sets
# ------------------------------------------------------------

MECHANISTIC_TERMS = {
    "sponge", "sponging", "cerna", "ceRNA", "bind", "binding", "recruit",
    "recruitment", "interact", "interaction", "stabilize", "stability",
    "scaffold", "transcription", "transcriptional", "methylation",
    "epigenetic", "chromatin", "target", "targeting", "ubiquitination"
}

PHENOTYPIC_TERMS = {
    "proliferation", "apoptosis", "migration", "invasion", "metastasis",
    "differentiation", "organogenesis", "autophagy", "ferroptosis",
    "glycolysis", "senescence", "resistance", "chemoresistance",
    "survival"
}

ASSOCIATIVE_TERMS = {
    "expression", "biomarker", "prognosis", "diagnostic", "diagnosis",
    "recurrence", "signature", "associated", "correlated"
}

PATHWAY_TERMS = {
    "pathway", "signaling"
}

REGULATORY_VERBS = {
    "promote", "suppress", "inhibit", "regulate", "mediate",
    "upregulate", "downregulate", "activate", "repress", "attenuate"
}

ALL_TERMS = {
    *(t.lower() for t in MECHANISTIC_TERMS),
    *(t.lower() for t in PHENOTYPIC_TERMS),
    *(t.lower() for t in ASSOCIATIVE_TERMS),
    *(t.lower() for t in PATHWAY_TERMS),
    *(t.lower() for t in REGULATORY_VERBS),
}

# ------------------------------------------------------------
# LncRNA patterns
# ------------------------------------------------------------

# More conservative regex patterns
LNCRNA_PATTERNS = [
    r"\bLINC\d{3,5}\b",
    r"\b[A-Z0-9]+-AS\d+\b",
    r"\b[A-Z0-9]+-DT\b",
    r"\b[A-Z0-9]+-IT\d*\b",
    r"\b[A-Z0-9]+P\d+\b",
    r"\bENSG\d{11,}(?:\.\d+)?\b",
    r"\b[A-Z][A-Z0-9]{2,}(?:-[A-Z0-9]+)+\b",
    r"\b[A-Z][A-Z0-9]{2,}\d+\b",
    r"\b[A-Z][A-Z0-9]{4,}\b",
]

LNCRNA_REGEX = re.compile("|".join(LNCRNA_PATTERNS))

IGNORE_EXACT = {
    "RNA", "DNA", "PCR", "WT", "MUT", "CI", "OS", "LNC", "LNCRNA", "LNCRNAS",
    "HCC", "NSCLC", "AML", "CRC", "LUAD", "GC", "EMT", "TCGA", "GEO", "RFS",
    "OSCC", "TMA", "RIP", "CHIP", "ISH", "RT-QPCR", "QRT-PCR", "MRNA", "MIRNA"
}

WEAK_LIST_CUES = {
    "identified", "selected", "included", "consists", "comprising", "signature",
    "model", "panel", "classifier", "risk", "nomogram", "screened"
}

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def normalize_candidate(name: str) -> str:
    return re.sub(r"\s+", "", name.strip().upper())

def is_valid_candidate(name: str) -> bool:
    x = normalize_candidate(name)

    if x in IGNORE_EXACT:
        return False

    # Very short uppercase tokens are often junk
    if len(x) <= 3 and not x.startswith("LINC"):
        return False

    # Reject pure Roman numeral-like things
    if re.fullmatch(r"[IVXLCDM]+", x):
        return False

    # Reject generic words ending in common biology suffixes only if too short
    if re.fullmatch(r"[A-Z]+", x) and len(x) <= 4:
        return False

    return True

def extract_candidates(sentence_text: str):
    raw = LNCRNA_REGEX.findall(sentence_text)
    valid = sorted({normalize_candidate(x) for x in raw if is_valid_candidate(x)})
    return valid

def collect_sentence_terms(sent):
    mech, pheno, assoc, pathway, verbs = set(), set(), set(), set(), set()

    for token in sent:
        lemma = token.lemma_.lower()
        text = token.text.lower()

        if lemma in {x.lower() for x in MECHANISTIC_TERMS} or text in {x.lower() for x in MECHANISTIC_TERMS}:
            mech.add(lemma)
        if lemma in {x.lower() for x in PHENOTYPIC_TERMS} or text in {x.lower() for x in PHENOTYPIC_TERMS}:
            pheno.add(lemma)
        if lemma in {x.lower() for x in ASSOCIATIVE_TERMS} or text in {x.lower() for x in ASSOCIATIVE_TERMS}:
            assoc.add(lemma)
        if lemma in {x.lower() for x in PATHWAY_TERMS} or text in {x.lower() for x in PATHWAY_TERMS}:
            pathway.add(lemma)
        if lemma in {x.lower() for x in REGULATORY_VERBS} or text in {x.lower() for x in REGULATORY_VERBS}:
            verbs.add(lemma)

    return {
        "mechanistic_terms": sorted(mech),
        "phenotypic_terms": sorted(pheno),
        "associative_terms": sorted(assoc),
        "pathway_terms": sorted(pathway),
        "regulatory_verbs": sorted(verbs),
    }

def detect_sentence_type(sentence_text: str) -> str:
    s = sentence_text.lower()

    if any(cue in s for cue in WEAK_LIST_CUES):
        return "list_or_signature_like"

    if "knockdown" in s or "silencing" in s or "overexpression" in s or "depletion" in s:
        return "functional_perturbation"

    if "bind" in s or "recruit" in s or "sponge" in s or "stability" in s or "transcription" in s:
        return "mechanistic_statement"

    if "associated with" in s or "correlated with" in s or "prognosis" in s:
        return "associative_statement"

    return "other"

def assign_evidence_tier(term_dict, sentence_type: str) -> str:
    has_mech = len(term_dict["mechanistic_terms"]) > 0
    has_pheno = len(term_dict["phenotypic_terms"]) > 0
    has_assoc = len(term_dict["associative_terms"]) > 0
    has_verbs = len(term_dict["regulatory_verbs"]) > 0

    if sentence_type == "mechanistic_statement" and (has_mech or (has_pheno and has_verbs)):
        return "Tier1_mechanistic"

    if sentence_type == "functional_perturbation" and (has_pheno or has_verbs):
        return "Tier2_phenotypic"

    if has_assoc and not (has_mech or has_pheno):
        return "Tier3_associative"

    if sentence_type == "list_or_signature_like":
        return "Tier4_weak"

    if has_mech:
        return "Tier1_mechanistic"
    if has_pheno:
        return "Tier2_phenotypic"
    if has_assoc:
        return "Tier3_associative"

    return "Tier4_weak"

def sentence_is_informative(candidates, term_dict, sentence_type: str) -> bool:
    if not candidates:
        return False

    total_terms = (
        len(term_dict["mechanistic_terms"]) +
        len(term_dict["phenotypic_terms"]) +
        len(term_dict["associative_terms"]) +
        len(term_dict["pathway_terms"]) +
        len(term_dict["regulatory_verbs"])
    )

    if total_terms == 0:
        return False

    # Reject weak list-like sentences unless they carry stronger evidence
    if sentence_type == "list_or_signature_like":
        if len(term_dict["mechanistic_terms"]) == 0 and len(term_dict["phenotypic_terms"]) == 0:
            return False

    return True

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

def extract_functions(
    input_csv="../results/lncrna_abstracts_dataset_full_no_function_filter.csv",
    output_csv="../results/lncrna_functions_extracted_v5.csv"
):
    print(f"Reading dataset from {input_csv}...")

    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print("Error: Could not find the input CSV.")
        return

    required_cols = {"PMID", "Year", "Abstract"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    extracted_rows = []

    print("\nAnalyzing abstracts...")
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Extracting Functions", unit="abstract"):
        pmid = row["PMID"]
        year = row["Year"]
        abstract = str(row["Abstract"])

        if abstract == "nan" or not abstract.strip():
            continue

        doc = nlp(abstract)

        for sent in doc.sents:
            sentence_text = sent.text.strip()

            if len(sentence_text) < 20:
                continue

            candidates = extract_candidates(sentence_text)
            if not candidates:
                continue

            term_dict = collect_sentence_terms(sent)
            sentence_type = detect_sentence_type(sentence_text)

            if not sentence_is_informative(candidates, term_dict, sentence_type):
                continue

            evidence_tier = assign_evidence_tier(term_dict, sentence_type)

            # One row per candidate lncRNA
            for candidate in candidates:
                extracted_rows.append({
                    "PMID": pmid,
                    "Year": year,
                    "Candidate_LncRNA": candidate,
                    "Mechanistic_Terms": "; ".join(term_dict["mechanistic_terms"]),
                    "Phenotypic_Terms": "; ".join(term_dict["phenotypic_terms"]),
                    "Associative_Terms": "; ".join(term_dict["associative_terms"]),
                    "Pathway_Terms": "; ".join(term_dict["pathway_terms"]),
                    "Regulatory_Verbs": "; ".join(term_dict["regulatory_verbs"]),
                    "Sentence_Type": sentence_type,
                    "Evidence_Tier": evidence_tier,
                    "Evidence_Sentence": sentence_text,
                    "n_candidates_in_sentence": len(candidates)
                })

    result_df = pd.DataFrame(extracted_rows)

    if result_df.empty:
        print("No relationships found.")
        return

    result_df = result_df.drop_duplicates()

    # Optional summary column similar to your previous output
    result_df["Functions_Detected"] = result_df.apply(
        lambda r: "; ".join(
            [x for x in [
                r["Mechanistic_Terms"],
                r["Phenotypic_Terms"],
                r["Associative_Terms"],
                r["Pathway_Terms"],
                r["Regulatory_Verbs"]
            ] if x]
        ),
        axis=1
    )

    result_df.to_csv(output_csv, index=False)

    print(f"\nExtraction complete! Found {len(result_df)} candidate relationships.")
    print(f"Results saved to {output_csv}")

# --- Run the script ---
extract_functions()