#!/usr/bin/env python3
"""
06_lncrna_function step 2 (frozen global order 39, effective order 40), LOCAL, OFFLINE:
extract lncRNA function statements from the frozen PubMed corpus (SRC16, D031).

REFACTOR_SCRIPT of analysis/scripts/scripts_categories/F_literature_mining/02_extract_lncrna_functions.py
(sha256 b3316531..., the producer of the canonical lncrna_functions_extracted_v5.csv). The keyword
sets, lncRNA patterns, IGNORE_EXACT / WEAK_LIST_CUES and every helper function are copied VERBATIM;
extract_functions() is unchanged apart from explicit paths (no silent return on a missing input or
an empty result) and a quieter progress bar. Release additions, marked RELEASE:
  * sockets disabled for the whole run (the corpus is frozen; nothing may be fetched);
  * the input must be the AVAILABLE SRC16 snapshot (sha256-checked); the SUPERSEDED earlier corpus
    is refused;
  * the output is written to a staging directory and published atomically.
Algorithm (historical v5, unchanged): every corpus row is processed independently (a PMID present in
two year slices is analysed twice); only the Abstract is mined (not the Title); spaCy en_core_web_sm
splits sentences and gives lemmas; sentences shorter than 20 characters are skipped; lncRNA candidates
are regex matches (9 patterns) filtered by is_valid_candidate(), normalised (upper-case, whitespace
removed), de-duplicated and sorted per sentence; a sentence is kept when it has >= 1 candidate and
>= 1 keyword (lemma or lower-cased token in the five keyword sets), except list-like sentences without
mechanistic or phenotypic terms; one row per candidate; exact duplicate rows are dropped (first kept).
No negation, species or disease filter exists. The environment is pinned (VALIDATED_RELEASE_ENVIRONMENT).
"""
# RELEASE: network prohibition guard, before any other import can open a connection
import socket


def _no_network(*_args, **_kwargs):
    raise RuntimeError("network access is disabled in the release build (frozen SRC16 corpus)")


# patch the connection entry points, not the socket class (ssl subclasses socket.socket)
socket.socket.connect = _no_network
socket.socket.connect_ex = _no_network
socket.create_connection = _no_network
socket.getaddrinfo = _no_network

import argparse
import csv
import hashlib
import os
import shutil
import tempfile

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
# Main function (historical, paths explicit)
# ------------------------------------------------------------

def extract_functions(input_csv, output_csv):
    print(f"Reading dataset from {input_csv}...")

    df = pd.read_csv(input_csv)

    required_cols = {"PMID", "Year", "Abstract"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    extracted_rows = []

    print("\nAnalyzing abstracts...")
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Extracting Functions", unit="abstract", mininterval=30):
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
        raise RuntimeError("No relationships found.")

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


# ------------------------------------------------------------
# RELEASE: frozen-source guard, CLI, atomic publication
# ------------------------------------------------------------
def verify_src16(corpus, sums_path, snapshots_path):
    name = os.path.basename(corpus)
    with open(snapshots_path, newline="", encoding="utf-8") as fh:
        status = [r["status"] for r in csv.DictReader(fh, delimiter="\t") if r["source_id"] == "SRC16" and r["snapshot_name"] == name]
    if status != ["AVAILABLE"]:
        raise ValueError(f"{name} is not the AVAILABLE SRC16 snapshot (status {status}); the superseded corpus is refused")
    with open(sums_path, newline="", encoding="utf-8") as fh:
        want = [r["sha256"] for r in csv.DictReader(fh, delimiter="\t") if r["source_id"] == "SRC16" and r["filename"] == name]
    h = hashlib.sha256()
    with open(corpus, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    if want != [h.hexdigest()]:
        raise ValueError(f"{corpus} does not match its SRC16 registration")
    print(f"[ok] SRC16 {name} sha256 {h.hexdigest()} (registered, AVAILABLE)", flush=True)


def main():
    ap = argparse.ArgumentParser(description="Extract lncRNA function statements (historical v5), offline.")
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--snapshot-checksums", required=True)
    ap.add_argument("--snapshots", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    verify_src16(a.corpus, a.snapshot_checksums, a.snapshots)
    print(f"[env] spacy {spacy.__version__}, model en_core_web_sm {nlp.meta.get('version')}, pandas {pd.__version__}", flush=True)
    parent = os.path.dirname(os.path.abspath(a.outdir)) or "."
    os.makedirs(parent, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".stage_function_statements.", dir=parent)
    extract_functions(a.corpus, os.path.join(stage, "lncrna_functions_extracted_v5.csv"))
    os.chmod(os.path.join(stage, "lncrna_functions_extracted_v5.csv"), 0o644)
    if os.path.exists(a.outdir):
        shutil.rmtree(a.outdir)
    os.rename(stage, a.outdir)
    os.chmod(a.outdir, 0o755)
    print(f"[OK] wrote {a.outdir}")


if __name__ == "__main__":
    main()
