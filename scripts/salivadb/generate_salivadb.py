"""
generate_salivadb.py
Generate the SalivaDB biomarkers dataset from downloaded TSV files.
"""

from __future__ import annotations

import argparse
import gzip
import io
import logging
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------

BASE_DIR = Path("/data/shared/repos/extraction")
INPUT_DIR = BASE_DIR / "downloads/salivadb"
MAPPING_DIR = BASE_DIR / "downloads/mapping_files"
OUTPUT_DIR = BASE_DIR / "generated/salivadb"
LOG_DIR = BASE_DIR / "logs"

OUTPUT_FILE = OUTPUT_DIR / "biomarkers_salivadb.tsv"

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)
HGNC_FILE = MAPPING_DIR / "hgnc_complete_set.txt"

NCBI_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
NCBI_FILE = MAPPING_DIR / "gene_info.gz"

DOID_URL = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
DOID_FILE = MAPPING_DIR / "doid.obo"

MIRBASE_URL = "https://mirbase.org/download/CURRENT/id_mapping.tsv"
MIRBASE_FILE = MAPPING_DIR / "mirbase_id_mapping.tsv"

OLS_API = "https://www.ebi.ac.uk/ols4/api/search"
PUBMED_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

OUTPUT_COLUMNS = [
    "biomarker_index",
    "component_index",
    "entity_index",
    "assessed_biomarker_entity",
    "assessed_biomarker_entity_id",
    "assessed_entity_type",
    "condition",
    "condition_id",
    "best_biomarker_role",
    "specimen",
    "specimen_id",
    "loinc_code",
    "evidence_source",
    "evidence",
    "biomarker",
    "exposure_agent",
    "exposure_agent_id",
    "tag",
    "biomarker_controlled_vocab",
    "change_type_vocab",
    "aspect_type_vocab",
    "mod_type_vocab",
    "legacy_biomarker_id",
]

SPECIMEN_MAP = {
    "na": "",
    "parotid": "parotid gland",
    "salivary gland": "saliva-secreting gland",
}
SPECIMEN_ID_MAP = {
    "parotid gland": "UBERON:0001831",
    "saliva": "UBERON:0001836",
    "saliva-secreting gland": "UBERON:0001044",
    "tongue": "UBERON:0001723",
}

CHANGE_TYPE_MAP = {
    "decrease": "decreased",
    "downregulated": "decreased",
    "downregulated in this study; otherwise known to be up": "decreased",
    "c0-expression upregulatedregulated": "increased",
    "increase": "increased",
    "upregulated": "increased",
}

GENE_FUSION_PATTERN = re.compile(
    r"^[A-Za-z]{3,}[0-9][A-Za-z]{3,}[0-9](\s+\S+)?$"
)

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

LOG_DIR.mkdir(parents=True, exist_ok=True)
_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
_log_file = LOG_DIR / f"generate_salivadb_{_timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.FileHandler(_log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)


# ===========================================================================
# Mapping table loaders
# ===========================================================================


def download_if_missing(url: str, dest: Path, decompress_gz: bool = False) -> None:
    """Download *url* to *dest* only if *dest* does not already exist."""
    if dest.exists():
        logger.info("Mapping file already present, skipping download: %s", dest)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s → %s", url, dest)
    try:
        resp = requests.get(url, timeout=120, stream=True)
        resp.raise_for_status()
        if decompress_gz:
            raw = gzip.decompress(resp.content)
            dest.write_bytes(raw)
        else:
            dest.write_bytes(resp.content)
        logger.info("Download complete: %s (%d bytes)", dest, dest.stat().st_size)
    except requests.RequestException:
        logger.exception("Failed to download %s", url)
        raise


def load_hgnc() -> dict[str, str]:
    """Return {symbol → 'HGNC:<id>'} from the HGNC complete set."""
    download_if_missing(HGNC_URL, HGNC_FILE)
    logger.info("Loading HGNC mapping from %s", HGNC_FILE)
    df = pd.read_csv(HGNC_FILE, sep="\t", dtype=str, low_memory=False)
    mapping: dict[str, str] = {}
    for _, row in df.iterrows():
        sym = str(row.get("symbol", "")).strip()
        hgnc_id = str(row.get("hgnc_id", "")).strip()
        if sym and hgnc_id and hgnc_id.lower() != "nan":
            mapping[sym.upper()] = hgnc_id  # hgnc_id already "HGNC:XXXX"
    logger.info("HGNC: loaded %d symbol→ID entries", len(mapping))
    return mapping


def load_ncbi_gene() -> dict[str, str]:
    """Return {symbol → 'NCBI:<gene_id>'} from NCBI gene_info (Homo sapiens only)."""
    download_if_missing(NCBI_URL, NCBI_FILE, decompress_gz=True)
    logger.info("Loading NCBI gene mapping from %s", NCBI_FILE)
    df = pd.read_csv(
        NCBI_FILE,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    # gene_info has a '#tax_id' column; keep only human (9606)
    tax_col = df.columns[0]
    human = df[df[tax_col] == "9606"]
    mapping: dict[str, str] = {}
    for _, row in human.iterrows():
        sym = str(row.get("Symbol", "")).strip()
        gene_id = str(row.get("GeneID", "")).strip()
        if sym and gene_id and gene_id.lower() != "nan":
            mapping[sym.upper()] = f"NCBI:{gene_id}"
        # Also add synonyms
        for syn in str(row.get("Synonyms", "")).split("|"):
            syn = syn.strip()
            if syn and syn != "-":
                mapping.setdefault(syn.upper(), f"NCBI:{gene_id}")
    logger.info("NCBI: loaded %d symbol→ID entries", len(mapping))
    return mapping


def load_doid() -> dict[str, str]:
    """
    Parse doid.obo and return {lowercase_name → 'DOID:XXXXXXX'}.
    Also maps synonyms.
    """
    download_if_missing(DOID_URL, DOID_FILE)
    logger.info("Parsing DOID OBO from %s", DOID_FILE)
    name_to_id: dict[str, str] = {}
    current_id: Optional[str] = None
    current_name: Optional[str] = None
    synonyms: list[str] = []
    is_obsolete = False

    with open(DOID_FILE, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                # Flush previous term
                if current_id and current_name and not is_obsolete:
                    name_to_id[current_name.lower()] = current_id
                    for syn in synonyms:
                        name_to_id.setdefault(syn.lower(), current_id)
                current_id = None
                current_name = None
                synonyms = []
                is_obsolete = False
            elif line.startswith("id: "):
                current_id = line[4:].strip()
            elif line.startswith("name: "):
                current_name = line[6:].strip()
            elif line.startswith("synonym: "):
                m = re.search(r'"([^"]+)"', line)
                if m:
                    synonyms.append(m.group(1))
            elif line.startswith("is_obsolete: true"):
                is_obsolete = True
    # flush last term
    if current_id and current_name and not is_obsolete:
        name_to_id[current_name.lower()] = current_id

    logger.info("DOID: loaded %d name→ID entries (inc. synonyms)", len(name_to_id))
    return name_to_id


def load_mirbase() -> dict[str, str]:
    """Return {mature_id_or_stem_loop_name → 'MRB:<accession>'}."""
    download_if_missing(MIRBASE_URL, MIRBASE_FILE)
    logger.info("Loading miRBase mapping from %s", MIRBASE_FILE)
    # id_mapping.tsv columns: accession, name, type, mature/hairpin flag
    df = pd.read_csv(MIRBASE_FILE, sep="\t", header=None, dtype=str)
    mapping: dict[str, str] = {}
    for _, row in df.iterrows():
        acc = str(row.iloc[0]).strip()
        name = str(row.iloc[1]).strip()
        if acc and name:
            mapping[name.upper()] = f"MRB:{acc}"
    logger.info("miRBase: loaded %d name→accession entries", len(mapping))
    return mapping


# ===========================================================================
# ID look-up helpers
# ===========================================================================


def map_gene_symbol(
    symbol: str,
    hgnc: dict[str, str],
    ncbi: dict[str, str],
) -> str:
    """Map a gene symbol to HGNC or NCBI ID; return empty string if not found."""
    key = symbol.strip().upper()
    if key in hgnc:
        return hgnc[key]
    if key in ncbi:
        return ncbi[key]
    logger.warning("Gene symbol not mapped: %s", symbol)
    return ""


def map_condition_doid(
    condition: str,
    doid: dict[str, str],
) -> str:
    """
    Attempt to find a DOID for *condition* using exact match, then fuzzy
    (substring) match, then OLS API fall-back.
    Returns 'DOID:XXXXXXX' or ''.
    """
    key = condition.lower().strip()
    if key in doid:
        return doid[key]

    # Fuzzy: look for the condition as a substring key in the dict
    for name, did in doid.items():
        if key == name:
            return did
    # Partial / substring match – log it
    for name, did in doid.items():
        if key in name or name in key:
            logger.info("Fuzzy DOID match: '%s' → '%s' (%s)", condition, name, did)
            return did

    # OLS fall-back
    return _ols_lookup(condition)


def _ols_lookup(condition: str) -> str:
    """Query OLS4 for a DOID for *condition*. Returns ID string or ''."""
    try:
        resp = requests.get(
            OLS_API,
            params={
                "q": condition,
                "ontology": "doid",
                "fieldList": "id,label",
                "rows": 1,
            },
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        docs = data.get("response", {}).get("docs", [])
        if docs:
            ols_id = docs[0].get("id", "")
            label = docs[0].get("label", "")
            logger.info("OLS fallback: '%s' → '%s' (%s)", condition, label, ols_id)
            # OLS returns IDs like "DOID:12345" or "http://…"
            m = re.search(r"DOID[_:](\d+)", ols_id)
            if m:
                return f"DOID:{m.group(1)}"
        logger.warning("OLS: no DOID found for '%s'; needs human review", condition)
    except requests.RequestException:
        logger.exception("OLS lookup failed for '%s'", condition)
    return ""


def map_pubmed_via_api(author_string: str, year: Optional[str]) -> str:
    """
    Given an author string (e.g. "Sugimoto., Masahiro et al.") and an optional
    year, attempt to look up the PMID via the NCBI ESearch API.
    Returns 'PubMed:<pmid>' or the original *author_string* on failure.
    """
    query = author_string
    if year and str(year).strip().isdigit():
        query += f"[Author] AND {year.strip()}[PDAT]"
    else:
        query += "[Author]"
    try:
        resp = requests.get(
            PUBMED_API,
            params={
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retmax": 1,
            },
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            pmid = ids[0]
            logger.info("PubMed resolved '%s' (%s) → %s", author_string, year, pmid)
            return f"PubMed:{pmid}"
        logger.warning(
            "PubMed: could not resolve '%s' (year=%s); keeping raw string",
            author_string,
            year,
        )
    except requests.RequestException:
        logger.exception("PubMed API error for '%s'", author_string)
    return author_string


def map_mirna_id(
    mirna_id: str,
    mirbase: dict[str, str],
) -> str:
    """Map a miRNA ID to 'MRB:<accession>'."""
    key = mirna_id.strip().upper()
    if key in mirbase:
        return mirbase[key]
    # Try stem-loop (e.g. strip the -5p/-3p suffix)
    stem = re.sub(r"[-_][35]P$", "", key)
    if stem in mirbase:
        logger.info("miRNA stem-loop fallback: '%s' → '%s'", mirna_id, mirbase[stem])
        return mirbase[stem]
    logger.warning("miRNA ID not mapped: %s", mirna_id)
    return mirna_id


# ===========================================================================
# Field transformation helpers
# ===========================================================================


def clean_biomarker_name(name: str) -> str:
    """Remove quotes, leading/trailing underscores, square-bracket suffixes."""
    name = str(name)
    name = name.replace('"', "")
    name = re.sub(r"\s*\[.*?\]", "", name)
    name = name.strip("_").strip()
    return name


def clean_protein_name(name: str) -> str:
    """Remove all parenthetical groups from protein name."""
    name = str(name)
    name = name.replace('"', "")
    name = re.sub(r"\s*\(.*?\)", "", name)
    # Remove any [bracket] blocks too
    name = re.sub(r"\s*\[.*?\]", "", name)
    name = name.strip("_").strip()
    return name


def clean_condition(raw: str) -> str:
    """Lowercase, remove quotes, remove parenthetical content."""
    s = str(raw).replace('"', "")
    s = re.sub(r"\s*\(.*?\)", "", s)
    return s.strip().lower()


def map_specimen(raw: str) -> str:
    """Apply specimen mapping rules."""
    s = str(raw).strip().lower()
    if s in SPECIMEN_MAP:
        return SPECIMEN_MAP[s]
    # any string where the last word is "saliva"
    if s.split()[-1] == "saliva":
        return "saliva"
    return s


def map_change_type(raw: str) -> str:
    """Lowercase and apply CHANGE_TYPE_MAP; return '' for unmapped."""
    s = str(raw).strip().lower()
    return CHANGE_TYPE_MAP.get(s, s)


def parse_evidence_source(pmid_val: str, year_val: Optional[str]) -> str:
    """
    Return 'PubMed:<id>' for numeric PMIDs; resolve author strings via API.
    """
    val = str(pmid_val).replace('"', "").strip()
    if val.isdigit():
        return f"PubMed:{val}"
    return map_pubmed_via_api(val, year_val)


def split_biomarker_role(raw: str) -> list[str]:
    """
    Split on '/' and return cleaned, lowercased roles.
    E.g. "diagnostic/ prognostic" → ["diagnostic", "prognostic"]
    """
    s = str(raw).strip().lower()
    if s == "na":
        return [""]
    parts = [p.strip() for p in re.split(r"/", s)]
    return [p for p in parts if p] or [""]


def build_biomarker_controlled_vocab(
    change_type: str,
    aspect_type: str,
    entity_name: str,
    entity_id: str,
) -> str:
    parts = [p for p in [change_type, aspect_type, "of", entity_name] if p]
    combined = " ".join(parts)
    if entity_id:
        combined += f"/{entity_id}"
    return combined.lower()


# ===========================================================================
# Source-specific extraction functions
# ===========================================================================


def process_gene(
    df: pd.DataFrame,
    hgnc: dict[str, str],
    ncbi: dict[str, str],
    doid: dict[str, str],
) -> pd.DataFrame:
    logger.info("Processing gene.tsv (%d rows before filtering)", len(df))
    df = df[df["Organism"].str.strip() == "Homo sapiens"].copy()
    df = df[~df["Disease_name"].str.strip().isin(["Healthy", "NA"])].copy()

    rows = []
    for _, src in df.iterrows():
        sal_id = str(src["SAL_ID"]).split("_")[-1]
        raw_id = str(src["ID"]).strip()
        is_fusion = bool(GENE_FUSION_PATTERN.match(raw_id))

        if is_fusion:
            change_type = "presence of"
            aspect_type = "sequence variation"
            # Skip regulation==NA check for fusions
            regulation_ok = True
        else:
            if str(src.get("Regulation", "")).strip() == "NA":
                continue
            change_type = map_change_type(str(src.get("Regulation", "")))
            aspect_type = "expression"
            regulation_ok = True

        if not regulation_ok:
            continue

        condition = clean_condition(str(src["Disease_name"]))
        if condition in ("healthy", "na", ""):
            continue

        condition_id = map_condition_doid(condition, doid)
        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(
            str(src.get("PMID", "")), str(src.get("Year_of_Publication", ""))
        )
        entity_type = str(src.get("Biomarker", "")).strip().lower()
        base_name = clean_biomarker_name(str(src["Biomarker_Name"]))

        if is_fusion:
            # Split fusion: gene A is first N chars before last all-uppercase run
            # We follow the regex: everything up to the second uppercase block
            m = re.match(
                r"^([A-Za-z0-9]+?)([A-Z0-9]{3,}[0-9](\s+\S+)?)$", raw_id
            )
            if m:
                gene_a_sym = m.group(1)
                gene_b_raw = m.group(2).split()[0]
            else:
                # fallback: split at midpoint
                mid = len(raw_id) // 2
                gene_a_sym = raw_id[:mid]
                gene_b_raw = raw_id[mid:]

            id_a = map_gene_symbol(gene_a_sym, hgnc, ncbi)
            id_b = map_gene_symbol(gene_b_raw, hgnc, ncbi)
            gene_names = [gene_a_sym, gene_b_raw]
            gene_ids = [id_a, id_b]
            entity_indices = [1, 2]
        else:
            gene_names = [base_name]
            gene_ids = [map_gene_symbol(raw_id, hgnc, ncbi)]
            entity_indices = [1]

        for role in split_biomarker_role(str(src.get("Biomarker_Type", ""))):
            for gene_name, gene_id, eidx in zip(gene_names, gene_ids, entity_indices):
                bvc = build_biomarker_controlled_vocab(
                    change_type, aspect_type, gene_name, gene_id
                )
                rows.append(
                    {
                        "biomarker_index": sal_id,
                        "component_index": "1",
                        "entity_index": str(eidx),
                        "assessed_biomarker_entity": gene_name,
                        "assessed_biomarker_entity_id": gene_id,
                        "assessed_entity_type": entity_type,
                        "condition": condition,
                        "condition_id": condition_id,
                        "best_biomarker_role": role,
                        "specimen": specimen,
                        "specimen_id": specimen_id,
                        "loinc_code": "",
                        "evidence_source": evidence_source,
                        "evidence": "",
                        "biomarker": bvc,
                        "exposure_agent": "",
                        "exposure_agent_id": "",
                        "tag": "",
                        "biomarker_controlled_vocab": bvc,
                        "change_type_vocab": change_type,
                        "aspect_type_vocab": aspect_type,
                        "mod_type_vocab": "not_specified",
                        "legacy_biomarker_id": "no_legacy_biomarker_id",
                    }
                )

    logger.info("gene.tsv produced %d output rows", len(rows))
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def process_metabolite(
    df: pd.DataFrame,
    doid: dict[str, str],
) -> pd.DataFrame:
    logger.info("Processing metabolite.tsv (%d rows before filtering)", len(df))
    df = df[df["Organism"].str.strip() == "Homo sapiens"].copy()
    df = df[~df["Disease_name"].str.strip().isin(["Healthy", "NA"])].copy()
    df = df[df["Regulation"].str.strip() != "NA"].copy()

    rows = []
    for _, src in df.iterrows():
        sal_id = str(src["SAL_ID"]).split("_")[-1]
        condition = clean_condition(str(src["Disease_name"]))
        if condition in ("healthy", "na", ""):
            continue

        condition_id = map_condition_doid(condition, doid)
        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(
            str(src.get("PMID", "")), str(src.get("Year_of_Publication", ""))
        )
        entity_type = str(src.get("Biomarker", "")).strip().lower()
        entity_name = clean_biomarker_name(str(src["Biomarker_Name"]))
        change_type = map_change_type(str(src.get("Regulation", "")))
        aspect_type = "level"

        # ID column: strip quotes; multiple IDs → multiple rows
        raw_ids = str(src["ID"]).replace('"', "").strip()
        ids = [i.strip() for i in raw_ids.split(",") if i.strip()]
        if not ids:
            ids = [""]

        for pid in ids:
            entity_id = f"PCCID:{pid}" if pid else ""
            bvc = build_biomarker_controlled_vocab(
                change_type, aspect_type, entity_name, entity_id
            )
            for role in split_biomarker_role(str(src.get("Biomarker_Type", ""))):
                rows.append(
                    {
                        "biomarker_index": sal_id,
                        "component_index": "1",
                        "entity_index": "1",
                        "assessed_biomarker_entity": entity_name,
                        "assessed_biomarker_entity_id": entity_id,
                        "assessed_entity_type": entity_type,
                        "condition": condition,
                        "condition_id": condition_id,
                        "best_biomarker_role": role,
                        "specimen": specimen,
                        "specimen_id": specimen_id,
                        "loinc_code": "",
                        "evidence_source": evidence_source,
                        "evidence": "",
                        "biomarker": bvc,
                        "exposure_agent": "",
                        "exposure_agent_id": "",
                        "tag": "",
                        "biomarker_controlled_vocab": bvc,
                        "change_type_vocab": change_type,
                        "aspect_type_vocab": aspect_type,
                        "mod_type_vocab": "not_specified",
                        "legacy_biomarker_id": "no_legacy_biomarker_id",
                    }
                )

    logger.info("metabolite.tsv produced %d output rows", len(rows))
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def process_mirna(
    df: pd.DataFrame,
    mirbase: dict[str, str],
    doid: dict[str, str],
) -> pd.DataFrame:
    logger.info("Processing mirna.tsv (%d rows before filtering)", len(df))
    df = df[df["Organism"].str.strip() == "Homo sapiens"].copy()
    df = df[~df["Disease_name"].str.strip().isin(["Healthy", "NA"])].copy()
    df = df[df["Regulation"].str.strip() != "NA"].copy()

    rows = []
    for _, src in df.iterrows():
        sal_id = str(src["SAL_ID"]).split("_")[-1]
        condition = clean_condition(str(src["Disease_name"]))
        if condition in ("healthy", "na", ""):
            continue

        condition_id = map_condition_doid(condition, doid)
        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(
            str(src.get("PMID", "")), str(src.get("Year_of_Publication", ""))
        )
        entity_type = str(src.get("Biomarker", "")).strip()  # keep "miRNA" as-is
        entity_name = clean_biomarker_name(str(src["Biomarker_Name"]))
        change_type = map_change_type(str(src.get("Regulation", "")))
        aspect_type = "expression"

        raw_id = str(src["ID"]).strip()
        entity_id = map_mirna_id(raw_id, mirbase)
        bvc = build_biomarker_controlled_vocab(
            change_type, aspect_type, entity_name, entity_id
        )

        for role in split_biomarker_role(str(src.get("Biomarker_Type", ""))):
            rows.append(
                {
                    "biomarker_index": sal_id,
                    "component_index": "1",
                    "entity_index": "1",
                    "assessed_biomarker_entity": entity_name,
                    "assessed_biomarker_entity_id": entity_id,
                    "assessed_entity_type": entity_type,
                    "condition": condition,
                    "condition_id": condition_id,
                    "best_biomarker_role": role,
                    "specimen": specimen,
                    "specimen_id": specimen_id,
                    "loinc_code": "",
                    "evidence_source": evidence_source,
                    "evidence": "",
                    "biomarker": bvc,
                    "exposure_agent": "",
                    "exposure_agent_id": "",
                    "tag": "",
                    "biomarker_controlled_vocab": bvc,
                    "change_type_vocab": change_type,
                    "aspect_type_vocab": aspect_type,
                    "mod_type_vocab": "not_specified",
                    "legacy_biomarker_id": "no_legacy_biomarker_id",
                }
            )

    logger.info("mirna.tsv produced %d output rows", len(rows))
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def process_protein(
    df: pd.DataFrame,
    doid: dict[str, str],
) -> pd.DataFrame:
    logger.info("Processing protein.tsv (%d rows before filtering)", len(df))
    df = df[df["Organism"].str.strip() == "Homo sapiens"].copy()
    df = df[~df["Disease_name"].str.strip().isin(["Healthy", "NA"])].copy()
    df = df[df["Regulation"].str.strip() != "NA"].copy()

    rows = []
    for _, src in df.iterrows():
        sal_id = str(src["SAL_ID"]).split("_")[-1]
        condition = clean_condition(str(src["Disease_name"]))
        if condition in ("healthy", "na", ""):
            continue

        condition_id = map_condition_doid(condition, doid)
        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(
            str(src.get("PMID", "")), str(src.get("Year_of_Publication", ""))
        )
        entity_type = str(src.get("Biomarker", "")).strip().lower()
        entity_name = clean_protein_name(str(src["Biomarker_Name"]))
        change_type = map_change_type(str(src.get("Regulation", "")))
        aspect_type = "level"

        raw_id = str(src["ID"]).strip()
        entity_id = f"UPKB:{raw_id}" if raw_id else ""
        bvc = build_biomarker_controlled_vocab(
            change_type, aspect_type, entity_name, entity_id
        )

        for role in split_biomarker_role(str(src.get("Biomarker_Type", ""))):
            rows.append(
                {
                    "biomarker_index": sal_id,
                    "component_index": "1",
                    "entity_index": "1",
                    "assessed_biomarker_entity": entity_name,
                    "assessed_biomarker_entity_id": entity_id,
                    "assessed_entity_type": entity_type,
                    "condition": condition,
                    "condition_id": condition_id,
                    "best_biomarker_role": role,
                    "specimen": specimen,
                    "specimen_id": specimen_id,
                    "loinc_code": "",
                    "evidence_source": evidence_source,
                    "evidence": "",
                    "biomarker": bvc,
                    "exposure_agent": "",
                    "exposure_agent_id": "",
                    "tag": "",
                    "biomarker_controlled_vocab": bvc,
                    "change_type_vocab": change_type,
                    "aspect_type_vocab": aspect_type,
                    "mod_type_vocab": "not_specified",
                    "legacy_biomarker_id": "no_legacy_biomarker_id",
                }
            )

    logger.info("protein.tsv produced %d output rows", len(rows))
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


# ===========================================================================
# Main
# ===========================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate SalivaDB biomarkers dataset."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=INPUT_DIR,
        help=f"Directory containing the four source TSVs (default: {INPUT_DIR})",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=OUTPUT_FILE,
        help=f"Path for the output TSV (default: {OUTPUT_FILE})",
    )
    parser.add_argument(
        "--mapping-dir",
        type=Path,
        default=MAPPING_DIR,
        help=f"Directory for mapping files (default: {MAPPING_DIR})",
    )
    args = parser.parse_args()

    # Allow CLI overrides to propagate to module-level globals used by loaders
    global INPUT_DIR, OUTPUT_FILE, MAPPING_DIR
    global HGNC_FILE, NCBI_FILE, DOID_FILE, MIRBASE_FILE
    INPUT_DIR = args.input_dir
    OUTPUT_FILE = args.output_file
    MAPPING_DIR = args.mapping_dir
    HGNC_FILE = MAPPING_DIR / "hgnc_complete_set.txt"
    NCBI_FILE = MAPPING_DIR / "gene_info.gz"
    DOID_FILE = MAPPING_DIR / "doid.obo"
    MIRBASE_FILE = MAPPING_DIR / "mirbase_id_mapping.tsv"

    logger.info("=== SalivaDB extraction started ===")
    logger.info("Input directory : %s", INPUT_DIR)
    logger.info("Output file     : %s", OUTPUT_FILE)
    logger.info("Mapping dir     : %s", MAPPING_DIR)

    # ------------------------------------------------------------------
    # Load mapping tables
    # ------------------------------------------------------------------
    hgnc = load_hgnc()
    ncbi = load_ncbi_gene()
    doid = load_doid()
    mirbase = load_mirbase()

    # ------------------------------------------------------------------
    # Load source files
    # ------------------------------------------------------------------
    def read_tsv(name: str) -> pd.DataFrame:
        path = INPUT_DIR / name
        logger.info("Reading %s", path)
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")

    gene_df = read_tsv("gene.tsv")
    metabolite_df = read_tsv("metabolite.tsv")
    mirna_df = read_tsv("mirna.tsv")
    protein_df = read_tsv("protein.tsv")

    # ------------------------------------------------------------------
    # Process each source
    # ------------------------------------------------------------------
    frames = [
        process_gene(gene_df, hgnc, ncbi, doid),
        process_metabolite(metabolite_df, doid),
        process_mirna(mirna_df, mirbase, doid),
        process_protein(protein_df, doid),
    ]

    result = pd.concat(frames, ignore_index=True)
    logger.info("Total output rows: %d", len(result))

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT_FILE, sep="\t", index=False)
    logger.info("Output written to %s", OUTPUT_FILE)
    logger.info("=== SalivaDB extraction complete ===")


if __name__ == "__main__":
    main()
