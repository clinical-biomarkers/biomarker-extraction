from __future__ import annotations

import argparse
import dataclasses
import gzip
import io
import json
import logging
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from dotenv import load_dotenv
from pathlib import Path
from typing import Optional

from scripts.utils.log_utils import setup_logging

import pandas as pd
import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------

BASE_DIR = Path("/data/shared/repos/extraction")
DEFAULT_INPUT_DIR = BASE_DIR / "downloads/salivadb"
DEFAULT_MAPPING_DIR = BASE_DIR / "downloads/mapping_files"
DEFAULT_OUTPUT_FILE = BASE_DIR / "generated/salivadb/biomarkers_salivadb.tsv"
DEFAULT_LOG_DIR = BASE_DIR / "logs"
DEFAULT_CURATED_MAPPING_FILE = BASE_DIR / "downloads/mapping_files/curated_condition_mappings.json"

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)

DOID_URL = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
MIRBASE_URL = "https://mirbase.org/download/CURRENT/id_mapping.tsv"

OLS_API = "https://www.ebi.ac.uk/ols4/api/search"

NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

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

# ===========================================================================
# Curated condition mappings
# ===========================================================================


@dataclasses.dataclass(frozen=True)
class ConditionMapping:
    """A single resolved (condition, condition_id) pair from the curated file."""
    condition: str
    condition_id: str


@dataclasses.dataclass(frozen=True)
class ExposureAgentMapping:
    """An exposure-agent resolution from the curated file."""
    exposure_agent: str
    exposure_agent_id: str


@dataclasses.dataclass
class CuratedMappings:
    """
    Parsed contents of curated_condition_mappings.json.

    Attributes
    ----------
    exclude : set[str]
        Lowercase disease names that should be dropped entirely.
    conditions : dict[str, list[ConditionMapping]]
        Lowercase disease name → one or more canonical (condition, condition_id)
        pairs.  A one-to-many entry produces multiple output rows (explode).
    exposure_agents : dict[str, ExposureAgentMapping]
        Lowercase disease name → exposure agent + ID.  These rows keep an empty
        condition/condition_id and fill the exposure_agent columns instead.
    """

    exclude: set[str]
    conditions: dict[str, list[ConditionMapping]]
    exposure_agents: dict[str, ExposureAgentMapping]


def load_curated_mappings(path: Path) -> CuratedMappings:
    """
    Load *path* (JSON) and return a :class:`CuratedMappings` instance.

    Expected JSON structure::

        {
          "exclude": ["physical frailty", ...],
          "conditions": {
            "salivary gland tumor": [
              {"condition": "salivary gland cancer", "condition_id": "DOID:8850"}
            ],
            "alcohol and tobacco dependent": [
              {"condition": "alcohol dependence",  "condition_id": "DOID:0050741"},
              {"condition": "nicotine dependence", "condition_id": "DOID:0050742"}
            ]
          },
          "exposure_agents": {
            "prenatal arsenic exposure": {
              "exposure_agent": "arsenic",
              "exposure_agent_id": "PCCID:5359596"
            }
          }
        }

    All keys are stored and matched in lowercase.
    """
    logger.info("Loading curated condition mappings from %s", path)
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)

    exclude: set[str] = {s.lower().strip() for s in raw.get("exclude", [])}

    conditions: dict[str, list[ConditionMapping]] = {}
    for disease_name, entries in raw.get("conditions", {}).items():
        key = disease_name.lower().strip()
        conditions[key] = [
            ConditionMapping(
                condition=e["condition"].strip(),
                condition_id=e["condition_id"].strip(),
            )
            for e in entries
        ]

    exposure_agents: dict[str, ExposureAgentMapping] = {}
    for disease_name, entry in raw.get("exposure_agents", {}).items():
        key = disease_name.lower().strip()
        exposure_agents[key] = ExposureAgentMapping(
            exposure_agent=entry["exposure_agent"].strip(),
            exposure_agent_id=entry["exposure_agent_id"].strip(),
        )

    logger.info(
        "Curated mappings loaded: %d excludes, %d condition remaps, %d exposure agents",
        len(exclude),
        len(conditions),
        len(exposure_agents),
    )
    return CuratedMappings(
        exclude=exclude,
        conditions=conditions,
        exposure_agents=exposure_agents,
    )


def resolve_condition(
    raw_disease_name: str,
    curated: CuratedMappings,
    doid: dict[str, str],
) -> list[tuple[str, str, str, str]]:
    """
    Resolve a raw disease name to one or more
    ``(condition, condition_id, exposure_agent, exposure_agent_id)`` tuples.

    Resolution order
    ----------------
    1. Curated exclude list → returns ``[]`` (row is dropped by the caller).
    2. Curated exposure-agent mapping → one tuple with empty condition fields.
    3. Curated condition mapping → one tuple per entry (explode for one-to-many).
    4. DOID exact / fuzzy / OLS fallback → one tuple, exposure fields empty.

    Parameters
    ----------
    raw_disease_name:
        The value straight from ``clean_condition()``.
    curated:
        Loaded curated mappings.
    doid:
        DOID name→ID dict from :func:`load_doid`.

    Returns
    -------
    list of (condition, condition_id, exposure_agent, exposure_agent_id)
        Empty list means the row should be skipped.
    """
    key = raw_disease_name.lower().strip()

    # 1. Exclude
    if key in curated.exclude:
        logger.info("Curated exclude: '%s'", raw_disease_name)
        return []

    # 2. Exposure agent
    if key in curated.exposure_agents:
        ea = curated.exposure_agents[key]
        logger.info(
            "Curated exposure agent: '%s' → %s (%s)",
            raw_disease_name,
            ea.exposure_agent,
            ea.exposure_agent_id,
        )
        return [("", "", ea.exposure_agent, ea.exposure_agent_id)]

    # 3. Curated condition remap (may be one-to-many)
    if key in curated.conditions:
        entries = curated.conditions[key]
        logger.info(
            "Curated condition remap: '%s' → %d entr%s",
            raw_disease_name,
            len(entries),
            "y" if len(entries) == 1 else "ies",
        )
        return [(cm.condition, cm.condition_id, "", "") for cm in entries]

    # 4. DOID fallback (existing fuzzy + OLS logic)
    condition_id = map_condition_doid(raw_disease_name, doid)
    return [(raw_disease_name, condition_id, "", "")]


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


def load_hgnc(hgnc_file: Path) -> dict[str, str]:
    """Return {symbol → 'HGNC:<id>'} from the HGNC complete set."""
    download_if_missing(HGNC_URL, hgnc_file)
    logger.info("Loading HGNC mapping from %s", hgnc_file)
    df = pd.read_csv(
        hgnc_file,
        sep="\t",
        dtype=str,
        low_memory=False,
        usecols=["hgnc_id", "symbol"],
    )
    df = df.dropna(subset=["hgnc_id", "symbol"])
    mapping = dict(zip(df["symbol"].str.strip().str.upper(), df["hgnc_id"].str.strip()))
    logger.info("HGNC: loaded %d symbol→ID entries", len(mapping))
    return mapping


def load_doid(doid_file: Path) -> dict[str, str]:
    """
    Parse doid.obo and return {lowercase_name → 'DOID:XXXXXXX'}.
    Also maps synonyms.
    """
    download_if_missing(DOID_URL, doid_file)
    logger.info("Parsing DOID OBO from %s", doid_file)
    name_to_id: dict[str, str] = {}
    current_id: Optional[str] = None
    current_name: Optional[str] = None
    synonyms: list[str] = []
    is_obsolete = False

    with open(doid_file, encoding="utf-8") as fh:
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


def load_mirbase(mirbase_file: Path) -> dict[str, str]:
    """Return {mature_id_or_stem_loop_name → 'MRB:<accession>'}."""
    download_if_missing(MIRBASE_URL, mirbase_file)
    logger.info("Loading miRBase mapping from %s", mirbase_file)
    # id_mapping.tsv columns: accession, name, type, mature/hairpin flag
    df = pd.read_csv(mirbase_file, sep="\t", header=None, dtype=str)
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

_ncbi_gene_cache: dict[str, str] = {}

def _ncbi_gene_lookup(symbol: str, api_key: str = "") -> str:
    key = symbol.strip().upper()
    if key in _ncbi_gene_cache:
        return _ncbi_gene_cache[key]

    params_base = {"api_key": api_key} if api_key else {}
    # Stay under NCBI's rate limit: 10 req/s with key, 3 req/s without.
    base_delay = 0.11 if api_key else 0.4

    for attempt in range(4):  # initial attempt + up to 3 retries
        sleep_time = base_delay if attempt == 0 else base_delay * (2 ** attempt)
        time.sleep(sleep_time)
        try:
            resp = requests.get(
                NCBI_ESEARCH_URL,
                params={
                    **params_base,
                    "db": "gene",
                    "term": f'{key}[sym] AND "Homo sapiens"[Organism]',
                    "retmode": "json",
                    "retmax": 1,
                },
                timeout=30,
            )
            if resp.status_code == 429:
                logger.warning(
                    "NCBI 429 for '%s'; backing off %.1fs (attempt %d/4)",
                    symbol, sleep_time * 2, attempt + 1,
                )
                continue  # sleep is doubled on the next iteration
            resp.raise_for_status()

            ids = resp.json().get("esearchresult", {}).get("idlist", [])
            if not ids:
                logger.warning("NCBI Gene: no result for symbol '%s'", symbol)
                _ncbi_gene_cache[key] = ""
                return ""

            ncbi_id = f"NCBI:{ids[0]}"
            logger.info("NCBI Gene fallback: '%s' → %s", symbol, ncbi_id)
            _ncbi_gene_cache[key] = ncbi_id
            return ncbi_id

        except requests.RequestException:
            logger.exception("NCBI Gene lookup failed for '%s'", symbol)
            _ncbi_gene_cache[key] = ""
            return ""

    logger.error("NCBI Gene: exhausted retries for '%s'; dropping", symbol)
    _ncbi_gene_cache[key] = ""
    return ""

def map_gene_symbol(
    symbol: str,
    hgnc: dict[str, str],
    ncbi_api_key: str = "",
) -> str:
    """Map a gene symbol to an HGNC ID; falls back to NCBI Gene API on miss."""
    key = symbol.strip().upper()
    if key in hgnc:
        return hgnc[key]
    # HGNC table miss → try NCBI Gene
    result = _ncbi_gene_lookup(key, api_key=ncbi_api_key)
    if not result:
        logger.warning("Gene symbol not mapped: %s", symbol)
    return result


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
    if s.split() and s.split()[-1] == "saliva":
        return "saliva"
    return s


def map_change_type(raw: str) -> str:
    """Lowercase and apply CHANGE_TYPE_MAP; return '' for unmapped."""
    s = str(raw).strip().lower()
    return CHANGE_TYPE_MAP.get(s, s)


def parse_evidence_source(pmid_val: str) -> Optional[str]:
    """
    Return 'PubMed:<id>' for numeric PMIDs; return None for non-numeric values.
    Rows with a None evidence_source are excluded by callers.
    """
    val = str(pmid_val).replace('"', "").strip()
    if val.isdigit():
        return f"PubMed:{val}"
    logger.info("Skipping row: non-numeric PMID '%s'", val)
    return None


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
    doid: dict[str, str],
    curated: CuratedMappings,
    ncbi_api_key: str = "",
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
            # Strip the " fusion" annotation before splitting
            fusion_id = re.sub(r'\s+fusion\b.*', '', raw_id, flags=re.IGNORECASE).strip()
            m = re.match(r'^([A-Za-z]+[0-9]+)([A-Za-z]+[0-9]+)$', fusion_id)
            if m:
                gene_a_sym = m.group(1)
                gene_b_raw = m.group(2)
            else:
                mid = len(fusion_id) // 2
                gene_a_sym = fusion_id[:mid]
                gene_b_raw = fusion_id[mid:]
            gene_names = [gene_a_sym, gene_b_raw]
            gene_ids = [map_gene_symbol(gene_a_sym, hgnc, ncbi_api_key),
                        map_gene_symbol(gene_b_raw, hgnc, ncbi_api_key)]
            entity_indices = [1, 2]
        else:
            if str(src.get("Regulation", "")).strip() == "NA":
                continue
            change_type = map_change_type(str(src.get("Regulation", "")))
            aspect_type = "expression"
            gene_names = [base_name]
            gene_ids = [map_gene_symbol(raw_id, hgnc, ncbi_api_key)]
            entity_indices = [1]

        for cond, cond_id, exp_agent, exp_agent_id in condition_tuples:
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
                            "condition": cond,
                            "condition_id": cond_id,
                            "best_biomarker_role": role,
                            "specimen": specimen,
                            "specimen_id": specimen_id,
                            "loinc_code": "",
                            "evidence_source": evidence_source,
                            "evidence": "",
                            "biomarker": bvc,
                            "exposure_agent": exp_agent,
                            "exposure_agent_id": exp_agent_id,
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
    curated: CuratedMappings,
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

        condition_tuples = resolve_condition(condition, curated, doid)
        if not condition_tuples:
            continue

        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(str(src.get("PMID", "")))
        if evidence_source is None:
            continue
        entity_type = str(src.get("Biomarker", "")).strip().lower()
        entity_name = clean_biomarker_name(str(src["Biomarker_Name"]))
        change_type = map_change_type(str(src.get("Regulation", "")))
        aspect_type = "level"

        # ID column: strip quotes; multiple IDs → multiple rows
        raw_ids = str(src["ID"]).replace('"', "").strip()
        ids = [i.strip() for i in raw_ids.split(",") if i.strip()]
        if not ids:
            ids = [""]

        for cond, cond_id, exp_agent, exp_agent_id in condition_tuples:
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
                            "condition": cond,
                            "condition_id": cond_id,
                            "best_biomarker_role": role,
                            "specimen": specimen,
                            "specimen_id": specimen_id,
                            "loinc_code": "",
                            "evidence_source": evidence_source,
                            "evidence": "",
                            "biomarker": bvc,
                            "exposure_agent": exp_agent,
                            "exposure_agent_id": exp_agent_id,
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
    curated: CuratedMappings,
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

        condition_tuples = resolve_condition(condition, curated, doid)
        if not condition_tuples:
            continue

        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(str(src.get("PMID", "")))
        if evidence_source is None:
            continue
        entity_type = str(src.get("Biomarker", "")).strip()  # keep "miRNA" as-is
        entity_name = clean_biomarker_name(str(src["Biomarker_Name"]))
        change_type = map_change_type(str(src.get("Regulation", "")))
        aspect_type = "expression"

        raw_id = str(src["ID"]).strip()
        entity_id = map_mirna_id(raw_id, mirbase)
        bvc = build_biomarker_controlled_vocab(
            change_type, aspect_type, entity_name, entity_id
        )

        for cond, cond_id, exp_agent, exp_agent_id in condition_tuples:
            for role in split_biomarker_role(str(src.get("Biomarker_Type", ""))):
                rows.append(
                    {
                        "biomarker_index": sal_id,
                        "component_index": "1",
                        "entity_index": "1",
                        "assessed_biomarker_entity": entity_name,
                        "assessed_biomarker_entity_id": entity_id,
                        "assessed_entity_type": entity_type,
                        "condition": cond,
                        "condition_id": cond_id,
                        "best_biomarker_role": role,
                        "specimen": specimen,
                        "specimen_id": specimen_id,
                        "loinc_code": "",
                        "evidence_source": evidence_source,
                        "evidence": "",
                        "biomarker": bvc,
                        "exposure_agent": exp_agent,
                        "exposure_agent_id": exp_agent_id,
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
    curated: CuratedMappings,
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

        condition_tuples = resolve_condition(condition, curated, doid)
        if not condition_tuples:
            continue

        specimen = map_specimen(str(src.get("Collection_site", "")))
        specimen_id = SPECIMEN_ID_MAP.get(specimen, "")
        evidence_source = parse_evidence_source(str(src.get("PMID", "")))
        if evidence_source is None:
            continue
        entity_type = str(src.get("Biomarker", "")).strip().lower()
        entity_name = clean_protein_name(str(src["Biomarker_Name"]))
        change_type = map_change_type(str(src.get("Regulation", "")))
        aspect_type = "level"

        raw_id = str(src["ID"]).strip()
        entity_id = f"UPKB:{raw_id}" if raw_id else ""
        bvc = build_biomarker_controlled_vocab(
            change_type, aspect_type, entity_name, entity_id
        )

        for cond, cond_id, exp_agent, exp_agent_id in condition_tuples:
            for role in split_biomarker_role(str(src.get("Biomarker_Type", ""))):
                rows.append(
                    {
                        "biomarker_index": sal_id,
                        "component_index": "1",
                        "entity_index": "1",
                        "assessed_biomarker_entity": entity_name,
                        "assessed_biomarker_entity_id": entity_id,
                        "assessed_entity_type": entity_type,
                        "condition": cond,
                        "condition_id": cond_id,
                        "best_biomarker_role": role,
                        "specimen": specimen,
                        "specimen_id": specimen_id,
                        "loinc_code": "",
                        "evidence_source": evidence_source,
                        "evidence": "",
                        "biomarker": bvc,
                        "exposure_agent": exp_agent,
                        "exposure_agent_id": exp_agent_id,
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
    load_dotenv()
    parser = argparse.ArgumentParser(
        description="Generate SalivaDB biomarkers dataset."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing the four source TSVs (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=DEFAULT_OUTPUT_FILE,
        help=f"Path for the output TSV (default: {DEFAULT_OUTPUT_FILE})",
    )
    parser.add_argument(
        "--mapping-dir",
        type=Path,
        default=DEFAULT_MAPPING_DIR,
        help=f"Directory for mapping files (default: {DEFAULT_MAPPING_DIR})",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help=f"Directory for log files (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--curated-mapping-file",
        type=Path,
        default=DEFAULT_CURATED_MAPPING_FILE,
        help=f"Path to curated condition mappings JSON (default: {DEFAULT_CURATED_MAPPING_FILE})",
    )
    parser.add_argument(
        "--ncbi-api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key (optional; raises rate limit from 3 to 10 req/s). "
            "Falls back to NCBI_API_KEY env var.",
    )
    args = parser.parse_args()

    setup_logging(args.log_dir)

    # Derive mapping file paths from the mapping directory
    hgnc_file = args.mapping_dir / "hgnc_complete_set.txt"
    doid_file = args.mapping_dir / "doid.obo"
    mirbase_file = args.mapping_dir / "mirbase_id_mapping.tsv"

    logger.info("=== SalivaDB extraction started ===")
    logger.info("Input directory : %s", args.input_dir)
    logger.info("Output file     : %s", args.output_file)
    logger.info("Mapping dir     : %s", args.mapping_dir)

    # ------------------------------------------------------------------
    # Load mapping tables
    # ------------------------------------------------------------------
    curated = load_curated_mappings(args.curated_mapping_file)
    hgnc = load_hgnc(hgnc_file)
    doid = load_doid(doid_file)
    mirbase = load_mirbase(mirbase_file)

    # ------------------------------------------------------------------
    # Load source files
    # ------------------------------------------------------------------
    def read_tsv(name: str) -> pd.DataFrame:
        path = args.input_dir / name
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
        process_gene(gene_df, hgnc, doid, curated, ncbi_api_key=args.ncbi_api_key),
        process_metabolite(metabolite_df, doid, curated),
        process_mirna(mirna_df, mirbase, doid, curated),
        process_protein(protein_df, doid, curated),
    ]

    result = pd.concat(frames, ignore_index=True)
    logger.info("Total output rows before ID filter: %d", len(result))

    before = len(result)
    result = result[result["assessed_biomarker_entity_id"] != ""].copy()
    logger.info(
        "Dropped %d rows with empty assessed_biomarker_entity_id; keeping %d",
        before - len(result),
        len(result),
    )

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_file, sep="\t", index=False)
    logger.info("Output written to %s", args.output_file)
    logger.info("=== SalivaDB extraction complete ===")


if __name__ == "__main__":
    main()
