# `generate_salivadb.py`

## Overview

`generate_salivadb.py` is the primary extraction and harmonisation script for the **SalivaDB** data. It reads four raw TSV files downloaded from the SalivaDB source database (genes, metabolites, miRNAs, and proteins), maps every entity to a controlled vocabulary identifier, resolves disease names to Disease Ontology (DOID), MONDO, or SNOMED terms, and writes a single harmonised biomarker table in TSV format suitable for downstream ingestion into BiomarkerKB.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Dependencies](#dependencies)
3. [CLI Reference](#cli-reference)
4. [Input Files](#input-files)
5. [Output Format](#output-format)
6. [Architecture Overview](#architecture-overview)
7. [Configuration & Constants](#configuration--constants)
8. [Mapping Tables](#mapping-tables)
9. [Curated Condition Mappings](#curated-condition-mappings)
10. [Source Processing Functions](#source-processing-functions)
11. [Field Transformation Helpers](#field-transformation-helpers)
12. [ID Lookup Functions](#id-lookup-functions)
13. [Logging & Diagnostics](#logging--diagnostics)
14. [External API Calls](#external-api-calls)
15. [Known Limitations & Caveats](#known-limitations--caveats)

---

## Quick Start

```bash
# Minimal invocation — uses all default paths
python generate_salivadb.py

# Specify custom paths
python generate_salivadb.py \
  --input-dir  /data/salivadb/raw \
  --output-file /data/salivadb/out/biomarkers_salivadb.tsv \
  --mapping-dir /data/salivadb/mappings \
  --curated-mapping-file /data/salivadb/curated_condition_mappings.json \
  --ncbi-api-key "$NCBI_API_KEY"
```

The script will auto-download any missing mapping files (HGNC, DOID, miRBase) on first run, so an internet connection is required unless all mapping files are pre-cached in `--mapping-dir`.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `pandas` | TSV I/O and tabular transformations |
| `requests` | HTTP downloads (mapping files, NCBI/OLS APIs) |
| `python-dotenv` | Load `NCBI_API_KEY` from a `.env` file |

All are standard pip-installable packages. The script also uses only stdlib modules: `argparse`, `dataclasses`, `gzip`, `json`, `logging`, `os`, `re`, `time`, `xml.etree.ElementTree`, `pathlib`, and `datetime`.

The internal helper `scripts.utils.log_utils.setup_logging` must be importable from the working directory.

---

## CLI Reference

| Argument | Default | Description |
|----------|---------|-------------|
| `--input-dir` | `/data/shared/repos/extraction/downloads/salivadb` | Directory containing the four source TSVs. |
| `--output-file` | `/data/shared/repos/extraction/generated/salivadb/biomarkers_salivadb.tsv` | Path for the harmonised output TSV. Created if it does not exist. |
| `--mapping-dir` | `/data/shared/repos/extraction/downloads/mapping_files` | Directory where downloaded mapping files are cached. |
| `--log-dir` | `/data/shared/repos/extraction/logs` | Directory for structured log files. |
| `--curated-mapping-file` | `<mapping-dir>/curated_condition_mappings.json` | Path to the hand-curated condition remapping JSON. |
| `--ncbi-api-key` | `$NCBI_API_KEY` env var | Optional NCBI API key; raises the NCBI rate limit from 3 req/s to 10 req/s. |

---

## Input Files

The script expects the following four TSV files inside `--input-dir`:

| File | Entity type | Key columns |
|------|-------------|-------------|
| `gene.tsv` | Protein-coding genes (and gene fusions) | `SAL_ID`, `ID`, `Biomarker_Name`, `Biomarker`, `Disease_name`, `Regulation`, `Collection_site`, `PMID`, `Biomarker_Type`, `Organism` |
| `metabolite.tsv` | Small-molecule metabolites | same as above, `ID` holds PubChem CID(s) (comma-separated) |
| `mirna.tsv` | microRNAs | same as above, `ID` holds the miRBase name (e.g. `hsa-miR-21-5p`) |
| `protein.tsv` | Proteins | same as above, `ID` holds a UniProtKB accession |

All files are read as `dtype=str` with `NaN` values replaced by empty strings.

---

## Output Format

The output is a tab-separated file with the following 23 columns:

| Column | Description |
|--------|-------------|
| `biomarker_index` | Numeric ID extracted from `SAL_ID` (the part after `_`) |
| `component_index` | Always `"1"` in the current implementation |
| `entity_index` | `"1"` for single-entity rows; `"1"` / `"2"` for gene fusion partners |
| `assessed_biomarker_entity` | Cleaned entity name |
| `assessed_biomarker_entity_id` | Ontology-prefixed ID (e.g. `HGNC:1100`, `UPKB:P04637`, `MRB:MIMAT0000076`, `PCCID:5793`) |
| `assessed_entity_type` | Lowercase entity type as recorded in the source (`gene`, `protein`, `mirna`, `metabolite`) |
| `condition` | Canonical disease name |
| `condition_id` | `DOID:XXXXXXX` identifier |
| `best_biomarker_role` | Biomarker role(s) (`diagnostic`, `prognostic`, etc.) — split on `/` for multi-role rows |
| `specimen` | Normalised specimen name |
| `specimen_id` | UBERON ID for the specimen |
| `loinc_code` | Empty in the current implementation |
| `evidence_source` | `PubMed:<PMID>` |
| `evidence` | Empty in the current implementation |
| `biomarker` | Biomarker controlled-vocab string (see below) |
| `exposure_agent` | Exposure agent name (populated from curated mappings; otherwise empty) |
| `exposure_agent_id` | Exposure agent ID (e.g. `PCCID:5359596`) |
| `tag` | Empty in the current implementation |
| `biomarker_controlled_vocab` | `"<change_type> <aspect_type> of <entity_name>/<entity_id>"` |
| `change_type_vocab` | Normalised change direction (`increased`, `decreased`, `presence of`) |
| `aspect_type_vocab` | `expression`, `level`, or `sequence variation` |
| `mod_type_vocab` | Always `"not_specified"` |
| `legacy_biomarker_id` | Always `"no_legacy_biomarker_id"` |

Rows where `assessed_biomarker_entity_id` is empty are dropped before writing.

---

## Architecture Overview

```
main()
 ├── load_curated_mappings()       # JSON → CuratedMappings dataclass
 ├── load_hgnc()                   # TSV → {symbol → HGNC:ID}
 ├── load_doid()                   # OBO → {name → DOID:ID}
 ├── load_mirbase()                # TSV + dat files → {NAME → MRB:accession}
 │
 ├── process_gene()                # gene.tsv → DataFrame
 ├── process_metabolite()          # metabolite.tsv → DataFrame
 ├── process_mirna()               # mirna.tsv → DataFrame
 ├── process_protein()             # protein.tsv → DataFrame
 │
 ├── pd.concat(frames)             # merge all entity types
 ├── drop rows with empty entity ID
 └── write TSV
```

Each `process_*` function follows the same pattern:
1. Filter to `Organism == "Homo sapiens"` and exclude `Healthy` / `NA` disease names.
2. For each row, call `resolve_condition()` to get one or more `(condition, condition_id, exposure_agent, exposure_agent_id)` tuples.
3. Map the entity identifier to its ontology ID.
4. Build the controlled-vocab string.
5. Explode multi-role and multi-condition rows into separate output rows.

---

## Configuration & Constants

### Default Paths

```python
BASE_DIR             = Path("/data/shared/repos/extraction")
DEFAULT_INPUT_DIR    = BASE_DIR / "downloads/salivadb"
DEFAULT_MAPPING_DIR  = BASE_DIR / "downloads/mapping_files"
DEFAULT_OUTPUT_FILE  = BASE_DIR / "generated/salivadb/biomarkers_salivadb.tsv"
DEFAULT_LOG_DIR      = BASE_DIR / "logs"
DEFAULT_CURATED_MAPPING_FILE = BASE_DIR / "downloads/mapping_files/curated_condition_mappings.json"
```

### Specimen Mapping

| Raw value | Normalised value | UBERON ID |
|-----------|-----------------|-----------|
| `na` | `""` | — |
| `parotid` | `parotid gland` | `UBERON:0001831` |
| `salivary gland` | `saliva-secreting gland` | `UBERON:0001044` |
| `saliva` (or `* saliva`) | `saliva` | `UBERON:0001836` |
| `tongue` | `tongue` | `UBERON:0001723` |

### Change-Type Mapping

Raw values are normalised to `increased` or `decreased`:

| Raw | Normalised |
|-----|-----------|
| `increase`, `upregulated`, `c0-expression upregulatedregulated` | `increased` |
| `decrease`, `downregulated`, `downregulated in this study; otherwise known to be up` | `decreased` |

Gene fusions use the special value `presence of`.

---

## Mapping Tables

### HGNC (`load_hgnc`)

Downloads `hgnc_complete_set.txt` from the HGNC Google Cloud bucket on first run. Returns a dict `{UPPERCASE_SYMBOL → "HGNC:<id>"}`.

### DOID (`load_doid`)

Downloads `doid.obo` from the Disease Ontology GitHub repo on first run. Parses all non-obsolete terms including synonyms. Returns `{lowercase_name → "DOID:XXXXXXX"}`.

### miRBase (`load_mirbase`)

Builds a three-layer lookup (highest priority first):

1. **CURRENT `id_mapping.tsv`** — downloaded from mirbase.org. Provides arm-annotated modern names (`hsa-miR-21-5p`).
2. **Previous-release `miRNA.dat` files** (releases 22, 21, 20) — parsed via `parse_mirbase_dat()`. Provides legacy names without arm suffixes (`hsa-miR-93`) and precursor names (`hsa-mir-93`).
3. **Synthesised stem entries** — for each modern `-5p` / `-3p` name in CURRENT, a no-arm stem entry is added if one does not already exist (preferring `-5p` over `-3p`).

All keys are stored **uppercase**. Values are `"MRB:<accession>"`.

---

## Curated Condition Mappings

The file `curated_condition_mappings.json` overrides automatic DOID resolution for problematic disease names. It has three top-level keys:

```json
{
  "exclude": ["physical frailty", "..."],
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
```

All keys are matched **case-insensitively** (lowercased at load time).

### Resolution order (`resolve_condition`)

1. **Exclude list** — row is dropped entirely.
2. **Exposure-agent mapping** — `condition` / `condition_id` are left empty; `exposure_agent` / `exposure_agent_id` are filled.
3. **Condition remap** — one-to-many entries produce multiple output rows (explode).
4. **Automatic DOID fallback** — exact match → fuzzy/substring match → OLS4 API query.

---

## Source Processing Functions

### `process_gene(df, hgnc, doid, curated, ncbi_api_key="")`

Handles the `gene.tsv` source. Special logic:

- **Gene fusions** are detected by `GENE_FUSION_PATTERN` (e.g. `BRCA1TP53`). Fusion rows produce **two entity rows** (one per partner gene), both linked to the same condition, with `change_type_vocab = "presence of"` and `aspect_type_vocab = "sequence variation"`.
- Rows with `Regulation == "NA"` are skipped for non-fusion genes.
- Entity IDs come from `map_gene_symbol()` → HGNC, with NCBI Gene API as a fallback.

### `process_metabolite(df, doid, curated)`

Handles the `metabolite.tsv` source.

- The `ID` column may contain comma-separated PubChem CIDs; each CID generates a separate output row.
- Entity IDs are prefixed `PCCID:`.
- `aspect_type_vocab` is always `"level"`.

### `process_mirna(df, mirbase, doid, curated)`

Handles the `mirna.tsv` source.

- Calls `_log_mirna_mapping_stats()` before filtering to report mapping coverage over the full input.
- **Mouse miRNA guard**: rows where the ID prefix is `mmu` (mouse) are dropped with a warning, as they represent data-entry errors.
- Other non-`hsa` prefixes (e.g. viral miRNAs like `ebv-`) are kept but flagged with a warning for curator review.
- Entity IDs are resolved via `map_mirna_id()` with two-tier arm-stripping fallback.
- `aspect_type_vocab` is always `"expression"`.

### `process_protein(df, doid, curated)`

Handles the `protein.tsv` source.

- Entity IDs are prefixed `UPKB:` (UniProtKB).
- Protein names are cleaned by `clean_protein_name()`, which removes all parenthetical groups.
- `aspect_type_vocab` is always `"level"`.

---

## Field Transformation Helpers

| Function | Input | Output | Notes |
|----------|-------|--------|-------|
| `clean_biomarker_name(name)` | Raw name string | Cleaned name | Removes quotes, `[bracket]` suffixes, leading/trailing underscores |
| `clean_protein_name(name)` | Raw protein name | Cleaned name | Same as above plus removes `(parenthetical)` groups |
| `clean_condition(raw)` | Raw disease string | Lowercase string | Removes quotes and parenthetical content |
| `map_specimen(raw)` | Raw collection site | Normalised specimen | Applies `SPECIMEN_MAP`; treats trailing word `saliva` as `"saliva"` |
| `map_change_type(raw)` | Raw regulation string | `"increased"` / `"decreased"` / original | Applies `CHANGE_TYPE_MAP` |
| `parse_evidence_source(pmid_val)` | Raw PMID field | `"PubMed:<id>"` or `None` | Returns `None` for non-numeric values; row is skipped by caller |
| `split_biomarker_role(raw)` | Raw biomarker type | `list[str]` | Splits on `/`; returns `[""]` for `"na"` |
| `build_biomarker_controlled_vocab(change_type, aspect_type, entity_name, entity_id)` | Four strings | Controlled-vocab string | Format: `"<change> <aspect> of <name>/<id>"` (lowercased) |

---

## ID Lookup Functions

### `map_gene_symbol(symbol, hgnc, ncbi_api_key="")`

1. Looks up the uppercased symbol in the HGNC dict.
2. On a miss, calls `_ncbi_gene_lookup()` via the NCBI ESearch API.
3. Returns `"HGNC:<id>"`, `"NCBI:<id>"`, or `""`.

### `_ncbi_gene_lookup(symbol, api_key="")`

Queries `NCBI ESearch` for `<symbol>[sym] AND "Homo sapiens"[Organism]`. Implements:
- An in-process cache (`_ncbi_gene_cache`) to avoid duplicate API calls within a run.
- Exponential back-off with up to 4 attempts on HTTP 429 responses.
- Base delay of 0.11 s (with API key) or 0.4 s (without) to stay within NCBI rate limits.

### `map_condition_doid(condition, doid)`

1. Exact lowercase match in the DOID dict.
2. Substring / fuzzy match (logs the match for review).
3. OLS4 API fallback via `_ols_lookup()`.

### `map_mirna_id(mirna_id, mirbase)`

1. Direct uppercase lookup in the merged miRBase dict.
2. Strip `-5p` / `-3p` arm suffix and retry.
3. Returns the original ID string unchanged if all tiers fail (so the row is preserved, but the entity ID will not start with `MRB:` and will be dropped by the final entity-ID filter).

---

## Logging & Diagnostics

The script uses Python's standard `logging` module with structured output configured by `setup_logging(log_dir)`. Key log messages include:

- Load counts for each mapping table.
- Per-source row counts before and after filtering.
- Every fuzzy DOID match and OLS fallback.
- Every NCBI Gene fallback lookup.
- miRNA mapping coverage summary (total / mapped / missing `-1`/`-2` suffix / other unmapped).
- Warnings for non-human miRNA prefixes and dropped mouse (`mmu`) rows.
- Final row count before and after the empty-entity-ID filter.

---

## External API Calls

| API | Used by | Rate limit | Caching |
|-----|---------|-----------|---------|
| NCBI ESearch (`eutils.ncbi.nlm.nih.gov`) | `_ncbi_gene_lookup` | 3 req/s (no key), 10 req/s (with key) | In-process dict `_ncbi_gene_cache` |
| OLS4 (`www.ebi.ac.uk/ols4/api/search`) | `_ols_lookup` | No hard limit | None (called at most once per unique unmapped condition) |
| HGNC TSV (Google Cloud) | `load_hgnc` | — | File on disk; skipped if present |
| DOID OBO (GitHub) | `load_doid` | — | File on disk; skipped if present |
| miRBase `id_mapping.tsv` (mirbase.org) | `load_mirbase` | — | File on disk; skipped if present |
| miRBase `miRNA.dat` releases 22/21/20 | `_download_mirbase_dat` | — | File on disk; skipped if present |

---

## Known Limitations & Caveats

- **`loinc_code`, `evidence`, and `tag` columns are always empty.** They are reserved for future enrichment steps.
- **`mod_type_vocab` is always `"not_specified"`** regardless of entity type.
- **Gene fusion detection is heuristic.** The regex `GENE_FUSION_PATTERN` (`^[A-Za-z]{3,}[0-9][A-Za-z]{3,}[0-9](\s+\S+)?$`) may produce false positives for unusual gene symbols.
- **OLS condition matching is not cached to disk.** Re-running the script will re-query OLS for any condition that is not in the DOID OBO file and not in the curated mappings.
- **Multi-role rows are exploded.** A row with `Biomarker_Type = "diagnostic/prognostic"` produces two output rows, one per role. The `biomarker_index` will be duplicated.
- **miRNA IDs that cannot be mapped remain as-is** in `assessed_biomarker_entity_id` and are subsequently dropped by the final empty-ID filter. Check the log for `"miRNA ID not mapped"` warnings to identify candidates for curated mappings.
- **The script requires `scripts.utils.log_utils`** to be importable. This is a local utility module not included in the script itself.
