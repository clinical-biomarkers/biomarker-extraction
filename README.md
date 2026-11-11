# biomarker-extraction
Raw data extraction and transformation scripts

---

## caDSR

### Downloads

#### Permissible Values via Browser Export (Recommended)

The caDSR OneData UI at `cadsr.cancer.gov/onedata` offers an "Export to Excel" option under each data element's Permissible Values tab. The downloaded file (`od_001.xls`) contains the following columns:

| Column | Description |
|---|---|
| Permissible Value | The actual value |
| VM Long Name | Value Meaning long name |
| VM Public ID | Value Meaning public identifier |
| Concept Codes | Associated NCI Thesaurus concept codes |
| VM Definition | Value Meaning definition |
| Begin Date | Date the value became active |
| End Date | Date the value was retired (if applicable) |

The export requires an active browser session and cannot be replicated with `wget` or `curl` directly, as the download URLs contain short-lived session tokens.

1. Go to: https://cadsr.cancer.gov/onedata/dmdirect/NIH/NCI/CO/CDEDD?filter=CDEDD.ITEM_ID=5473%20and%20ver_nr=14
2. Click the notepad icon under "View" (2nd column)

   <img width="488" height="214" alt="1" src="https://github.com/user-attachments/assets/72a9bf90-0a1c-4c98-8e09-d99dc6c9b691" />
4. Select "6. Permissible Values"
   <img width="1879" height="807" alt="2" src="https://github.com/user-attachments/assets/43670058-fedc-4f8e-9799-dc7dcca1fadc" />
5. For bulk download go to "Delivery Options" and select "Export to Excel"
   <img width="1960" height="912" alt="3" src="https://github.com/user-attachments/assets/1501a860-001d-4fe7-93e8-95cb65227629" />

The downloaded sheet can be found in `downloads/cadsr/od_001.xls` on the BiomarkerKB dev server. The downloads folder is not git tracked.

#### Permissible Values via caDSR API (TBA)

More research is needed to determine the correct API endpoint.

### Generated files
Location: `generated/cadsr`

`cadsr.tsv` is a manually curated selection of rows from `od_001.xls`.

---

## CarboCurator
### Download
Cyrus's data: 

https://github.com/clinical-biomarkers/carbocurator-data/

https://github.com/clinical-biomarkers/bGSL-data/

The data files live in these repositories and present the point of truth for this data:

https://github.com/clinical-biomarkers/carbocurator-data/blob/main/releases/2026-10-08/biomarkers_llm_glycan_2026-10-08.tsv 

https://github.com/clinical-biomarkers/bGSL-data/blob/main/releases/2026-10-08/bGSL_2026-10-08.tsv 

---

## SalivaDB
### Download
```bash
cd scripts/salivadb
./download_salivadb.sh -d path/to/download/dir
```
### generate_salivadb.py

Generates the SalivaDB biomarkers TSV dataset from the four raw download files
(`gene.tsv`, `metabolite.tsv`, `mirna.tsv`, `protein.tsv`).

#### Requirements

- Python 3.9+
- pip packages:

```bash
pip install pandas requests
```

#### Usage

```bash
python generate_salivadb.py
```

or with explicit paths:

```bash
python generate_salivadb.py \
  --input-dir  /data/shared/repos/extraction/downloads/salivadb \
  --output-file /data/shared/repos/extraction/generated/salivadb/biomarkers_salivadb.tsv \
  --mapping-dir /data/shared/repos/extraction/downloads/mapping_files
```

#### Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--input-dir` | No | `/data/shared/repos/extraction/downloads/salivadb` | Directory containing the four source TSVs |
| `--output-file` | No | `/data/shared/repos/extraction/generated/salivadb/biomarkers_salivadb.tsv` | Path for the output TSV |
| `--mapping-dir` | No | `/data/shared/repos/extraction/downloads/mapping_files` | Directory where mapping files are stored (downloaded on first run) |

#### Output

A single TSV at the path specified by `--output-file` with these columns:

```
biomarker_index, component_index, entity_index, assessed_biomarker_entity,
assessed_biomarker_entity_id, assessed_entity_type, condition, condition_id,
best_biomarker_role, specimen, specimen_id, loinc_code, evidence_source,
evidence, biomarker, exposure_agent, exposure_agent_id, tag,
biomarker_controlled_vocab, change_type_vocab, aspect_type_vocab,
mod_type_vocab, legacy_biomarker_id
```

#### Logging

Timestamped log files are written to:

```
/data/shared/repos/extraction/logs/generate_salivadb_<YYYYMMDD_HHMMSS>.log
```

Each run creates a new log file. The log captures:
- Download decisions (skipped vs. fetched) for each mapping file
- Fuzzy DOID matches (name used vs. original query)
- OLS fallback attempts and results
- Unmapped gene symbols, miRNA IDs, and conditions flagged for human review
- PubMed resolution of author-string citations (success and failure)
- Row counts at each pipeline stage

#### Mapping files downloaded automatically

| File | Source | Cached at |
|---|---|---|
| `hgnc_complete_set.txt` | HGNC (Google Storage) | `mapping_files/hgnc_complete_set.txt` |
| `gene_info.gz` | NCBI FTP | `mapping_files/gene_info.gz` |
| `doid.obo` | Disease Ontology GitHub | `mapping_files/doid.obo` |
| `mirbase_id_mapping.tsv` | miRBase | `mapping_files/mirbase_id_mapping.tsv` |

Re-running the script will **not** re-download these files if they already exist.
Delete the file manually to force a fresh download.

#### Processing pipeline overview

```
gene.tsv        ─┐
metabolite.tsv   ├─► filter (Homo sapiens, exclude Healthy/NA conditions)
mirna.tsv        │   → field transformations
protein.tsv     ─┘   → ID mapping (HGNC/NCBI/UPKB/PCCID/MRB)
                     → condition → DOID (exact → fuzzy → OLS → human review)
                     → evidence_source: numeric PMID or PubMed API resolution
                     → role explosion (diagnostic/prognostic → 2 rows)
                     → biomarker_controlled_vocab construction
                     → concat → write TSV
```

##### Source-specific notes

**gene.tsv**
- Gene symbols are mapped to HGNC IDs first, then NCBI as fallback.
- Gene fusions (IDs matching `^[A-Za-z]{3,}[0-9][A-Za-z]{3,}[0-9](\s+\S+)?$`) produce
  two entity rows (entity_index 1 and 2), `change_type_vocab = "presence of"`,
  `aspect_type_vocab = "sequence variation"`.
- Non-fusion rows where `Regulation == "NA"` are excluded.

**metabolite.tsv**
- Multiple IDs in the `ID` column (comma-separated) produce multiple rows.
- IDs are prefixed with `PCCID:`.

**mirna.tsv**
- IDs are mapped via miRBase `id_mapping.tsv` to `MRB:<accession>`.
- Falls back to stem-loop accession (stripping -5p/-3p suffix) when the
  mature ID is not found directly.

**protein.tsv**
- UniProtKB IDs prefixed with `UPKB:`.
- Parenthetical content stripped from `Biomarker_Name`.

#### Condition ID mapping (three-pass)

1. **Exact / substring match** against DOID OBO (names + synonyms). Fuzzy matches
   are logged at INFO level.
2. **OLS4 API fallback** (`https://www.ebi.ac.uk/ols4/api/search`).
3. **Human review** — rows where neither pass found a DOID emit a WARNING in the log
   with the condition string. Search the log for `needs human review` to find them.

#### Evidence source resolution

- Numeric `PMID` values → `PubMed:<id>` directly.
- Author string values → queried against NCBI ESearch using the author string
  and `Year_of_Publication`. On success → `PubMed:<id>`; on failure → the raw
  author string is kept and a WARNING is logged.

#### Notes / known limitations

1. **Internet access required on first run** to download mapping files and to
   call the OLS4 and PubMed APIs at run time.
2. **Gene fusion splitting** uses a heuristic regex; unusual fusion notations
   may not split correctly. Review WARNING log lines for `Gene symbol not mapped`.
3. **DOID OBO fuzzy matching** is a simple substring check; it may occasionally
   map to a broader term. Review INFO lines for `Fuzzy DOID match`.
4. **`condition_id` for human review** — search the log for `needs human review`
   after each run and populate those values manually or via a curated override file.
5. The NCBI `gene_info.gz` file is several hundred MB; the first download may take
   a few minutes.
---
