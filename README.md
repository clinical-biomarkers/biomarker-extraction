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

See documentation here: https://github.com/clinical-biomarkers/biomarker-extraction/blob/main/scripts/salivadb/README.md
---
