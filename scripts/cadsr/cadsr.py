import pandas as pd
import requests
import sys
import time
from pathlib import Path

ROOT_DIR = Path(__file__).parents[2]
DL_DIR = ROOT_DIR / "downloads"
#OUT_DIR = '/data/shared/biomarkerdb/downloads/cadsr/current'

INPUT_FILE = DL_DIR / "cadsr" / "od_001.xls"
OUTPUT_FILE = "od_001_omim.tsv"
API_BASE = "https://api-evsrest.nci.nih.gov/api/v1/concept/ncit"

def has_omim(concept_code):
    if not concept_code or pd.isna(concept_code):
        return False
    code = str(concept_code).strip()
    try:
        r = requests.get(f"{API_BASE}/{code}?limit=100&include=summary", timeout=15)
        if r.status_code != 200:
            return False
        data = r.json()
        return any(p.get("type") == "OMIM_Number" for p in data.get("properties", []))
    except Exception as e:
        print(f"Error fetching {code}: {e}", file=sys.stderr)
        return False

df = pd.read_csv(INPUT_FILE.with_suffix(".csv"), dtype=str, encoding="latin-1")

print(df.columns.tolist())
print(df["Concept Codes"].head(10).tolist())
print(has_omim("C104175"))  # should print True

results = []
for i, row in df.iterrows():
    code = row.get("Concept Codes")
    print(f"Checking row {i}: {code}", file=sys.stderr)
    if has_omim(code):
        results.append(row)
    time.sleep(0.1)

out = pd.DataFrame(results)
out.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"Done. {len(results)} rows written to {OUTPUT_FILE}")
