import time
import json
from pathlib import Path

import requests
import pandas as pd

# =========================================================
# CONFIG
# =========================================================

EMAIL = "your_email@example.com"   # βάλε το email σου
API_KEY = ""                       # προαιρετικό
REQUEST_DELAY = 0.12 if API_KEY else 0.40

BASE_DIR = Path("EXPERIMENT_200_GENES")
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SEARCH_DATE = "2026/03/15"

# Αν count >= 1, θεωρούμε ότι το gene εμφανίζεται στο domain
MIN_COUNT_FOR_PRESENCE = 1

# =========================================================
# DOMAINS
# =========================================================

DOMAINS = {
    "EEG": [
        '"Electroencephalography"[MeSH Terms]',
        '"Electroencephalography"[Majr]',
        '"EEG"[Title/Abstract]'
    ],
    "Cardiology": [
        '"Cardiology"[MeSH Terms]',
        '"Cardiovascular Diseases"[MeSH Terms]'
    ],
    "Endocrinology": [
        '"Endocrinology"[MeSH Terms]',
        '"Diabetes Mellitus"[MeSH Terms]',
        '"Metabolic Diseases"[MeSH Terms]'
    ],
    "Immunology": [
        '"Immunology"[MeSH Terms]',
        '"Autoimmune Diseases"[MeSH Terms]',
        '"Inflammation"[MeSH Terms]'
    ],
    "Oncology": [
        '"Medical Oncology"[MeSH Terms]',
        '"Neoplasms"[MeSH Terms]'
    ]
}

DOMAIN_NAMES = ["EEG", "Cardiology", "Endocrinology", "Immunology", "Oncology"]

# =========================================================
# LOAD GENES
# =========================================================

GENE_FILE = Path("genes.txt")
if not GENE_FILE.exists():
    raise FileNotFoundError("Δεν βρέθηκε genes.txt στον ίδιο φάκελο με το script.")

with open(GENE_FILE, "r", encoding="utf-8") as f:
    GENES = [line.strip() for line in f if line.strip()]

if not GENES:
    raise ValueError("Το genes.txt είναι άδειο.")

print(f"\nGenes loaded: {len(GENES)}\n")

# =========================================================
# HELPERS
# =========================================================

def ncbi_get_json(endpoint: str, params: dict) -> dict:
    all_params = {
        "tool": "experiment_200_genes_overlap",
        "email": EMAIL,
        "retmode": "json"
    }
    if API_KEY:
        all_params["api_key"] = API_KEY
    all_params.update(params)

    last_error = None
    for attempt in range(5):
        try:
            time.sleep(REQUEST_DELAY)
            r = requests.get(NCBI_EUTILS + endpoint, params=all_params, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_error = e
            print(f"NCBI error, retry {attempt+1}/5 ...")
            time.sleep(2)

    raise RuntimeError(f"NCBI API failed after retries: {last_error}")


def build_or_block(terms):
    return "(" + " OR ".join(terms) + ")"


def build_query(gene: str, domain: str) -> str:
    gene_block = f'("{gene}"[Title/Abstract])'
    domain_block = build_or_block(DOMAINS[domain])

    publication_filter = (
        '("journal article"[Publication Type] OR "review"[Publication Type]) '
        'NOT ("editorial"[Publication Type] OR "letter"[Publication Type])'
    )

    date_filter = f'("1800/01/01"[Date - Publication] : "{SEARCH_DATE}"[Date - Publication])'

    return f"({gene_block} AND {domain_block}) AND ({publication_filter}) AND ({date_filter})"


def esearch_count(query: str) -> int:
    data = ncbi_get_json("esearch.fcgi", {
        "db": "pubmed",
        "term": query,
        "retmax": 0
    })
    return int(data["esearchresult"]["count"])


# =========================================================
# MAIN
# =========================================================

rows = []
query_log_rows = []

print("Collecting PubMed counts...\n")

for i, gene in enumerate(GENES, start=1):
    print(f"[{i}/{len(GENES)}] {gene}")

    row = {"Gene": gene}

    for domain in DOMAIN_NAMES:
        q = build_query(gene, domain)
        c = esearch_count(q)

        row[f"{domain}_count"] = c
        row[f"{domain}_present"] = int(c >= MIN_COUNT_FOR_PRESENCE)

        query_log_rows.append({
            "Gene": gene,
            "Domain": domain,
            "Query": q,
            "Count": c
        })

    row["Domains_present"] = int(sum(row[f"{d}_present"] for d in DOMAIN_NAMES))
    row["NonEEG_domains_present"] = int(sum(row[f"{d}_present"] for d in DOMAIN_NAMES if d != "EEG"))

    rows.append(row)

df = pd.DataFrame(rows)
query_df = pd.DataFrame(query_log_rows)

# =========================================================
# OVERLAPS
# =========================================================

def gene_set(domain):
    return set(df.loc[df[f"{domain}_present"] == 1, "Gene"])

eeg_set = gene_set("EEG")
card_set = gene_set("Cardiology")
endo_set = gene_set("Endocrinology")
immun_set = gene_set("Immunology")
onco_set = gene_set("Oncology")

overlap_eeg_card = sorted(eeg_set & card_set)
overlap_eeg_endo = sorted(eeg_set & endo_set)
overlap_eeg_immun = sorted(eeg_set & immun_set)
overlap_eeg_onco = sorted(eeg_set & onco_set)
overlap_eeg_all4 = sorted(eeg_set & card_set & endo_set & immun_set & onco_set)

summary_df = pd.DataFrame([
    {"Comparison": "EEG genes", "Count": len(eeg_set)},
    {"Comparison": "Cardiology genes", "Count": len(card_set)},
    {"Comparison": "Endocrinology genes", "Count": len(endo_set)},
    {"Comparison": "Immunology genes", "Count": len(immun_set)},
    {"Comparison": "Oncology genes", "Count": len(onco_set)},
    {"Comparison": "EEG ∩ Cardiology", "Count": len(overlap_eeg_card)},
    {"Comparison": "EEG ∩ Endocrinology", "Count": len(overlap_eeg_endo)},
    {"Comparison": "EEG ∩ Immunology", "Count": len(overlap_eeg_immun)},
    {"Comparison": "EEG ∩ Oncology", "Count": len(overlap_eeg_onco)},
    {"Comparison": "EEG ∩ all 4 other domains", "Count": len(overlap_eeg_all4)},
])

# =========================================================
# SAVE
# =========================================================

df = df.sort_values(
    ["Domains_present", "NonEEG_domains_present", "EEG_count"],
    ascending=[False, False, False]
).reset_index(drop=True)

df.to_csv(RESULTS_DIR / "gene_domain_table.csv", index=False)
query_df.to_csv(RESULTS_DIR / "query_log.csv", index=False)
summary_df.to_csv(RESULTS_DIR / "overlap_summary.csv", index=False)

pd.DataFrame({"Gene": overlap_eeg_card}).to_csv(RESULTS_DIR / "overlap_EEG_Cardiology.csv", index=False)
pd.DataFrame({"Gene": overlap_eeg_endo}).to_csv(RESULTS_DIR / "overlap_EEG_Endocrinology.csv", index=False)
pd.DataFrame({"Gene": overlap_eeg_immun}).to_csv(RESULTS_DIR / "overlap_EEG_Immunology.csv", index=False)
pd.DataFrame({"Gene": overlap_eeg_onco}).to_csv(RESULTS_DIR / "overlap_EEG_Oncology.csv", index=False)
pd.DataFrame({"Gene": overlap_eeg_all4}).to_csv(RESULTS_DIR / "overlap_EEG_all4domains.csv", index=False)

with open(RESULTS_DIR / "overlap_report.txt", "w", encoding="utf-8") as f:
    f.write(summary_df.to_string(index=False))

with open(RESULTS_DIR / "config.json", "w", encoding="utf-8") as f:
    json.dump({
        "email": EMAIL,
        "search_date": SEARCH_DATE,
        "presence_threshold": MIN_COUNT_FOR_PRESENCE,
        "domains": DOMAINS,
        "n_genes": len(df)
    }, f, indent=2, ensure_ascii=False)

# =========================================================
# OUTPUT
# =========================================================

print("\n========================================")
print("EXPERIMENT COMPLETED")
print("========================================")
print(summary_df.to_string(index=False))
print("\nSaved files:")
print("-", RESULTS_DIR / "gene_domain_table.csv")
print("-", RESULTS_DIR / "overlap_summary.csv")
print("-", RESULTS_DIR / "overlap_report.txt")
print("-", RESULTS_DIR / "query_log.csv")
print("========================================")