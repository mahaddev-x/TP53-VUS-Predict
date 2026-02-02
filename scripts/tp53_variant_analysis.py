"""
Variant Effect Prediction Analysis - TP53 Gene
================================================
Searches ClinVar for TP53 variants of uncertain significance,
filters for SNV / Missense types, and exports protein-change data to CSV.
Also fetches the wild-type TP53 protein sequence from UniProt.
"""

import csv
import time
import re
import requests
from Bio import Entrez

# -- Configuration ------------------------------------------------------------
Entrez.email = "variant_analysis@example.com"
GENE = "TP53"
OUTPUT_CSV = "tp53_uncertain_variants.csv"
OUTPUT_SEQ = "tp53_wildtype_sequence.txt"
UNIPROT_ID = "P04637"
BATCH_SIZE = 500

# Known RefSeq protein for TP53 canonical transcript NM_000546
TP53_REFSEQ_PROTEIN = "NP_000537.3"

# -- 1. Search ClinVar (usehistory keeps results server-side) -----------------
print(f"[1/5] Searching ClinVar for {GENE} uncertain-significance variants ...")

query = f'{GENE}[gene] AND "uncertain significance"[clinical_significance]'

handle = Entrez.esearch(db="clinvar", term=query, retmax=0, usehistory="y")
record = Entrez.read(handle, validate=False)
handle.close()

total = int(record["Count"])
webenv = record["WebEnv"]
query_key = record["QueryKey"]
print(f"       Found {total} records (stored on NCBI history server)")

# -- 2. Fetch summaries via history server ------------------------------------
print(f"[2/5] Fetching summaries in batches of {BATCH_SIZE} ...")

all_summaries = []
for start in range(0, total, BATCH_SIZE):
    for attempt in range(1, 5):
        try:
            handle = Entrez.esummary(
                db="clinvar",
                query_key=query_key,
                WebEnv=webenv,
                retstart=start,
                retmax=BATCH_SIZE,
            )
            summaries = Entrez.read(handle, validate=False)
            handle.close()
            break
        except Exception as e:
            wait = 2 ** attempt
            print(f"       [retry {attempt}] {e} -- waiting {wait}s")
            time.sleep(wait)
    else:
        raise RuntimeError("Max retries exceeded fetching summaries")

    doc_sums = summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    all_summaries.extend(doc_sums)
    fetched = min(start + BATCH_SIZE, total)
    print(f"       ... {fetched}/{total}")
    time.sleep(0.4)

print(f"       Total summaries: {len(all_summaries)}")

# -- 3. Filter for SNV / Missense and extract protein changes -----------------
print("[3/5] Filtering for SNV / Missense with protein changes ...")

# Three-letter amino acid codes (excludes Ter = stop codon)
AA3 = (
    "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|"
    "Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
)

# Matches true missense: e.g. p.Arg175His (AA -> AA, not Ter/stop)
missense_3letter_re = re.compile(
    rf'p\.({AA3})(\d+)({AA3})'
)

# One-letter: e.g. R175H (single uppercase letter + digits + single uppercase letter)
# Excludes * and X which denote stop codons
VALID_AA1 = set("ACDEFGHIKLMNPQRSTVWY")

TYPE_KW = {"single nucleotide variant", "missense", "missense variant"}

rows = []
for doc in all_summaries:
    # ---- variant type check ----
    var_type = str(doc.get("obj_type", "")).strip().lower()
    title = str(doc.get("title", ""))
    title_lower = title.lower()

    variation_set = doc.get("variation_set", [])
    mol_con = ""
    for vs in (variation_set if isinstance(variation_set, list) else []):
        if isinstance(vs, dict):
            mol_con += " " + str(vs.get("molecular_consequence", ""))

    combined = f"{var_type} {title_lower} {mol_con}".lower()
    if not any(kw in combined for kw in TYPE_KW):
        continue

    # ---- extract protein change from title (three-letter form) ----
    m = missense_3letter_re.search(title)
    if m:
        protein_change = f"p.{m.group(1)}{m.group(2)}{m.group(3)}"
    else:
        # Fall back to one-letter form from protein_change field
        pc_field = str(doc.get("protein_change", ""))
        # e.g. "R324W, R231W, R363W, R204W" -- pick last (canonical isoform)
        one_letter = None
        for token in reversed(pc_field.replace(",", " ").split()):
            token = token.strip()
            match = re.match(r'^([A-Z])(\d+)([A-Z])$', token)
            if match and match.group(1) in VALID_AA1 and match.group(3) in VALID_AA1:
                one_letter = token
                break
        if one_letter:
            protein_change = one_letter
        else:
            continue  # no valid missense change found

    # ---- RefSeq protein ID ----
    # If the title references NM_000546 (TP53 canonical transcript), use known NP
    refseq_protein = ""
    if "NM_000546" in title:
        refseq_protein = TP53_REFSEQ_PROTEIN

    uid = str(doc.attributes.get("uid", "")) if hasattr(doc, "attributes") else ""

    rows.append({
        "ClinVar_ID": uid,
        "Gene": GENE,
        "Variant_Type": var_type.title() if var_type else "Single Nucleotide Variant",
        "Clinical_Significance": "Uncertain significance",
        "Protein_Change": protein_change,
        "RefSeq_Protein_ID": refseq_protein,
        "Title": title,
    })

print(f"       Matched: {len(rows)}")

# Deduplicate by Protein_Change
seen = set()
unique = []
for r in rows:
    if r["Protein_Change"] not in seen:
        seen.add(r["Protein_Change"])
        unique.append(r)
rows = unique
print(f"       Unique protein changes: {len(rows)}")

# -- 4. Export CSV -------------------------------------------------------------
print(f"[4/5] Writing {OUTPUT_CSV} ...")

fieldnames = [
    "ClinVar_ID", "Gene", "Variant_Type",
    "Clinical_Significance", "Protein_Change",
    "RefSeq_Protein_ID", "Title",
]
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"       {len(rows)} rows saved.")

# -- 5. Fetch wild-type TP53 from UniProt -------------------------------------
print(f"[5/5] Fetching wild-type TP53 from UniProt ({UNIPROT_ID}) ...")

url = f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ID}.fasta"
resp = requests.get(url, timeout=30)
resp.raise_for_status()

fasta_lines = resp.text.strip().split("\n")
header = fasta_lines[0]
sequence = "".join(fasta_lines[1:])

with open(OUTPUT_SEQ, "w", encoding="utf-8") as f:
    f.write(f"UniProt ID : {UNIPROT_ID}\n")
    f.write(f"Gene       : {GENE}\n")
    f.write(f"Organism   : Homo sapiens (Human)\n")
    f.write(f"Header     : {header}\n")
    f.write(f"Length     : {len(sequence)} amino acids\n")
    f.write(f"\n{'='*60}\n")
    f.write("Protein Sequence (FASTA):\n")
    f.write(f"{'='*60}\n\n")
    f.write(f"{header}\n")
    for i in range(0, len(sequence), 60):
        f.write(sequence[i : i + 60] + "\n")

print(f"       Saved ({len(sequence)} aa).")

# -- Summary -------------------------------------------------------------------
print("\n" + "=" * 60)
print("  ANALYSIS COMPLETE")
print("=" * 60)
print(f"  Gene                  : {GENE}")
print(f"  ClinVar hits (total)  : {total}")
print(f"  SNV/Missense filtered : {len(rows)}")
print(f"  CSV output            : {OUTPUT_CSV}")
print(f"  Wild-type sequence    : {OUTPUT_SEQ}")
print("=" * 60)

if rows:
    print(f"\n  Preview (first 10 of {len(rows)}):")
    print(f"  {'Protein_Change':<20} {'RefSeq_Protein_ID':<20} {'ClinVar_ID'}")
    print(f"  {'-'*20} {'-'*20} {'-'*12}")
    for r in rows[:10]:
        print(f"  {r['Protein_Change']:<20} {r['RefSeq_Protein_ID']:<20} {r['ClinVar_ID']}")
    if len(rows) > 10:
        print(f"  ... and {len(rows) - 10} more in the CSV")
