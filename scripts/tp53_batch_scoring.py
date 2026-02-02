"""
Batch ESM-2 Scoring of TP53 Uncertain-Significance Variants
=============================================================
Reads tp53_uncertain_variants.csv, scores every variant with ESM-2,
appends the LLR column, sorts by most damaging first, and saves
tp53_ai_predictions.csv.  Also generates a histogram of scores.

Optimization: mutations at the same position share a single forward pass.
"""

import csv
import re
import sys
import time
from collections import defaultdict

import torch
import matplotlib
matplotlib.use("Agg")  # non-interactive backend (no GUI needed)
import matplotlib.pyplot as plt
from transformers import AutoTokenizer, EsmForMaskedLM

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
INPUT_CSV  = "tp53_uncertain_variants.csv"
OUTPUT_CSV = "tp53_ai_predictions.csv"
HISTOGRAM  = "tp53_llr_distribution.png"

TP53_SEQUENCE = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_protein_change(notation: str):
    """Return (position, wt_1letter, mut_1letter) or None."""
    # p.Arg175His
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', notation)
    if m:
        wt  = AA3_TO_1.get(m.group(1))
        mut = AA3_TO_1.get(m.group(3))
        if wt and mut:
            return int(m.group(2)), wt, mut

    # R175H
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', notation)
    if m:
        return int(m.group(2)), m.group(1), m.group(3)

    return None

# ---------------------------------------------------------------------------
# Device + model
# ---------------------------------------------------------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device : {device}"
      + (f" ({torch.cuda.get_device_name(0)})" if torch.cuda.is_available() else ""))

print(f"Loading {MODEL_NAME} ...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
model = EsmForMaskedLM.from_pretrained(MODEL_NAME).to(device).eval()
print("Model ready.\n")

# ---------------------------------------------------------------------------
# 1. Read CSV and parse mutations
# ---------------------------------------------------------------------------
print(f"Reading {INPUT_CSV} ...")
with open(INPUT_CSV, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
print(f"  {len(rows)} variants loaded.")

# Parse each row's protein change
parsed = []  # (index, position, wt, mut)  -- index into rows[]
skipped = 0
for i, row in enumerate(rows):
    result = parse_protein_change(row["Protein_Change"])
    if result is None:
        skipped += 1
        rows[i]["ESM2_LLR"] = ""
        rows[i]["ESM2_Interpretation"] = "parse_error"
        continue
    pos, wt, mut = result
    # Verify wild-type residue matches the sequence
    if pos < 1 or pos > len(TP53_SEQUENCE) or TP53_SEQUENCE[pos - 1] != wt:
        skipped += 1
        rows[i]["ESM2_LLR"] = ""
        rows[i]["ESM2_Interpretation"] = "wt_mismatch"
        continue
    parsed.append((i, pos, wt, mut))

if skipped:
    print(f"  Skipped {skipped} rows (parse error or WT mismatch)")

# ---------------------------------------------------------------------------
# 2. Group by position for efficient batching
# ---------------------------------------------------------------------------
pos_groups = defaultdict(list)  # position -> [(row_index, wt, mut), ...]
for i, pos, wt, mut in parsed:
    pos_groups[pos].append((i, wt, mut))

unique_positions = sorted(pos_groups.keys())
print(f"  {len(parsed)} scorable mutations across {len(unique_positions)} unique positions.\n")

# ---------------------------------------------------------------------------
# 3. Score -- one forward pass per unique position
# ---------------------------------------------------------------------------
print("Scoring with ESM-2 (one forward pass per unique position) ...")
t0 = time.time()
scored = 0

for idx, pos in enumerate(unique_positions, 1):
    # Mask the target position
    masked_seq = (
        TP53_SEQUENCE[: pos - 1]
        + tokenizer.mask_token
        + TP53_SEQUENCE[pos:]
    )
    inputs = tokenizer(masked_seq, return_tensors="pt", truncation=True)
    inputs = {k: v.to(device) for k, v in inputs.items()}

    with torch.no_grad():
        logits = model(**inputs).logits  # (1, seq_len, vocab)

    # Locate the <mask> token in the input
    mask_id = tokenizer.mask_token_id
    input_ids = inputs["input_ids"][0]
    mask_idx = (input_ids == mask_id).nonzero(as_tuple=True)[0].item()

    log_probs = torch.log_softmax(logits[0, mask_idx], dim=-1)

    # Score every mutation at this position
    for row_idx, wt, mut in pos_groups[pos]:
        wt_tid  = tokenizer.convert_tokens_to_ids(wt)
        mut_tid = tokenizer.convert_tokens_to_ids(mut)

        lp_wt  = log_probs[wt_tid].item()
        lp_mut = log_probs[mut_tid].item()
        llr    = lp_mut - lp_wt

        if llr > -0.5:
            verdict = "Likely neutral"
        elif llr > -2.0:
            verdict = "Possibly damaging"
        elif llr > -4.0:
            verdict = "Likely damaging"
        else:
            verdict = "Strongly damaging"

        rows[row_idx]["ESM2_LLR"] = f"{llr:.4f}"
        rows[row_idx]["ESM2_Interpretation"] = verdict
        scored += 1

    # Progress every 50 positions
    if idx % 50 == 0 or idx == len(unique_positions):
        elapsed = time.time() - t0
        rate = idx / elapsed
        eta = (len(unique_positions) - idx) / rate if rate > 0 else 0
        print(f"  [{idx:>4}/{len(unique_positions)}] positions done"
              f"  |  {scored} mutations scored"
              f"  |  {elapsed:.0f}s elapsed"
              f"  |  ETA {eta:.0f}s")

elapsed = time.time() - t0
print(f"\nScoring complete: {scored} mutations in {elapsed:.1f}s "
      f"({scored/elapsed:.1f} mut/s)\n")

# ---------------------------------------------------------------------------
# 4. Sort by LLR (most negative = most damaging first)
# ---------------------------------------------------------------------------
def sort_key(row):
    try:
        return float(row["ESM2_LLR"])
    except (ValueError, KeyError):
        return 0.0  # unparseable rows go to the end

rows.sort(key=sort_key)

# ---------------------------------------------------------------------------
# 5. Write output CSV
# ---------------------------------------------------------------------------
print(f"Writing {OUTPUT_CSV} ...")
fieldnames = [
    "ClinVar_ID", "Gene", "Variant_Type", "Clinical_Significance",
    "Protein_Change", "RefSeq_Protein_ID", "Title",
    "ESM2_LLR", "ESM2_Interpretation",
]
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"  {len(rows)} rows saved.\n")

# ---------------------------------------------------------------------------
# 6. Histogram
# ---------------------------------------------------------------------------
print(f"Generating histogram -> {HISTOGRAM} ...")

llr_values = []
for row in rows:
    try:
        llr_values.append(float(row["ESM2_LLR"]))
    except (ValueError, KeyError):
        pass

fig, ax = plt.subplots(figsize=(12, 6))

# Main histogram
n, bins, patches = ax.hist(
    llr_values, bins=80, edgecolor="black", linewidth=0.4, color="#4C8BF5",
    alpha=0.85,
)

# Color-code the bins by severity
for patch, left_edge in zip(patches, bins[:-1]):
    if left_edge <= -4.0:
        patch.set_facecolor("#D32F2F")   # red -- strongly damaging
    elif left_edge <= -2.0:
        patch.set_facecolor("#F57C00")   # orange -- likely damaging
    elif left_edge <= -0.5:
        patch.set_facecolor("#FBC02D")   # yellow -- possibly damaging
    else:
        patch.set_facecolor("#388E3C")   # green -- likely neutral

# Vertical threshold lines
for thresh, label, color in [
    (-4.0, "Strongly\ndamaging", "#D32F2F"),
    (-2.0, "Likely\ndamaging",   "#F57C00"),
    (-0.5, "Possibly\ndamaging", "#FBC02D"),
]:
    ax.axvline(x=thresh, color=color, linestyle="--", linewidth=1.5, alpha=0.9)
    ax.text(thresh - 0.15, ax.get_ylim()[1] * 0.92, label,
            fontsize=7.5, color=color, ha="right", va="top", fontweight="bold")

ax.set_xlabel("ESM-2 Log-Likelihood Ratio (LLR)", fontsize=12)
ax.set_ylabel("Number of Variants", fontsize=12)
ax.set_title(
    f"TP53 Uncertain-Significance Variants -- ESM-2 Mutation Scores  (n={len(llr_values)})",
    fontsize=13, fontweight="bold",
)

# Summary stats annotation
import statistics
mean_llr   = statistics.mean(llr_values)
median_llr = statistics.median(llr_values)
sd_llr     = statistics.stdev(llr_values)
strongly   = sum(1 for v in llr_values if v <= -4.0)
likely     = sum(1 for v in llr_values if -4.0 < v <= -2.0)
possibly   = sum(1 for v in llr_values if -2.0 < v <= -0.5)
neutral    = sum(1 for v in llr_values if v > -0.5)

stats_text = (
    f"Mean:   {mean_llr:.2f}\n"
    f"Median: {median_llr:.2f}\n"
    f"SD:     {sd_llr:.2f}\n"
    f"\n"
    f"Strongly damaging: {strongly}\n"
    f"Likely damaging:   {likely}\n"
    f"Possibly damaging: {possibly}\n"
    f"Likely neutral:    {neutral}"
)
ax.text(
    0.98, 0.95, stats_text, transform=ax.transAxes,
    fontsize=9, verticalalignment="top", horizontalalignment="right",
    fontfamily="monospace",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="gray", alpha=0.9),
)

ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(HISTOGRAM, dpi=200)
plt.close()
print(f"  Saved.\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("=" * 70)
print("  BATCH SCORING COMPLETE")
print("=" * 70)
print(f"  Variants scored       : {scored}")
print(f"  Output CSV            : {OUTPUT_CSV}")
print(f"  Histogram             : {HISTOGRAM}")
print(f"  Strongly damaging     : {strongly}  (LLR <= -4.0)")
print(f"  Likely damaging       : {likely}  (-4.0 < LLR <= -2.0)")
print(f"  Possibly damaging     : {possibly}  (-2.0 < LLR <= -0.5)")
print(f"  Likely neutral        : {neutral}  (LLR > -0.5)")
print("=" * 70)

# Top 10 most damaging
print(f"\n  Top 10 most damaging variants:")
print(f"  {'Protein_Change':<20} {'LLR':>8}  Interpretation")
print(f"  {'-'*20} {'-'*8}  {'-'*24}")
for row in rows[:10]:
    print(f"  {row['Protein_Change']:<20} {row['ESM2_LLR']:>8}  {row['ESM2_Interpretation']}")
