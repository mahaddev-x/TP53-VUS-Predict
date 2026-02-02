"""
ESM-2 Zero-Shot Mutation Scoring for TP53
==========================================
Uses Meta's esm2_t33_650M_UR50D protein language model to compute
Log-Likelihood Ratios (LLR) for single amino-acid substitutions.

    LLR = log P(mutant | context) - log P(wildtype | context)

  - LLR ~ 0    : mutation is roughly as probable as wild-type (neutral)
  - LLR << 0   : model considers the mutation unnatural (likely damaging)
  - LLR > 0    : mutation is *more* probable than wild-type (rare, possible gain)

Runs on GPU when available, falls back to CPU otherwise.
"""

import torch
import math
from transformers import AutoTokenizer, EsmForMaskedLM

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MODEL_NAME = "facebook/esm2_t33_650M_UR50D"

# Wild-type TP53 protein sequence (UniProt P04637, 393 aa)
TP53_SEQUENCE = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

# Standard amino acid alphabet
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Three-letter to one-letter mapping
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# ---------------------------------------------------------------------------
# Device selection
# ---------------------------------------------------------------------------
if torch.cuda.is_available():
    DEVICE = torch.device("cuda")
    print(f"[GPU]  Using {torch.cuda.get_device_name(0)}")
    print(f"       VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
else:
    DEVICE = torch.device("cpu")
    print("[CPU]  CUDA not available -- running on CPU (slower)")

# ---------------------------------------------------------------------------
# Load model and tokenizer
# ---------------------------------------------------------------------------
print(f"\nLoading ESM-2 model: {MODEL_NAME} ...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
model = EsmForMaskedLM.from_pretrained(MODEL_NAME)
model = model.to(DEVICE)
model.eval()
print("Model loaded.\n")


def get_mutation_score(
    sequence: str,
    position: int,
    mutant_residue: str,
) -> dict:
    """Score a single amino-acid substitution using ESM-2 masked marginals.

    Parameters
    ----------
    sequence : str
        Wild-type protein sequence (one-letter codes).
    position : int
        1-based residue position to mutate.
    mutant_residue : str
        One-letter code of the substituted amino acid.

    Returns
    -------
    dict with keys:
        wt_residue        - original amino acid at *position*
        mutant_residue    - the substitution
        position          - the queried position (1-based)
        log_prob_wt       - log probability of the wild-type residue
        log_prob_mut      - log probability of the mutant residue
        llr               - log-likelihood ratio (mut - wt)
        interpretation    - human-readable verdict
    """
    # --- validate inputs ---
    if position < 1 or position > len(sequence):
        raise ValueError(
            f"Position {position} out of range [1, {len(sequence)}]"
        )
    mutant_residue = mutant_residue.upper()
    if mutant_residue not in AMINO_ACIDS:
        raise ValueError(f"Invalid amino acid: {mutant_residue}")

    wt_residue = sequence[position - 1]  # convert to 0-based
    if mutant_residue == wt_residue:
        return {
            "wt_residue": wt_residue,
            "mutant_residue": mutant_residue,
            "position": position,
            "log_prob_wt": 0.0,
            "log_prob_mut": 0.0,
            "llr": 0.0,
            "interpretation": "Synonymous (no amino-acid change)",
        }

    # --- mask the target position and tokenize ---
    masked_seq = sequence[: position - 1] + tokenizer.mask_token + sequence[position:]

    inputs = tokenizer(masked_seq, return_tensors="pt", truncation=True)
    inputs = {k: v.to(DEVICE) for k, v in inputs.items()}

    # --- forward pass (no gradient needed) ---
    with torch.no_grad():
        logits = model(**inputs).logits  # (1, seq_len, vocab_size)

    # The tokenizer prepends <cls> so the mask is at token index = position
    mask_token_id = tokenizer.mask_token_id
    input_ids = inputs["input_ids"][0]
    mask_idx = (input_ids == mask_token_id).nonzero(as_tuple=True)[0].item()

    # --- extract log-probabilities over the vocabulary at the masked position ---
    log_probs = torch.log_softmax(logits[0, mask_idx], dim=-1)

    wt_token_id = tokenizer.convert_tokens_to_ids(wt_residue)
    mut_token_id = tokenizer.convert_tokens_to_ids(mutant_residue)

    log_prob_wt = log_probs[wt_token_id].item()
    log_prob_mut = log_probs[mut_token_id].item()
    llr = log_prob_mut - log_prob_wt

    # --- interpret ---
    if llr > -0.5:
        verdict = "Likely neutral / tolerated"
    elif llr > -2.0:
        verdict = "Possibly damaging"
    elif llr > -4.0:
        verdict = "Likely damaging"
    else:
        verdict = "Strongly predicted damaging"

    return {
        "wt_residue": wt_residue,
        "mutant_residue": mutant_residue,
        "position": position,
        "log_prob_wt": round(log_prob_wt, 4),
        "log_prob_mut": round(log_prob_mut, 4),
        "llr": round(llr, 4),
        "interpretation": verdict,
    }


def parse_hgvs(notation: str) -> tuple:
    """Parse 'p.Arg175His' or 'R175H' into (position, mutant_one_letter).

    Returns (position: int, mutant: str) or raises ValueError.
    """
    import re

    # Three-letter form: p.Arg175His
    m = re.match(
        r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', notation
    )
    if m:
        pos = int(m.group(2))
        mut = AA3_TO_1.get(m.group(3))
        if mut is None:
            raise ValueError(f"Unknown amino acid: {m.group(3)}")
        return pos, mut

    # One-letter form: R175H
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', notation)
    if m:
        return int(m.group(2)), m.group(3)

    raise ValueError(f"Cannot parse mutation notation: {notation}")


# ---------------------------------------------------------------------------
# Demo: score a few well-known TP53 mutations
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    demo_mutations = [
        "p.Arg175His",   # R175H -- classic hotspot, structural
        "p.Arg248Trp",   # R248W -- DNA-contact hotspot
        "p.Arg273His",   # R273H -- DNA-contact hotspot
        "p.Gly245Ser",   # G245S -- structural hotspot
        "p.Arg282Trp",   # R282W -- structural hotspot
        "p.Ala74Thr",    # from our uncertain-significance CSV
        "p.Pro77Leu",    # from our uncertain-significance CSV
        "p.Arg363Trp",   # from our uncertain-significance CSV
    ]

    print("=" * 78)
    print("  ESM-2 Zero-Shot Mutation Scoring -- TP53 (P04637)")
    print("=" * 78)
    print(f"  Sequence length : {len(TP53_SEQUENCE)} aa")
    print(f"  Model           : {MODEL_NAME}")
    print(f"  Device          : {DEVICE}")
    print("=" * 78)

    print(f"\n  {'Mutation':<16} {'WT':>3} {'Pos':>5} {'MUT':>3}"
          f"  {'log P(wt)':>10} {'log P(mut)':>10} {'LLR':>8}  Interpretation")
    print(f"  {'-'*16} {'-'*3} {'-'*5} {'-'*3}"
          f"  {'-'*10} {'-'*10} {'-'*8}  {'-'*28}")

    for notation in demo_mutations:
        pos, mut = parse_hgvs(notation)
        result = get_mutation_score(TP53_SEQUENCE, pos, mut)

        print(
            f"  {notation:<16} {result['wt_residue']:>3} {result['position']:>5} "
            f"{result['mutant_residue']:>3}  {result['log_prob_wt']:>10.4f} "
            f"{result['log_prob_mut']:>10.4f} {result['llr']:>8.4f}  "
            f"{result['interpretation']}"
        )

    print("\n" + "=" * 78)
    print("  LLR interpretation guide:")
    print("    LLR >  -0.5  : Likely neutral / tolerated")
    print("    LLR >  -2.0  : Possibly damaging")
    print("    LLR >  -4.0  : Likely damaging")
    print("    LLR <= -4.0  : Strongly predicted damaging")
    print("=" * 78)
