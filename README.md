# TP53-VUS-Predict

**Cross-Validating Variants of Uncertain Significance via ESM-2 & AlphaMissense**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)

---

## Abstract

The tumor suppressor TP53 harbors over 1,200 missense variants classified as
"Uncertain Significance" (VUS) in ClinVar, leaving clinicians without
actionable guidance for carrier patients. This project applies two independent
AI-based pathogenicity predictors -- Meta's **ESM-2** protein language model
(zero-shot log-likelihood ratios) and Google DeepMind's **AlphaMissense**
(proteome-wide deep learning scores) -- to systematically reclassify these
variants. We retrieved 1,211 TP53 VUS from ClinVar via the NCBI Entrez API,
scored each with ESM-2 (model `esm2_t33_650M_UR50D`), and cross-referenced
against 7,467 AlphaMissense predictions for UniProt P04637. Of the 1,199
variants successfully matched, both models showed strong concordance (Pearson
*r* = -0.706, Spearman *rho* = -0.714), with 349 variants flagged as
high-risk by both methods simultaneously. Structural analysis on PDB 1TUP
confirmed that the top five candidates occupy functionally critical positions:
two at direct DNA-contact residues (R248, R280), one at a zinc-coordination
site (C176), and two within the hydrophobic core of the DNA-binding domain
(L257, V157). These results provide computational evidence supporting
reclassification of five ClinVar VUS as likely pathogenic, with potential
implications for cancer risk assessment and clinical reporting.

---

## Methodology

```
ClinVar (NCBI Entrez API)
       |
       v
[1] tp53_variant_analysis.py
    Query TP53 + "uncertain significance" + missense SNVs
    Output: 1,211 VUS variants
       |
       v
[2] tp53_batch_scoring.py
    Load ESM-2 (esm2_t33_650M_UR50D, 650M params)
    Zero-shot LLR scoring per variant
    Output: tp53_ai_predictions.csv (1,211 scored variants)
       |
       v
[3] tp53_domain_analysis.py
    Map variants to TP53 functional domains:
      - Transactivation Domain (1-92)
      - DNA-Binding Domain (102-292)
      - Tetramerization Domain (325-356)
    Statistical tests (Welch's t, Cohen's d)
    Output: domain comparison plots
       |
       v
[4] alphamissense_crossref.py                    AlphaMissense
    Download proteome predictions (1.12 GB)  <--- (Google Cloud)
    Stream-filter for TP53 (UniProt P04637)
    Merge on protein_variant key
    Output: tp53_esm2_alphamissense_merged.csv (1,199 matched)
       |
       v
[5] pymol_variant_views.py
    Structural analysis on PDB 1TUP (p53-DNA complex)
    BioPython NeighborSearch: DNA contacts, Zn2+ sites
    Generate PyMOL .pml script for ray-traced renders
    Output: 6 publication-quality PNGs per method
```

### Scoring Methods

| Method | Model | Metric | Threshold |
|--------|-------|--------|-----------|
| **ESM-2** | `esm2_t33_650M_UR50D` (650M params) | Log-Likelihood Ratio (LLR) | <= -4.0 = Strongly damaging |
| **AlphaMissense** | DeepMind proteome-wide predictor | Pathogenicity score (0-1) | > 0.564 = Pathogenic |

### Concordance

- **1,199** variants matched between ESM-2 and AlphaMissense
- **349** flagged high-risk by both models (ESM-2 LLR <= -4.0 AND AM > 0.564)
- **438** flagged low-risk by both models
- Pearson *r* = **-0.706** | Spearman *rho* = **-0.714**

---

## Key Findings: Top 5 Confirmed Candidates

All five variants are classified as **"Uncertain Significance"** in ClinVar
but scored as **strongly damaging** by both AI models and confirmed by
structural analysis on PDB 1TUP.

| Rank | Variant | ClinVar ID | ESM-2 LLR | AlphaMissense | Structural Role | Domain |
|:----:|---------|:----------:|:----------:|:-------------:|-----------------|--------|
| 1 | **p.Leu257Arg** (L257R) | 142134 | **-13.88** | **0.9938** (pathogenic) | Buried hydrophobic core | DNA-Binding |
| 2 | **p.Val157Asp** (V157D) | 482231 | **-13.69** | **0.9992** (pathogenic) | Buried hydrophobic core | DNA-Binding |
| 3 | **p.Arg248Pro** (R248P) | 237954 | **-12.63** | **0.9994** (pathogenic) | Direct DNA contact (3.4 A) | DNA-Binding |
| 4 | **p.Cys176Arg** (C176R) | 376573 | **-12.53** | **0.9999** (pathogenic) | Zinc coordination site (2.3 A) | DNA-Binding |
| 5 | **p.Arg280Ile** (R280I) | 161517 | **-12.47** | **0.9996** (pathogenic) | Direct DNA contact (2.8 A) | DNA-Binding |

### Structural Context

- **R248P** and **R280I** make direct hydrogen bonds with DNA backbone in the
  p53-DNA co-crystal (PDB 1TUP, chain B). Proline substitution at R248
  eliminates the guanidinium group required for major groove recognition.
- **C176R** coordinates the structural zinc ion (Zn2+) that stabilizes the
  L2-L3 loop architecture. Loss of the thiolate ligand is expected to
  collapse the DNA-binding surface.
- **L257R** and **V157D** are buried in the hydrophobic core of the
  beta-sandwich fold. Introducing a charged residue into these tightly packed
  positions is predicted to destabilize the domain.

---

## Project Structure

```
TP53-VUS-Predict/
|-- data/
|   |-- tp53_uncertain_variants.csv         # Raw ClinVar VUS (1,211 variants)
|   |-- tp53_ai_predictions.csv             # ESM-2 scored variants
|   |-- tp53_alphamissense.csv              # AlphaMissense TP53 predictions (7,467)
|   |-- tp53_esm2_alphamissense_merged.csv  # Merged dataset (1,199 matched)
|   |-- tp53_wildtype_sequence.txt          # TP53 reference sequence (NP_000537.3)
|
|-- scripts/
|   |-- tp53_variant_analysis.py            # [Step 1] ClinVar data acquisition
|   |-- esm2_mutation_scorer.py             # ESM-2 scoring engine
|   |-- tp53_batch_scoring.py               # [Step 2] Batch LLR scoring
|   |-- tp53_domain_analysis.py             # [Step 3] Domain mapping & statistics
|   |-- alphamissense_crossref.py           # [Step 4] AlphaMissense cross-reference
|   |-- pymol_variant_views.py              # [Step 5] Structural analysis & PyMOL
|   |-- tp53_variant_views.pml              # PyMOL render script (generated)
|   |-- requirements.txt                    # Libraries and Dependencies for the Code to work
|
|-- figures/
|   |-- main/
|   |   |-- fig1_llr_distribution.png.png   # ESM-2 LLR score distribution
|   |   |-- fig2_model_concordance.png.png  # Cross-validation scatter plot
|   |   |-- fig3_domain_enrichment.png      # Domain-level statistical comparison
|   |   |-- fig4a_L257R_structural.png
|   |   |-- fig4b_V157D_structural.png
|   |   |-- fig4c_R248P_structural.png
|   |   |-- fig4d_C176R_structural.png
|   |   |-- fig4e_R280I_structural.png
|   |   |-- fig4f_overview_structural.png
|   |
|   |-- supplementary/
|       |-- figS1a_L257R_3d_context.png
|       |-- figS1b_V157D_3d_context.png
|       |-- figS1c_R248P_3d_context.png
|       |-- figS1d_C176R_3d_context.png
|       |-- figS1e_R280I_3d_context.png
|       |-- figS1f_overview_3d_context.png
|       |-- figS2a_L257R_pymol_v1.png
|       |-- figS2b_V157D_pymol_v1.png
|       |-- figS2c_R248P_pymol_v1.png
|       |-- figS2d_C176R_pymol_v1.png
|       |-- figS2e_R280I_pymol_v1.png
|       |-- figS2f_overview_pymol_v1.png
|
|-- vizualization/                          # Vizualization and Render Data Files
|   |-- 1tup.pdb
|   |-- pdb1tup.ent
|   |-- tp53_publication.pse
|   |-- tp53_variant_views.pml
|   |-- tp53_variants.pse
|
|-- manuscript/                             # LaTeX manuscript
|   |-- manuscript.tex
|   |-- manuscript.pdf
|
|-- article/                                # Article For Medium & Other Sites.
|   |-- article.md
|
|-- LICENSE                                 # MIT
|-- README.md
```

---

## Reproduction

### Prerequisites

- Python 3.10+
- CUDA-capable GPU recommended for ESM-2 (CPU works but is slower)
- ~2 GB disk space for AlphaMissense download
- [PyMOL](https://pymol.org/) for structural renders (conda recommended)

### 1. Clone and install dependencies

```bash
git clone https://github.com/mahaddev-x/TP53-VUS-Predict.git
cd TP53-VUS-Predict
pip install -r requirements.txt
```

### 2. Acquire ClinVar variants

```bash
python scripts/tp53_variant_analysis.py
```

Queries NCBI Entrez for TP53 missense VUS. Outputs `data/tp53_uncertain_variants.csv`.

### 3. Score variants with ESM-2

```bash
python scripts/tp53_batch_scoring.py
```

Loads `esm2_t33_650M_UR50D` (downloads ~2.5 GB on first run), scores all
1,211 variants by log-likelihood ratio. Outputs `data/tp53_ai_predictions.csv`.

### 4. Domain analysis

```bash
python scripts/tp53_domain_analysis.py
```

Maps each variant to TP53 functional domains, runs statistical comparisons
(Welch's t-test, Cohen's d), and generates domain comparison plots.

### 5. Cross-reference with AlphaMissense

```bash
python scripts/alphamissense_crossref.py
```

Downloads the AlphaMissense proteome file (~1.12 GB compressed) from Google
Cloud Storage, stream-filters for TP53 (UniProt P04637), merges with ESM-2
results, and produces the cross-validation scatter plot. The AlphaMissense
file is cached locally after the first download.

### 6. Structural visualization

```bash
python scripts/pymol_variant_views.py
```

Performs BioPython-based structural analysis on PDB 1TUP (downloaded
automatically) and generates a PyMOL `.pml` script. To produce ray-traced
renders, run the generated script in PyMOL:

```bash
# Using conda-installed PyMOL (recommended on Windows):
conda install -c conda-forge pymol-open-source
pymol -cq scripts/tp53_variant_views.pml
```

---

## Data Sources

| Source | Description | Access |
|--------|-------------|--------|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Variant clinical classifications | NCBI Entrez API |
| [ESM-2](https://github.com/facebookresearch/esm) | Protein language model | HuggingFace / fair-esm |
| [AlphaMissense](https://alphamissense.hegelab.org/) | Proteome-wide pathogenicity | Google Cloud Storage |
| [PDB 1TUP](https://www.rcsb.org/structure/1TUP) | p53 core domain-DNA complex | RCSB PDB |

---

## Citation

If you use this pipeline or findings in your work, please cite the underlying
models and databases:

- Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein
  structure with a language model. *Science*, 379(6637), 1123-1130.
- Cheng, J. et al. (2023). Accurate proteome-wide missense variant effect
  prediction with AlphaMissense. *Science*, 381(6664), eadg7492.
- Landrum, M.J. et al. (2024). ClinVar: improvements to accessing data.
  *Nucleic Acids Research*, 52(D1), D1265-D1273.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
