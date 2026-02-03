# I Used Google and Meta's AI to Hunt for Cancer Mutations — And Found Five.

*An open-source pipeline that cross-references two deep-learning models to reclassify 349 TP53 variants of uncertain significance.*

---

There are over 1,200 missense mutations in the TP53 gene — the single most mutated gene in human cancer — that sit in the NCBI ClinVar database labeled "Variant of Uncertain Significance." That label means clinicians can't act on them. Can't use them for genetic counseling. Can't use them to select patients for clinical trials. They're stuck in limbo.

I built a Python pipeline that takes two independent AI models — Meta's **ESM-2** protein language model and Google DeepMind's **AlphaMissense** pathogenicity predictor — and cross-references their scores across 1,211 TP53 variants of uncertain significance from ClinVar. When both models independently agree that a variant is damaging, that's not one opinion. That's convergent evidence from two fundamentally different architectures trained on fundamentally different data.

The result: **349 variants** flagged as high-risk by both models. Five of those variants scored so extreme that structural analysis confirms, residue by residue, exactly how they break the protein.

This article walks through the method, the code, the models, the results, and the structural biology behind the findings.

---

## The Problem: 1,211 Diagnostic Dead Ends

TP53 encodes the p53 protein, often called the "guardian of the genome" (Lane, 1992). When a cell's DNA is damaged, p53 activates repair, halts division, or triggers cell death. When TP53 is mutated, that checkpoint fails and tumors grow. TP53 loss-of-function mutations appear in over 50% of all solid tumors (Kandoth et al., 2013) and are the molecular basis of Li-Fraumeni syndrome, a hereditary condition where lifetime cancer risk exceeds 90% by age 60 (Villani et al., 2016).

The clinical problem is classification. Under the ACMG/AMP framework (Richards et al., 2015), a variant needs convergent evidence from multiple sources — population data, functional assays, segregation studies, computational prediction — to be classified as pathogenic. For rare variants seen only once or twice in clinical databases, that evidence may never accumulate through observation alone. Functional assays like the TP53 saturation mutagenesis screen by Giacomelli et al. (2018) help, but they're labor-intensive and expensive.

That leaves over 1,200 TP53 VUS in a diagnostic grey zone where no clinical action is taken.

---

## The Two Models: Different Data, Same Answer

The core idea of this project is **orthogonal validation** — using two AI models that derive their predictions from entirely different information sources, so that agreement between them carries stronger evidence than either model alone.

### Meta's ESM-2: A Language Model for Proteins

ESM-2 (`esm2_t33_650M_UR50D`) is a 650-million-parameter transformer model trained on 250 million protein sequences from the UniRef database (Lin et al., 2023). It works exactly like a large language model, but instead of predicting the next word in a sentence, it predicts which amino acid belongs at each position in a protein.

A protein is a chain of amino acids — a sequence written in a 20-letter alphabet. Evolution has been writing and editing these sequences for billions of years. Amino acid combinations that break a protein's function get selected out. Combinations that work get conserved. ESM-2 learns these evolutionary patterns from raw sequence data alone, with no structural information and no labels.

To score a mutation, you give ESM-2 the full wild-type TP53 sequence and compute the **log-likelihood ratio (LLR)**:

```
LLR = log P(mutant amino acid | sequence context) - log P(wild-type amino acid | sequence context)
```

A negative LLR means the mutant amino acid is less likely than the wild-type at that position — the model thinks evolution would reject it. The more negative the score, the more damaging the prediction. I stratified variants into four tiers: strongly damaging (LLR ≤ -4.0), likely damaging (-4.0 to -2.0), possibly damaging (-2.0 to -0.5), and likely neutral (> -0.5).

### Google DeepMind's AlphaMissense: Structure Meets Sequence

AlphaMissense (Cheng et al., 2023) takes a different approach. It combines **AlphaFold2-derived structural features** with sequence context and population frequency data to classify every possible human missense variant on a 0-to-1 pathogenicity scale. The published thresholds are: pathogenic (> 0.564), ambiguous (0.340-0.564), and benign (< 0.340).

Where ESM-2 sees only the sequence, AlphaMissense also sees the protein's predicted 3D shape. These are genuinely different information sources, which is what makes their agreement meaningful.

---

## The Pipeline: Five Scripts, End to End

The entire analysis runs as five Python scripts executed sequentially. Everything is in the GitHub repo.

### Step 1 — Variant Retrieval (`tp53_variant_analysis.py`)

Queries the NCBI Entrez API for all TP53 missense variants classified as "Uncertain Significance" in ClinVar (Landrum et al., 2018). Filters to single-nucleotide substitutions mapped to the canonical TP53 isoform (UniProt P04637, 393 amino acids). After deduplication by HGVS protein notation: **1,211 unique VUS**.

**Libraries:** `BioPython` (Entrez API access), `csv`, `urllib`

### Step 2 — ESM-2 Scoring (`tp53_batch_scoring.py`)

Loads the ESM-2 model via HuggingFace Transformers and scores all 1,211 variants in batch mode. Mutations at the same residue position share a single forward pass through the model for efficiency.

**Libraries:** `torch`, `transformers` (HuggingFace), `numpy`
**Hardware:** Consumer laptop with GPU acceleration (CUDA). Throughput: ~3 variants/second.
**Output:** `tp53_ai_predictions.csv` — each row contains the variant ID, ESM-2 LLR score, and severity tier.

12 variants were excluded where the ClinVar reference amino acid didn't match the UniProt wild-type sequence — likely annotation errors in ClinVar.

### Step 3 — AlphaMissense Cross-Reference (`alphamissense_crossref.py`)

Downloads the precomputed AlphaMissense proteome-wide predictions from the Zenodo repository (record 10813168, ~1.12 GB compressed). Stream-filters for TP53 (UniProt P04637) on the fly using `gzip` decompression to avoid loading the entire file into memory. Yields **7,467 scored substitutions** — every possible single-amino-acid change across all 393 TP53 residue positions.

Merges ESM-2 scores with AlphaMissense scores on the protein-level variant key. **1,199 of 1,211 VUS matched (99.0%).**

**Libraries:** `pandas`, `gzip`, `io`, `matplotlib`, `scipy` (for correlation statistics)

### Step 4 — Domain Analysis (`tp53_domain_analysis.py`)

Maps each variant to its TP53 functional domain — transactivation (residues 1-92), DNA-binding (102-292), or tetramerization (325-356) — and computes domain-level score enrichment using Welch's t-test and Cohen's d effect size.

**Libraries:** `pandas`, `scipy.stats`, `matplotlib`

### Step 5 — Structural Validation (`pymol_variant_views.py` and `pymol_publication_renders.py`)

Two complementary structural analysis scripts:

1. **BioPython spatial analysis** — Loads the PDB 1TUP crystal structure (Cho et al., 1994) of the p53 core domain bound to DNA at 2.2 angstrom resolution. Uses `BioPython.PDB.NeighborSearch` with a 4.0 angstrom contact distance threshold to identify neighboring residues, DNA contacts, and zinc ion proximity for each variant site. Generates 3D context plots.

2. **PyMOL headless rendering** — Runs inside PyMOL's headless mode (`pymol -cqr`) using the `cmd` API directly. For each variant, the wild-type residue is shown as green sticks and a computationally modeled mutant rotamer (applied via PyMOL's mutagenesis wizard) is overlaid as salmon sticks. All renders are ray-traced at 2,400 x 2,400 pixels with antialiasing level 4.

**Libraries:** `BioPython` (PDB parsing, NeighborSearch), `PyMOL` (`cmd` API), `matplotlib`, `numpy`

---

## Results

### ESM-2 Score Distribution

ESM-2 LLR scores for the 1,199 matched VUS ranged from **-13.88 to +4.21**. The distribution was left-skewed — consistent with the DNA-binding domain harboring a disproportionate fraction of damaging variants.

- Strongly damaging (LLR ≤ -4.0): **380 variants (31.7%)**
- Likely damaging: **232 variants (19.3%)**
- Possibly damaging: **282 variants (23.5%)**
- Likely neutral: **305 variants (25.4%)**

> **Figure 1** — *Distribution of ESM-2 log-likelihood ratio scores across 1,199 TP53 VUS. Dashed vertical lines mark classification thresholds.* (`figures/main/fig1_llr_distribution.png`)

### Model Concordance

The two models showed strong anti-correlation: **Pearson r = -0.706, Spearman rho = -0.714**. The negative correlation is expected — ESM-2 scores go more negative for damaging variants while AlphaMissense scores go higher. The key point: two models trained on different data representations independently converge on the same variants.

Applying dual thresholds:
- **349 variants (29.1%)** — high-risk by both models
- **438 variants (36.5%)** — low-risk by both models
- **412 variants (34.4%)** — discordant (models disagreed)

> **Figure 2** — *Scatter plot of ESM-2 LLR vs. AlphaMissense pathogenicity score for 1,199 matched TP53 VUS. Strong anti-correlation supports complementary predictive capacity.* (`figures/main/fig2_model_concordance.png`)

### Domain Enrichment

The DNA-binding domain (residues 102-292) showed significantly more negative mean ESM-2 scores than the transactivation or tetramerization domains — consistent with stronger evolutionary constraint on the DNA-binding interface. All five top-ranked variants mapped to this domain.

> **Figure 3** — *ESM-2 LLR score distribution by TP53 functional domain. The DNA-binding domain harbors the most severely scored variants.* (`figures/main/fig3_domain_enrichment.png`)

---

## The Five Variants

These are the five VUS with the most extreme concordant pathogenicity scores. All reside in the DNA-binding domain of p53 and exhibit structurally validated damage mechanisms.

| Variant | ESM-2 LLR | AlphaMissense | ClinVar ID | Structural Mechanism |
|---------|-----------|---------------|------------|---------------------|
| **L257R** | -13.88 | 0.9938 | 142134 | Hydrophobic core disruption |
| **V157D** | -13.69 | 0.9992 | 482231 | Hydrophobic core disruption |
| **R248P** | -12.63 | 0.9994 | 237954 | DNA minor groove contact lost |
| **C176R** | -12.53 | 0.9999 | 376573 | Zinc coordination abolished |
| **R280I** | -12.47 | 0.9996 | 161517 | DNA major groove contact lost |

Three distinct molecular mechanisms of p53 loss of function:

**1. Direct DNA contact loss (R248P, R280I).** Arginine 248 inserts directly into the DNA minor groove at ~3.4 angstroms, forming hydrogen bonds with the DNA backbone. Proline at this position eliminates all hydrogen-bonding capacity and introduces a rigid backbone kink. R248 is already the most frequently mutated residue in cancer — R248W and R248Q are among the six most common TP53 mutations (Bouaoun et al., 2016) — but R248P has not been classified as pathogenic in ClinVar despite affecting the same critical contact residue. Arginine 280 forms a direct hydrogen bond with a guanine base in the DNA major groove at ~2.8 angstroms; isoleucine is hydrophobic and abolishes this contact entirely.

**2. Zinc coordination abolished (C176R).** Cysteine 176 is one of four residues (C176, H179, C238, C242) that coordinate the structural Zn2+ ion essential for folding the L2 and L3 loops of the DNA-binding domain (Cho et al., 1994; Bullock et al., 2000). The cysteine thiolate provides a zinc ligand at ~2.3 angstroms. Arginine cannot coordinate zinc. AlphaMissense assigned its highest score in the entire dataset — **0.9999** — to this variant.

**3. Hydrophobic core disruption (L257R, V157D).** Both leucine 257 and valine 157 are buried within the hydrophobic beta-sandwich core. Substitution to charged residues (arginine, aspartate) introduces electrostatic repulsion and packing cavities in a tightly packed nonpolar environment. This is a well-established mechanism of p53 thermodynamic destabilization — such mutations reduce the melting temperature of the DNA-binding domain and accelerate unfolding at body temperature (Bullock et al., 2000).

> **Figure 4** — *Publication-quality PyMOL renders of all five variants. Wild-type residues (green sticks) overlaid with computationally modeled mutant rotamers (salmon sticks). Neighboring residues within 4.0 angstroms shown as light-blue sticks. DNA contacts: red dashes. Zinc coordination: slate dashes. All panels ray-traced at 2,400 x 2,400 pixels, 300 DPI.* (`figures/main/fig4a-f`)
>
> (a) L257R — charged arginine disrupts the hydrophobic core
> (b) V157D — aspartate introduces charge and a packing void
> (c) R248P — proline eliminates DNA minor groove contact
> (d) C176R — arginine abolishes Zn2+ coordination
> (e) R280I — isoleucine abolishes DNA major groove hydrogen bond
> (f) Overview — all five sites mapped on the p53-DNA complex

---

## Clinical Implications

**Li-Fraumeni syndrome.** Germline TP53 pathogenic variants trigger intensive cancer surveillance — annual whole-body MRI, breast MRI, abdominal ultrasound (Villani et al., 2016). When a variant is labeled VUS, patients and their families can't access this surveillance. Reclassifying high-confidence VUS would directly enable cascade testing and early detection.

**Somatic tumor profiling.** TP53 status predicts chemotherapy resistance in chronic lymphocytic leukemia (Zenz et al., 2010) and adverse prognosis in myelodysplastic syndromes (Bejar et al., 2011). VUS labels reduce the precision of molecular stratification.

**Targeted therapy selection.** New drugs — APR-246 (eprenetapopt) and PC14586 (rezatapopt) — are in clinical trials to reactivate mutant p53 (Chen et al., 2021). Patient enrollment requires classified mutations. Variants that destabilize the protein fold (L257R, V157D) may respond to fold-stabilizing compounds; those that abolish DNA contact (R248P, R280I) may require different strategies.

---

## Limitations

These are computational predictions — they do not constitute functional evidence under ACMG/AMP criteria and cannot alone reclassify a variant beyond "supporting" evidence strength. The structural analysis is based on a single crystal structure (PDB 1TUP) capturing one conformational state. The ESM-2 and AlphaMissense thresholds were adopted from their original publications and have not been independently calibrated on a TP53-specific truth set. All 349 high-risk classifications should be treated as provisional until validated by functional assays — yeast-based transcription assays, thermal stability measurements, or DNA-binding electrophoretic mobility shift assays.

---

## Tech Stack Summary

| Component | Tool / Library | Purpose |
|-----------|---------------|---------|
| Data retrieval | BioPython + NCBI Entrez API | ClinVar VUS download |
| Sequence scoring | ESM-2 (650M params) via HuggingFace Transformers + PyTorch CUDA | Zero-shot LLR variant effect prediction |
| Structure scoring | AlphaMissense via Zenodo | Precomputed pathogenicity scores |
| Cross-referencing | pandas, scipy | Merge, correlation, concordance |
| Domain mapping | pandas, scipy.stats | Welch's t-test, Cohen's d |
| Spatial analysis | BioPython PDB + NeighborSearch | 3D contact geometry |
| 3D rendering | PyMOL 3.1.0 headless (`cmd` API) | Publication ray-traced renders |
| Visualization | matplotlib | Score distributions, scatter plots |
| Manuscript | LaTeX (Tectonic engine) | Full scientific paper |
| Environment | Python 3.10+, micromamba, CUDA | Package and environment management |

---

## The Repository

The full pipeline — every script, every scored variant, every figure, and the complete LaTeX manuscript — is open source:

**[github.com/mahaddev-x/TP53-AI-Predict](https://github.com/mahaddev-x/TP53-AI-Predict)**

What's inside:
- `scripts/tp53_variant_analysis.py` — ClinVar data retrieval
- `scripts/tp53_batch_scoring.py` — ESM-2 batch scoring engine
- `scripts/alphamissense_crossref.py` — AlphaMissense cross-reference and concordance
- `scripts/tp53_domain_analysis.py` — Domain-level enrichment statistics
- `scripts/pymol_variant_views.py` — BioPython structural analysis + PyMOL script generation
- `scripts/pymol_publication_renders.py` — PyMOL headless publication renders
- `data/` — All scored variant CSVs
- `figures/main/` — Publication figures (Figs 1-4)
- `figures/supplementary/` — BioPython 3D context plots (Fig S1) and initial PyMOL renders (Fig S2)
- `manuscript/` — Full LaTeX manuscript with methods, results, and references

### Where You Can Contribute

- **Other cancer genes** — extend the pipeline to BRCA1, BRCA2, PTEN, RB1
- **Functional assay integration** — add MAVE/DMS datasets as a third validation axis
- **Threshold calibration** — benchmark against known pathogenic/benign TP53 variants
- **Web interface** — build a lookup tool where clinicians can query any TP53 variant
- **Additional predictors** — integrate EVE (Frazer et al., 2021), REVEL (Ioannidis et al., 2016), or molecular dynamics simulations

If you work in computational biology, clinical genetics, or cancer genomics — or if you're a developer who wants to contribute to open-source cancer research — I'd like to hear from you.

**Email:** mahaddevx@gmail.com
**GitHub:** [@mahaddev-x](https://github.com/mahaddev-x)

---

## References

- Bejar, R. et al. (2011) Clinical effect of point mutations in myelodysplastic syndromes. *N. Engl. J. Med.* 364(26):2496-2506.
- Bouaoun, L. et al. (2016) TP53 variations in human cancers: new lessons from the IARC TP53 Database. *Hum. Mutat.* 37(9):865-876.
- Bullock, A. N., Henckel, J. & Fersht, A. R. (2000) Quantitative analysis of residual folding and DNA binding in mutant p53 core domain. *Oncogene* 19(10):1245-1256.
- Chen, S. et al. (2021) Small molecule therapeutics for TP53-mutant cancers. *Trends Pharmacol. Sci.* 42(12):1049-1062.
- Cheng, J. et al. (2023) Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381(6664):eadg7492.
- Cho, Y. et al. (1994) Crystal structure of a p53 tumor suppressor-DNA complex. *Science* 265(5170):346-355.
- Cock, P. J. A. et al. (2009) Biopython: freely available Python tools for computational molecular biology. *Bioinformatics* 25(11):1422-1423.
- Frazer, J. et al. (2021) Disease variant prediction with deep generative models of evolutionary data. *Nature* 599(7883):91-95.
- Giacomelli, A. O. et al. (2018) Mutational processes shape the landscape of TP53 mutations in human cancer. *Nat. Genet.* 50(10):1381-1387.
- Ioannidis, N. M. et al. (2016) REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am. J. Hum. Genet.* 99(4):877-885.
- Kandoth, C. et al. (2013) Mutational landscape and significance across 12 major cancer types. *Nature* 502(7471):333-339.
- Landrum, M. J. et al. (2018) ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res.* 46(D1):D1062-D1067.
- Lane, D. P. (1992) p53, guardian of the genome. *Nature* 358(6381):15-16.
- Lin, Z. et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379(6637):1123-1130.
- Richards, S. et al. (2015) Standards and guidelines for the interpretation of sequence variants. *Genet. Med.* 17(5):405-424.
- Villani, A. et al. (2016) Biochemical and imaging surveillance in germline TP53 mutation carriers with Li-Fraumeni syndrome. *Lancet Oncol.* 17(9):1295-1305.
- Zenz, T. et al. (2010) TP53 mutation and survival in chronic lymphocytic leukemia. *J. Clin. Oncol.* 28(29):4473-4479.

---

*Mahad Asif — 9th-grade computer science student. This project was developed as independent computational biology research using publicly available data and open-source AI models.*

---

**Tags:** `#ArtificialIntelligence` `#Cancer` `#Genomics` `#DeepLearning` `#TP53` `#Bioinformatics` `#OpenSource` `#ProteinLanguageModels` `#PrecisionMedicine` `#ClinVar` `#ESM2` `#AlphaMissense` `#ComputationalBiology`
