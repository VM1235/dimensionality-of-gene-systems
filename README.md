# Transcriptomic_Analysis — Figure 1A/1B reproduction (Hari et al., iScience 2025)

This folder reproduces **Figure 1A and 1B (panels i–iii)** from:

- Hari et al., iScience (2025), *“Low dimensionality of phenotypic space as an emergent property of coordinated teams in biological regulatory networks”*

It generates, for two binary phenotype axes:

- **(i)** gene–gene Pearson correlation heatmap (signature genes only)
- **(ii)** PC1 loadings bar plot
- **(iii)** PC1 vs PC2 scatter of cell lines, colored by a phenotype score

Outputs are saved as a single PNG.

---

## What’s in this folder

- **`reproduce_fig1AB.py`**: main script that loads data, subsets cell lines, computes correlations + PCA, and saves the figure.
- **`Genesets.csv`**: signature gene lists used in the figure.
  - `GeneSetID == SCLC` has **50 genes** in the order: first 25 **NE**, last 25 **non-NE**.
  - `GeneSetID == EMT` has **22 genes** (split into 11 epithelial + 11 mesenchymal by the script).
- **`Model.csv`**: DepMap model metadata used to filter cell lines for panel A/B.
- **`requirements.txt`**: Python dependencies.

---
## What is NOT there:
- **`OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`**: DepMap expression matrix (TPM log1p; protein-coding genes).
- **`virtual environment codes`**: python environment to ensure all libraries are available.
## Data → panel definitions (high level)

The code links files by **`ModelID`**:

- `Model.csv` provides annotations (e.g., `DepmapModelType`, `OncotreePrimaryDisease`).
- `OmicsExpression…csv` provides expression rows keyed by `ModelID` and gene columns named like `SYMBOL (EntrezID)`.

Panel subsets:

- **Figure 1A (SCLC)**: models where `DepmapModelType == "SCLC"`.
- **Figure 1B (EMT; non-SCLC cancers)**: models where:
  - `OncotreePrimaryDisease != "Non-Cancerous"` and
  - `DepmapModelType != "SCLC"`.

Expression filtering:

- Keeps only rows where `IsDefaultEntryForModel == "Yes"` (one profile per model).

Phenotype score (per cell line):

- `score = mean(team1 genes) − mean(team2 genes)`

PCA preprocessing:

- Uses **center-only** preprocessing (subtract per-gene mean; no variance scaling) to align with typical **R `prcomp` default** behavior.

Gene symbol compatibility:

- The script includes a small alias map for HGNC symbol updates:
  - `CYR61 → CCN1`
  - `FAM57B → CFAP97`

---

## Setup

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
pip install -r requirements.txt
```

Quick dependency check:

```bash
python -c "import pandas, numpy, matplotlib, seaborn, sklearn; print('deps ok')"
```

---

## Run

Matplotlib may need a writable cache directory. Use the included `.mplconfig/`:

```bash
source .venv/bin/activate
export MPLCONFIGDIR="$PWD/.mplconfig"
export MPLBACKEND=Agg
python reproduce_fig1AB.py
```

Expected output:

- `figure1AB_reproduction_final.png`

---

## Notes / caveats

- **Exact numeric match vs the paper**: the paper reports specific sample counts (e.g., 49 SCLC, 848 non-SCLC). Your results may differ because DepMap/CCLE releases and filters differ, which changes PCA variance explained (PC1%).
- **Gene symbols**: if additional genes go missing in other DepMap releases, you may need to add more alias mappings.
- **Runtime/memory**: the expression CSV is large; loading it can take time and requires enough RAM.

---

## Citation

If you use this code/output, please cite the iScience paper above and the DepMap/CCLE data sources used to obtain the expression and model metadata files.

