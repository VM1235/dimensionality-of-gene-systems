# Intrinsic Dimension of Gene Expression During Cell Differentiation

## Quick Overview

This project reproduces **Figures 1A and 1B** from the paper *"The intrinsic dimension of gene expression during cell differentiation"*. 

It analyzes how gene expression complexity changes during embryonic development in **zebrafish** and **mouse** using single-cell RNA-seq data. The key finding: **as cells differentiate, the intrinsic dimension (ID) of gene expression decreases**, indicating cells shift from pluripotent (many genes expressed) to specialized (few genes expressed) states.

---

## What This Does

### Main Analysis
The notebook (`Panel1_figures.ipynb`) contains:

1. **Zebrafish Embryogenesis Analysis**
   - Dataset: 36,749 cells across 7 developmental stages (4-24 hours post-fertilization)
   - Measures how gene expression complexity decreases during early development
   - Generates Figure 1A

2. **Mouse Gastrulation Analysis**
   - Dataset: 108,839 cells across 9 developmental stages (E6.5-E8.5 embryonic days)
   - Tracks complexity changes during gastrulation and early organogenesis
   - Generates Figure 1B

### The Method
- **Algorithm**: 2-Nearest Neighbor (2-NN) method for intrinsic dimension calculation
- **Validation**: 3 independent subsamples per stage for robustness
- **Metric**: ID-score normalized to 0-1 scale (0 = low complexity, 1 = high complexity)

---

## File Structure

```
.
├── README.md                                      # This file
├── EXPLANATION.md                                 # Detailed technical explanation
├── Panel1_figures.ipynb                          # Main Jupyter notebook
├── requirements.txt                              # Python dependencies
├── .venv/                                        # Virtual environment
├── My_libs/
│   ├── IDmeter.py                               # Intrinsic dimension calculation
│   ├── Plot_figures.py                          # Plotting functions
│   └── File_handler.py                          # File I/O utilities
├── ZebrafishEmbryo_Wagner_formatted.h5ad        # Zebrafish data (preprocessed)
└── MouseGastrulation_formatted.h5ad             # Mouse data (preprocessed)
```

---

## Data Sources

### Preprocessed Data (h5ad format)
The formatted single-cell datasets are available here:
📥 **[Download preprocessed data](https://drive.google.com/drive/folders/1bm69GFaq8lcXRjAtxbgQIi2j_H_cX6bi)**
*(Provided by the authors in their GitHub repository)*

### Original Raw Data
- **Zebrafish**: GEO repository [GSE112294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112294)
- **Mouse**: E-MTAB repository [6967](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6967)

---

## Quick Start

### 1. Install Dependencies
```bash
cd /Users/apple/Downloads/Project/ID/self
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Run the Notebook
Open `Panel1_figures.ipynb` in Jupyter or VS Code and run cells in order:
- Cell 1: Install packages
- Cell 2: Import libraries
- Cells 3-13: Zebrafish analysis and plot
- Cells 14-32: Mouse gastrulation analysis and plot

### 3. View Results
Plots are generated inline in the notebook showing ID-score trends during development.

---


## Results Summary

### Zebrafish (04hpf → 24hpf)
```
ID-score: 0.65 → 0.42 → 0.35 → 0.30 → 0.28 → 0.15 → 0.04
Trend: Steady decrease in complexity during early embryogenesis
```

### Mouse (E6.5 → E8.5)
```
ID-score: 0.90 → 0.88 → 0.67 → 0.55 → 0.52 → 0.17 → 0.33 → 0.38 → 0.02
Trend: Overall decrease with dynamic changes during gastrulation
```

Both show: **Differentiation = Dimensionality Reduction**

---

## Main Libraries Used

| Library | Purpose |
|---------|---------|
| **scanpy** | Single-cell RNA-seq analysis, data loading |
| **dadapy** | Intrinsic dimension calculation (2-NN method) |
| **numpy** | Numerical computing |
| **pandas** | Data manipulation |
| **matplotlib/seaborn** | Publication-quality plotting |

---


### What Each Step Does

1. **Load data**: Read h5ad file and extract expression matrix + metadata
2. **Count cells**: Determine how many cells in each stage
3. **Calculate ID**: Apply 2-NN algorithm to measure complexity
4. **Plot results**: Visualize temporal trend with error bars

---


## Key Findings

✅ **Intrinsic dimension decreases with differentiation** (both zebrafish and mouse)

✅ **Principle is evolutionarily conserved** (works across species)

✅ **Method is robust** (consistent across 3 independent subsamples)

✅ **Quantitative measure of development** (provides numerical metric for developmental stage)

