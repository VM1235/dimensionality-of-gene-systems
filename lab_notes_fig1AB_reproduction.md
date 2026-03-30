# Lab Notes: Reproducing Figure 1A–B from Hari et al. iScience 2025
**"Low dimensionality of phenotypic space as an emergent property of coordinated teams in biological regulatory networks"**

Author of these notes: [Your name]  
Date: [Today]  
Script: `reproduce_fig1AB.py`  
Output: `figure1AB_reproduction_final.png`

---

## 1. Workflow (Step-by-Step)

### Step 0 — Environment setup
Install dependencies from `requirements.txt`:
```
pandas>=2.0, numpy>=1.24, matplotlib>=3.7,
seaborn>=0.13, scikit-learn>=1.3, scipy>=1.11
```
No custom or specialized package is required beyond standard scientific Python. The `sklearn.decomposition.PCA` is used for the actual PCA calculation. Notably, this is Python-based; the paper used R's `prcomp`.

---

### Step 1 — Inputs and what each file contributes

| File | Role |
|---|---|
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | Main expression matrix from DepMap. Rows = samples (one per model after filtering); columns = genes in `"SYMBOL (EntrezID)"` format; values = log₁p(TPM). |
| `Model.csv` | DepMap sample metadata. Provides `DepmapModelType` (for SCLC/non-SCLC classification) and `OncotreePrimaryDisease` (for cancer/non-cancer filtering). |
| `Genesets.csv` | Curated gene lists with a `GeneSetID` column. Contains 50 SCLC genes (Lissa et al. 2022), 22 EMT genes (Aiello et al. 2018), plus Alveolar and FUGUE sets not used in these two panels. |
| `reproduce_fig1AB.py` | End-to-end pipeline: loads all inputs, filters samples, builds matrices, runs PCA, draws all six sub-panels, and saves the PNG. |

---

### Step 2 — Gene list construction

**SCLC signature (Panel A):**  
`load_genesets()` reads all rows in `Genesets.csv` where `GeneSetID == "SCLC"`. It expects exactly 50 genes. The function enforces a hard split by row position: rows 0–24 are treated as NE genes; rows 25–49 as non-NE genes. This ordering is assumed to reflect the Lissa et al. 2022 gene ordering — the code does not validate this biologically, it enforces it structurally.

**EMT signature (Panel B):**  
All 22 rows where `GeneSetID == "EMT"` are loaded. The code then validates that this set equals `set(EMT_E) | set(EMT_M)` where both lists are hardcoded in the script based on Aiello et al. 2018. The E/M split for plotting uses those hardcoded lists — the ordering within the heatmap block is therefore determined by the `EMT_E` and `EMT_M` Python lists in the script, not by the CSV row order.

**Gene symbol alias mapping:**  
Two SCLC genes required alias lookup due to HGNC symbol updates:
- `FAM57B` → `CFAP97` (approved symbol in current DepMap)
- `CYR61` → `CCN1` (approved symbol in current DepMap)

These are silently resolved inside `build_matrix()` via a hardcoded `ALIASES` dict. No other aliases are handled; if any other signature gene had been renamed, it would appear in the "genes not found" warning and be silently dropped.

---

### Step 3 — Sample (cell line) filtering

**Panel A — SCLC cell lines:**  
`get_model_ids(meta, "SCLC")` returns all `ModelID` values where `DepmapModelType == "SCLC"`. In our `Model.csv`, this yields **81 SCLC models** (vs. 49 in the paper).

**Panel B — Non-SCLC cancer cell lines:**  
The filter is:  
`OncotreePrimaryDisease != "Non-Cancerous"` **AND** `DepmapModelType != "SCLC"`.  
This yields **1,900 non-SCLC cancer models** (vs. 848 in the paper; the paper used only breast cancer + non-SCLC cell lines based on the STAR Methods text).

**Default-entry deduplication:**  
Inside `load_ccle_expression()`, before any panel-specific filtering, the script keeps only rows where `IsDefaultEntryForModel == "Yes"`. This removes technical replicates or multi-pass profiling entries, retaining one expression profile per DepMap model ID. The number of rows kept depends on the DepMap release.

**Intersection with expression data:**  
`build_matrix()` intersects the model ID set from metadata with the rows actually present in the expression file. Any model in `Model.csv` that is absent from the expression file is silently dropped.

---

### Step 4 — Expression matrix construction

For each panel, `build_matrix()` subsets the full expression matrix to:
- Rows: cell lines passing the group filter (above)
- Columns: signature genes found in the expression file (after alias resolution)

Column names in the expression file follow `"SYMBOL (EntrezID)"` format. The function strips the ` (EntrezID)` suffix and renames columns to clean gene symbols for all downstream use. Rows that are entirely NaN after subsetting are dropped (this would only happen if a sample had no expression values at all for the selected genes, which is rare).

The resulting values are log₁p-transformed TPM, as provided by DepMap — no further normalization is applied at this step.

---

### Step 5 — Phenotypic score calculation

For each cell line:
```
score = mean(expression of team1 genes) − mean(expression of team2 genes)
```
- Panel A: team1 = NE genes, team2 = non-NE genes → positive score = more NE-like
- Panel B: team1 = Epithelial genes, team2 = Mesenchymal genes → positive score = more E-like

The score is computed on the **raw (unnormalized, un-mean-centered) log₁p TPM values**. This is used only for colorizing the PCA scatter plot (panel iii); it does not enter the PCA computation itself.

---

### Step 6 — PCA preprocessing (critical choice)

The script intentionally uses **center-only** PCA:
```python
X = mat_ord.values - mat_ord.values.mean(axis=0, keepdims=True)
```
This subtracts the per-gene mean across samples, but does **not** divide by the per-gene standard deviation. This matches R's `prcomp(data, center=TRUE, scale.=FALSE)` default, which the paper states was used.

If `StandardScaler` (mean=0, std=1) had been applied instead, each gene would be treated as equally variable regardless of its actual expression variance, which would give different PC1 variance percentages and potentially different gene loading magnitudes.

`sklearn.decomposition.PCA` is then run on this mean-centered matrix with `n_components = min(10, n_samples-1, n_genes)`.

---

### Step 7 — Visualization: three sub-panels per row

**Panel (i) — Gene-gene Pearson correlation heatmap:**  
`corr = mat_ord.corr(method="pearson")` computes pairwise gene–gene correlations across all samples in that panel. Genes are ordered as `team1 block + team2 block` (NE then non-NE, or E then M), which pre-arranges the expected 2-block structure. The colormap is `RdBu_r` (diverging, centered at 0, range −1 to +1). Team bracket labels are annotated outside the y-axis.

**Panel (ii) — PC1 loadings bar chart:**  
`pca.components_[0]` gives the first PC loading vector (one value per gene). Bars are colored red (team1) or blue (team2) by hardcoded membership. The y-axis is ordered to match the heatmap (team1 top, team2 bottom) and then inverted so that team1 appears at the top of the chart. The title reports the % variance explained by PC1.

**Panel (iii) — PCA scatter plot:**  
Cell lines are plotted at their (PC1 score, PC2 score) coordinates and colored by phenotypic score using the `RdBu_r` diverging colormap. Axis labels include variance percentages. No clustering, contour lines, or marginal distributions are added — this is a plain scatter.

---

### Step 8 — Output

The figure is saved as `figure1AB_reproduction_final.png` at 200 dpi using `bbox_inches="tight"`. The layout is a 2-row × 3-column grid with explicit `hspace`, `wspace`, `left`, `right`, `top`, and `bottom` spacing parameters set manually.

---

## 2. My Observations

- **Two-block correlation structure is clearly visible in both panels.** In Panel A, the NE genes (upper-left block) show strong positive mutual correlation (red) and strong negative cross-correlation with non-NE genes (lower-right block, also internally correlated). Panel B shows the same E vs. M block structure. This is qualitatively consistent with Fig. 1A(i) and 1B(i) in the paper.

- **PC1 variance in our Panel A (SCLC) is ~46.74%**, considerably lower than the paper's reported 56.97%. This is the most numerically significant discrepancy and is likely attributable to our larger, more heterogeneous SCLC sample (81 models vs. 49), a different DepMap release, or both. A more heterogeneous sample dilutes the binary NE/non-NE signal and reduces PC1 dominance.

- **PC1 variance in our Panel B (EMT) is ~48.78%**, compared to the paper's 46.33%. This is actually slightly higher than the paper, which is plausible given that using a much larger and broader non-SCLC cancer set (1,900 lines) might still preserve the strong E/M dichotomy across diverse cancer types, since the EMT axis is near-universal.

- **PC2 variance is ~10.07% for Panel A and ~17.95% for Panel B in our run.** The paper reports 8.89% and 19.79% respectively. The relative ordering is preserved (PC2 is proportionally larger in the EMT panel), but our absolute numbers differ slightly.

- **Gene loading signs in Panel A correctly separate NE (positive PC1 loadings) from non-NE (negative PC1 loadings).** The bar chart shows NE genes as rightward (positive) bars in red and non-NE genes as leftward (negative) bars in blue. This matches the paper's Fig. 1A(ii) qualitatively. The specific gene magnitudes and ordering differ because our sample set differs.

- **Gene loading signs in Panel B correctly separate Epithelial (positive) from Mesenchymal (negative) genes.** CDH1, CLDN2, EPCAM etc. load positively; TNC, SPARC, POSTN etc. load negatively. This is consistent with the paper's Fig. 1B(ii).

- **The PCA scatter for Panel B (EMT, ~1,900 lines) shows a clear gradient from deep blue (M-like) to deep red (E-like) along PC1.** The spread along PC2 likely captures additional cancer-type-level variation not specific to EMT. The paper's scatter for the same panel, with ~848 lines, appeared qualitatively similar in gradient structure.

- **The PCA scatter for Panel A (SCLC, 81 lines) shows loose clustering**, with some deeply blue (non-NE) and some deeply red (NE) points, consistent with the binary NE/non-NE classification, but with intermediate/hybrid points visible. Our scatter is more populated than the paper's (49 lines) and thus looks slightly denser.

- **The alias resolution for CYR61→CCN1 and FAM57B→CFAP97 worked silently and correctly.** No "genes not found" warning was raised for the SCLC panel, indicating that all 50 SCLC genes were resolved either directly or via alias. Without this mapping, 2 out of 50 SCLC genes (4%) would be silently missing, potentially biasing loadings and correlation structure.

- **The EMT gene set passes the set equality check exactly.** All 22 genes in Genesets.csv match the union of the hardcoded `EMT_E` and `EMT_M` lists, so no mismatches or typos are present in our EMT gene list.

- **The model count discrepancy is large for SCLC (81 vs. 49).** Our `Model.csv` contains 81 SCLC-annotated entries, which is 65% more than the paper reports. This strongly suggests our `Model.csv` is from a newer DepMap release that has added SCLC lines since the paper's data freeze.

- **The non-SCLC cancer count is also substantially larger (1,900 vs. 848).** The paper specifically states it used "non-SCLC and breast cancer cell lines," implying a more restricted subset definition that our simple `DepmapModelType != "SCLC"` filter does not capture. The paper likely filtered to specific cancer lineages, not all non-SCLC cancers globally.

---

## 3. My Questions

- **Which exact DepMap release did the paper use?** The STAR Methods say "CCLE" (Barretina et al. 2012 cite) but CCLE expression data has been hosted and versioned through DepMap since ~2019. The specific DepMap release version (e.g., 22Q2, 23Q2) would tell us exactly how many SCLC and non-SCLC lines were available. Can the authors share the release tag, or is it in the GitHub repository at `aashnasaxena/Teams-PC1`?

- **The paper says 49 SCLC lines and 848 non-SCLC lines — what exact filters generated those numbers?** For non-SCLC, did they filter to specific Oncotree lineages (e.g., lung + breast)? The STAR Methods say "non-SCLC (848) and breast cancer (62) cell lines," which reads as if those are two separate subsets — did the paper use the union, or analyze them separately?

- **Is `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` from DepMap equivalent to the CCLE TPM matrix?** DepMap reprocesses CCLE RNA-seq periodically. The log₁p transformation should be standard, but the alignment pipeline, gene model, and TPM computation method may differ between releases. Can we confirm the expression matrix provenance matches?

- **Was `prcomp` called with `scale.=FALSE` confirmed?** The STAR Methods say "PCA analysis was performed using the prcomp function from R 4.1" but do not explicitly state the `scale.` parameter. We assumed `FALSE` (center-only) to match the R default, but `scale.=TRUE` would give substantially different PC1 variance. Is this documented in the GitHub code?

- **Why does our center-only PCA in Python yield a different result than R's prcomp on the same centered matrix?** `sklearn.PCA` uses SVD under the hood, as does `prcomp`. They should be numerically identical for the same input matrix. If the numbers still differ, it must be sample composition, not the algorithm. Are the per-gene means computed over the same set of samples in both cases?

- **How were samples with any missing gene expression values handled?** Our script drops rows that are entirely NaN, but partial NaN within a row (a gene missing for one sample) would propagate through Pearson correlation. Did the paper impute, drop, or mask missing values? Is there a meaningful fraction of missing values in the SCLC expression data for these 50 genes?

- **Why does `IsDefaultEntryForModel` exist, and how many rows does it remove?** If a model has multiple expression profiles (e.g., from different assay batches), choosing only the default entry makes results batch-dependent. How many models in the current release have non-default entries?

- **The SCLC geneset contains `FAM57B` (an alias), but should the gene list have been updated to `CFAP97` before submission?** Is the Genesets.csv meant to store historical aliases or current HGNC symbols? If someone reproduces this on a future DepMap release where `CCN1` has again been renamed, the alias mapping would silently fail.

- **What does the PC2 axis capture in Panel B?** With ~1,900 non-SCLC cancer lines spanning all cancer types, PC2 likely captures cancer-lineage identity (e.g., carcinoma vs. sarcoma) rather than anything EMT-specific. Did the paper comment on this, or did the restricted breast cancer + non-SCLC subset in the paper reduce this confound?

- **Is the phenotypic score formula in the paper identical to ours?** The STAR Methods define it as `score = (Σ_{i∈P1} e_i / #P1) − (Σ_{i∈P2} e_i / #P2)`, which is mean(P1) − mean(P2). Our implementation matches this. But the paper computes this on raw TPM log₁p, not on mean-centered values. Is that intended, or was it supposed to be on z-scored expression?

- **Could the 2-alias issue in SCLC genes indicate that more genes have been renamed and we just don't know about them?** We only added aliases for the two we happened to discover. A systematic check against the current HGNC database would be more rigorous. How many of the 50 SCLC genes are at risk of symbol changes?

- **How sensitive are the PC1 variance numbers to the exact set of samples included?** Given that going from 49 to 81 SCLC lines changes PC1 variance from 56.97% to ~46.74%, this is highly sensitive to sample composition. Does the paper's GitHub code include a sample list or model ID file for exact reproducibility?

---

## 4. Caveats / Limitations

- **DepMap release mismatch is unverified and likely.** We have no documentation of which DepMap release was used in the paper. Our `Model.csv` contains 2,132 models, and our SCLC count (81) substantially exceeds the paper's (49). Comparing figures is informative, but numerical results should not be treated as validated reproductions until the release is confirmed.

- **The non-SCLC panel uses a much broader sample pool than the paper.** The paper used ~848 cell lines (likely restricted to specific lineages); we use ~1,900. This inflates sample size and includes cancer types with highly variable EMT status, which may add noise to the correlation structure and PCA while paradoxically still yielding high PC1 variance due to the universality of the EMT axis.

- **`IsDefaultEntryForModel` filtering is release-dependent.** The set of "default" entries may change between DepMap releases if technical QC criteria are updated. Our filtering is consistent within a single run, but is not guaranteed to reproduce the exact sample set across releases.

- **Two gene aliases are hardcoded without a systematic audit.** `CYR61→CCN1` and `FAM57B→CFAP97` were identified manually. No automated cross-check against HGNC was performed. If any of the remaining 48 SCLC genes or 22 EMT genes have been renamed, they would be silently dropped with a console warning that is easy to miss.

- **Gene order within the NE and non-NE blocks is determined by row order in Genesets.csv, which is not validated biologically.** The heatmap block structure looks correct visually because genes are pre-sorted by team identity, but the within-block ordering (which gene appears first in the NE block, etc.) is arbitrary and will differ from the paper's reordering (which appears to sort by loading magnitude or hierarchical clustering within teams).

- **PCA preprocessing (center-only) is inferred from R documentation, not confirmed from the paper's code.** If the paper's R script used `scale.=TRUE` (rare for expression data but possible), our center-only approach would underestimate genes with low variance and overestimate PC1 variance for high-variance genes, or vice versa.

- **The phenotypic score is computed on raw log₁p TPM values, not on PCA-preprocessed (mean-centered) values.** This means the colorbar in panel (iii) reflects absolute expression differences between teams, which are affected by differences in mean expression level across gene groups unrelated to phenotypic state. This may explain why our colorbar range differs from the paper's.

- **sklearn's PCA and R's prcomp can give sign-flipped loadings.** PCA eigenvectors are defined up to a sign flip: PC1 could have NE genes loading positively or negatively depending on the implementation. Our script does not enforce a canonical sign convention. If the sign flips, the scatter coloring and loading bar chart would appear inverted. Visual inspection confirms that our signs appear consistent with the paper for this run, but this should be explicitly enforced in a production script.

- **No statistical testing is performed.** The paper computes Spearman correlations, p-values, and ANOVA tests. Our reproduction computes only the descriptive visualizations (correlation heatmap, PCA loadings, scatter). We do not validate whether the team structure is statistically significant or whether PC1 variance is significantly higher than for random genesets.

- **The `dropna(how="all")` row filter in `build_matrix()` is permissive.** A sample with partial NaN values (some genes missing) is retained and would introduce NaN-driven artifacts in Pearson correlation (Pandas `corr()` computes pairwise-complete-case correlations, which is correct behavior but may reduce effective sample sizes per gene pair unpredictably).

- **No versioning of Genesets.csv.** The gene list file has no metadata about its source, creation date, or the exact Lissa et al. 2022 supplementary table it was derived from. If the file was constructed manually, transcription errors cannot be excluded without independently regenerating it from the source paper.

- **Figures are visual reproductions, not quantitative validation.** We can confirm that the qualitative structure (2-block correlation, opposite-sign loadings, gradient scatter) is reproduced. We cannot confirm numerical equivalence (exact PC1 %, exact sample counts) without the paper's exact data and code environment.

---

## 5. Paper-Alignment Check

| What the paper says / implies | What our pipeline does |
|---|---|
| Uses CCLE data (Barretina et al. 2012); cites DepMap for download | Uses DepMap `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` — consistent in format, release unverified |
| 49 SCLC cell lines for Panel A | Our pipeline yields **81** SCLC lines (DepmapModelType == "SCLC"); likely a newer DepMap release |
| 848 non-SCLC + 62 breast cancer lines for Panel B | Our pipeline yields **~1,900** non-SCLC cancer lines; we do not restrict to specific lineages |
| SCLC geneset: 50 genes from Lissa et al. 2022, split NE/non-NE | 50 genes from Genesets.csv; first 25 = NE, last 25 = non-NE by row order |
| EMT geneset: 22 genes from Aiello et al. 2018 | 22 genes from Genesets.csv; E/M split enforced by hardcoded lists matching Aiello et al. |
| Gene-gene Pearson correlation matrix (heatmap, panel i) | Computed identically: `DataFrame.corr(method="pearson")` |
| PCA with R `prcomp`, center=TRUE (implied default) | `sklearn.PCA` on manually mean-centered matrix; numerically equivalent to prcomp(scale.=FALSE) |
| PC1 variance Panel A: 56.97% | Our Panel A: **~46.74%** (discrepancy attributed to sample set difference) |
| PC1 variance Panel B: 46.33% | Our Panel B: **~48.78%** (slightly higher; plausible given broader sample pool) |
| Phenotypic score = mean(team1) − mean(team2) per sample | Implemented identically in `phenotypic_score()` |
| Gene loadings show opposite signs for the two teams | Confirmed in our output: NE genes load positively, non-NE negatively (Panel A); E positive, M negative (Panel B) |
| 2-block positive/negative correlation structure | Confirmed in both heatmaps; qualitatively matches paper |
| Genes `CYR61` and `FAM57B` used in SCLC signature | Resolved to `CCN1` and `CFAP97` respectively via alias mapping; not documented in the paper's methods |
| Paper's GitHub: `aashnasaxena/Teams-PC1` | Not accessed in this pipeline; would be required to confirm exact sample lists and R preprocessing parameters |
