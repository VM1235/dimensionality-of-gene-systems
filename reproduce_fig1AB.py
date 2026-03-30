"""
Reproduction of Figure 1A and 1B from:
Hari et al., iScience 2025
"Low dimensionality of phenotypic space as an emergent property
 of coordinated teams in biological regulatory networks"

Figure 1A: SCLC gene signature in CCLE SCLC cell lines
  (i)   Gene-gene Pearson correlation heatmap  →  should show 2-block structure
  (ii)  PC1 loadings bar chart                 →  NE and non-NE bars on opposite sides
  (iii) PCA scatter of cell lines              →  colored by phenotypic score (viridis)

Figure 1B: EMT gene signature in CCLE non-SCLC cancer cell lines
  (same three panels, E vs M genes)

SCLC gene signature source:
  Lissa et al. 2022 (Nat Comms) — 50-gene signature.
  25 NE genes (positively associated with neuroendocrine differentiation)
  25 non-NE genes (negatively associated).
  These are the exact genes used in Hari et al. Fig 1A.

EMT gene signature source:
  Aiello et al. 2018 — 22 genes split into Epithelial (11) and Mesenchymal (11).

Data:
  - OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv  (DepMap)
      Rows: one row per sample (with ModelID and IsDefaultEntryForModel columns)
      Cols: gene expression columns named "SYMBOL (EntrezID)"
  - Model.csv (DepMap)
      Contains ModelID, DepmapModelType, OncotreePrimaryDisease columns

Usage:
  Place all files in the same directory as this script, then run:
      python reproduce_fig1AB_final.py
"""

# ─────────────────────────────────────────────────────────────────────────────
# 0.  IMPORTS
# ─────────────────────────────────────────────────────────────────────────────
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

BASE = pathlib.Path(__file__).resolve().parent

# ─────────────────────────────────────────────────────────────────────────────
# 1.  GENE SIGNATURES
#
# IMPORTANT:
# - To match the paper figure, we must use the same gene lists (and ordering)
#   used in the analysis. In this project, those lists are provided in
#   Genesets.csv.
# - For SCLC: Genesets.csv stores the 50-gene signature in the order
#   [first 25 = NE] + [last 25 = non-NE].
# - For EMT: Genesets.csv stores the 22-gene signature (we split into
#   Epithelial (11) and Mesenchymal (11) using the standard Aiello list).
# ─────────────────────────────────────────────────────────────────────────────

GENESETS_FILE = BASE / "Genesets.csv"

# EMT signature split (Aiello et al. 2018) — used for ordering/splitting.
EMT_E = [
    "CDH1", "CLDN2", "CLDN4", "EPCAM", "ESRP1",
    "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "MUC1",
]
EMT_M = [
    "CDH11", "DDR2", "FAP", "GLI1", "LUM",
    "PALLD", "PDGFRB", "PDPN", "POSTN", "SPARC", "TNC",
]


def load_genesets(filepath: pathlib.Path) -> tuple[list[str], list[str], list[str], list[str]]:
    """
    Load gene lists from Genesets.csv.
    Returns: (SCLC_NE, SCLC_NONNE, EMT_E, EMT_M)
    """
    gs = pd.read_csv(filepath)
    if not {"GeneSetID", "Genes"}.issubset(set(gs.columns)):
        raise ValueError(f"{filepath.name} must have columns: GeneSetID, Genes")

    sclc = gs.loc[gs["GeneSetID"] == "SCLC", "Genes"].tolist()
    if len(sclc) != 50:
        raise ValueError(f"Expected 50 SCLC genes in {filepath.name}, found {len(sclc)}")
    sclc_ne, sclc_nonne = sclc[:25], sclc[25:]

    emt = gs.loc[gs["GeneSetID"] == "EMT", "Genes"].tolist()
    if len(emt) != 22:
        raise ValueError(f"Expected 22 EMT genes in {filepath.name}, found {len(emt)}")
    # sanity: EMT_E + EMT_M should match EMT set (order in file may differ)
    if set(emt) != (set(EMT_E) | set(EMT_M)):
        missing = sorted((set(EMT_E) | set(EMT_M)) - set(emt))
        extra = sorted(set(emt) - (set(EMT_E) | set(EMT_M)))
        raise ValueError(
            "EMT genes in Genesets.csv don't match expected Aiello EMT list.\n"
            f"Missing: {missing}\nExtra: {extra}"
        )

    return sclc_ne, sclc_nonne, list(EMT_E), list(EMT_M)

# ─────────────────────────────────────────────────────────────────────────────
# 2.  DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_ccle_expression(filepath: pathlib.Path) -> tuple[pd.DataFrame, dict]:
    """
    Load DepMap OmicsExpression file.
    - Keep only one row per model (IsDefaultEntryForModel == Yes).
    - Build a gene-symbol → full-column-name mapping.
    Returns:
        expr : DataFrame  (index = ModelID, columns = full "SYMBOL (ID)" names)
        gmap : dict       symbol → full column name
    """
    print(f"Loading expression data from {filepath.name} ...")
    expr = pd.read_csv(filepath, low_memory=False)

    # Keep default entry per model
    if "IsDefaultEntryForModel" in expr.columns:
        expr = expr[expr["IsDefaultEntryForModel"] == "Yes"].copy()

    # Set ModelID as index
    expr = expr.set_index("ModelID")

    # Drop non-expression metadata columns
    drop_cols = [c for c in expr.columns
                 if not (" (" in c and c.endswith(")"))]
    expr = expr.drop(columns=drop_cols, errors="ignore")

    # Build symbol → column name map
    gmap = {}
    for col in expr.columns:
        sym = col.split(" (")[0]
        gmap[sym] = col

    print(f"  Expression matrix: {expr.shape[0]} samples × {expr.shape[1]} genes")
    return expr, gmap


def load_model_info(filepath: pathlib.Path) -> pd.DataFrame:
    """Load Model.csv and return with ModelID as index."""
    print(f"Loading model info from {filepath.name} ...")
    meta = pd.read_csv(filepath, index_col="ModelID")
    print(f"  {len(meta)} models in metadata")
    return meta


def get_model_ids(meta: pd.DataFrame, group: str) -> set:
    """
    Return ModelID set for 'SCLC' or 'non_SCLC_cancer' groups.
    Uses DepmapModelType and OncotreePrimaryDisease columns.
    """
    if group == "SCLC":
        ids = meta[meta["DepmapModelType"] == "SCLC"].index
        print(f"  SCLC cell lines found: {len(ids)}")
        return set(ids)

    elif group == "non_SCLC_cancer":
        # All cancer lines that are NOT SCLC
        is_cancer = meta["OncotreePrimaryDisease"] != "Non-Cancerous"
        is_not_sclc = meta["DepmapModelType"] != "SCLC"
        ids = meta[is_cancer & is_not_sclc].index
        print(f"  Non-SCLC cancer cell lines found: {len(ids)}")
        return set(ids)

    else:
        raise ValueError(f"Unknown group: {group}")


# ─────────────────────────────────────────────────────────────────────────────
# 3.  MATRIX BUILDER
# ─────────────────────────────────────────────────────────────────────────────

def build_matrix(
    expr: pd.DataFrame,
    gmap: dict,
    model_ids: set,
    symbols: list
) -> tuple[pd.DataFrame, list, list]:
    """
    Subset expression matrix to:
      - rows  : model_ids that are present in expr
      - columns: genes in symbols that are present in gmap

    Returns:
        mat     : DataFrame (samples × resolved genes), numeric
        resolved: list of gene symbols actually found
        missing : list of gene symbols not found in expression data
    """
    resolved, cols, missing = [], [], []
    
    ALIASES = {
  "CYR61": "CCN1",
  "FAM57B": "CFAP97",
    }

    for sym in symbols:
        if sym in gmap:
            resolved.append(sym)
            cols.append(gmap[sym])
        elif sym in ALIASES and ALIASES[sym] in gmap:
            resolved.append(ALIASES[sym])
            cols.append(gmap[ALIASES[sym]])
        else:
            missing.append(sym)

    if missing:
        print(f"  Warning – {len(missing)} genes not found: {missing}")

    # Filter to the requested cell lines
    rows = [m for m in model_ids if m in expr.index]
    sub = expr.loc[rows, cols].apply(pd.to_numeric, errors="coerce")
    sub.columns = resolved          # rename to clean symbols
    sub = sub.dropna(how="all")     # drop rows that are all NaN

    print(f"  Matrix built: {sub.shape[0]} samples × {sub.shape[1]} genes")
    return sub, resolved, missing


# ─────────────────────────────────────────────────────────────────────────────
# 4.  PHENOTYPIC SCORE
# ─────────────────────────────────────────────────────────────────────────────

def phenotypic_score(mat: pd.DataFrame,
                     team1: list, team2: list) -> np.ndarray:
    """
    score = mean(team1 expression) − mean(team2 expression)  per sample.
    Positive → team1-like.  Negative → team2-like.
    """
    t1 = [g for g in team1 if g in mat.columns]
    t2 = [g for g in team2 if g in mat.columns]
    return mat[t1].mean(axis=1).values - mat[t2].mean(axis=1).values


# ─────────────────────────────────────────────────────────────────────────────
# 5.  SINGLE PANEL DRAWING
# ─────────────────────────────────────────────────────────────────────────────

def draw_panel(
    ax_hm:  plt.Axes,
    ax_bar: plt.Axes,
    ax_sc:  plt.Axes,
    mat:    pd.DataFrame,
    team1:  list,           # e.g. NE genes  / Epithelial genes
    team2:  list,           # e.g. non-NE    / Mesenchymal genes
    team1_label: str,       # axis bracket label
    team2_label: str,
    panel_title: str,
    colorbar_label: str,
) -> float:
    """
    Draw sub-panels (i), (ii), (iii) for one row (A or B).
    Returns the PC1 variance % for printing.
    """

    # ── filter to genes present in mat ───────────────────────────────────────
    t1 = [g for g in team1 if g in mat.columns]
    t2 = [g for g in team2 if g in mat.columns]
    ordered = t1 + t2          # team1 block first, then team2 block
    mat_ord = mat[ordered]

    # ── (i) Gene-gene Pearson correlation heatmap ────────────────────────────
    corr = mat_ord.corr(method="pearson")

    sns.heatmap(
        corr,
        ax=ax_hm,
        cmap="RdBu_r",
        center=0,
        vmin=-1, vmax=1,
        square=True,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Correlation", "shrink": 0.8},
        linewidths=0,
    )
    ax_hm.set_xticklabels(ax_hm.get_xticklabels(), fontsize=5, rotation=90)
    ax_hm.set_yticklabels(ax_hm.get_yticklabels(), fontsize=5, rotation=0)
    ax_hm.set_title(panel_title, fontsize=9, pad=4)

    # Draw bracket lines separating the two teams on the heatmap y-axis
    n1, n2 = len(t1), len(t2)
    total  = n1 + n2
    # team1 occupies top rows (0..n1-1), team2 occupies (n1..total-1)
    ax_hm.annotate(
        team1_label,
        xy=(-0.18, 1 - (n1 / total) / 2),
        xycoords="axes fraction",
        fontsize=7, rotation=90, va="center", ha="center",
    )
    ax_hm.annotate(
        team2_label,
        xy=(-0.18, (n2 / total) / 2),
        xycoords="axes fraction",
        fontsize=7, rotation=90, va="center", ha="center",
    )

    # ── PCA (used for both panels ii and iii) ────────────────────────────────
   # Center-only (like R prcomp default: center=TRUE, scale.=FALSE)
    X = mat_ord.values - mat_ord.values.mean(axis=0, keepdims=True)

    n_components = min(10, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(X)                   # shape: (n_samples, n_components)

    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    # ── (ii) PC1 loadings bar chart ──────────────────────────────────────────
    loadings = pd.Series(pca.components_[0], index=ordered)

    # Assign colors by team membership
    bar_colors = ["#c44e52" if g in t1 else "#4c72b0" for g in loadings.index]

    # Keep the gene order aligned with the heatmap (team1 block then team2 block).
    ax_bar.barh(
        range(len(loadings)),
        loadings.values,
        color=bar_colors,
        height=0.85,
        edgecolor="none",
    )
    ax_bar.set_yticks(range(len(loadings)))
    ax_bar.set_yticklabels(loadings.index, fontsize=6)
    ax_bar.axvline(0, color="gray", linewidth=0.6)
    ax_bar.set_xlabel("PC1 loading", fontsize=8)
    ax_bar.set_title(
        f"PC1 loadings ({pc1_var:.2f}% var.)",
        fontsize=9, pad=4,
    )
    ax_bar.invert_yaxis()

    # Add team bracket labels on right side of bar chart
    n_t2_sorted = len(t2)
    n_t1_sorted = len(t1)
    ax_bar.annotate(
        team2_label,
        xy=(1.02, (n_t2_sorted / 2) / len(loadings)),
        xycoords="axes fraction",
        fontsize=7, rotation=270, va="center",
    )
    ax_bar.annotate(
        team1_label,
        xy=(1.02, (n_t2_sorted + n_t1_sorted / 2) / len(loadings)),
        xycoords="axes fraction",
        fontsize=7, rotation=270, va="center",
    )

    # ── (iii) PCA scatter plot ───────────────────────────────────────────────
    pheno = phenotypic_score(mat_ord, t1, t2)

    sc = ax_sc.scatter(
        scores[:, 0], scores[:, 1],
        c=pheno,
        cmap="RdBu_r",
        s=20,
        alpha=0.85,
        edgecolors="none",
    )
    plt.colorbar(sc, ax=ax_sc, label=colorbar_label, shrink=0.8)
    ax_sc.set_xlabel(f"PC1 ({pc1_var:.2f}%)", fontsize=9)
    ax_sc.set_ylabel(f"PC2 ({pc2_var:.2f}%)", fontsize=9)
    ax_sc.tick_params(labelsize=7)
    ax_sc.set_title("PC1 vs PC2", fontsize=9, pad=4)

    return pc1_var


# ─────────────────────────────────────────────────────────────────────────────
# 6.  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():

    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"Loading gene sets from {GENESETS_FILE.name} ...")
    sclc_ne, sclc_nonne, emt_e, emt_m = load_genesets(GENESETS_FILE)

    expr_path  = BASE / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    model_path = BASE / "Model.csv"

    if not expr_path.exists():
        raise FileNotFoundError(
            f"Expression file not found: {expr_path}\n"
            "Download from https://depmap.org/portal/download/"
        )

    expr, gmap = load_ccle_expression(expr_path)
    meta       = load_model_info(model_path)

    # ── Get cell line subsets ─────────────────────────────────────────────────
    print("\n--- Identifying cell line subsets ---")
    sclc_ids     = get_model_ids(meta, "SCLC")
    non_sclc_ids = get_model_ids(meta, "non_SCLC_cancer")

    # ── Build expression matrices ─────────────────────────────────────────────
    print("\n--- Building expression matrices ---")
    all_sclc_genes = sclc_ne + sclc_nonne
    all_emt_genes  = emt_e + emt_m

    print("SCLC panel:")
    mat_sclc, _, _ = build_matrix(expr, gmap, sclc_ids, all_sclc_genes)

    print("EMT panel:")
    mat_emt, _, _  = build_matrix(expr, gmap, non_sclc_ids, all_emt_genes)

    # ── Build figure ──────────────────────────────────────────────────────────
    print("\n--- Drawing figure ---")

    fig = plt.figure(figsize=(16, 11))
    fig.suptitle(
        "Hari et al. 2025 — Figure 1A-B reproduction\n"
        "(SCLC signature from Lissa et al. 2022 · EMT signature from Aiello et al. 2018)",
        fontsize=10, y=0.99,
    )

    # 2-row × 3-column grid
    gs = gridspec.GridSpec(
        2, 3, figure=fig,
        hspace=0.5, wspace=0.42,
        left=0.06, right=0.97,
        top=0.94, bottom=0.06,
    )

    # Panel A axes (row 0)
    ax_A_hm  = fig.add_subplot(gs[0, 0])
    ax_A_bar = fig.add_subplot(gs[0, 1])
    ax_A_sc  = fig.add_subplot(gs[0, 2])

    # Panel B axes (row 1)
    ax_B_hm  = fig.add_subplot(gs[1, 0])
    ax_B_bar = fig.add_subplot(gs[1, 1])
    ax_B_sc  = fig.add_subplot(gs[1, 2])

    # Large panel labels
    fig.text(0.01, 0.96, "A", fontsize=14, fontweight="bold")
    fig.text(0.01, 0.48, "B", fontsize=14, fontweight="bold")

    # ── Draw Panel A ──────────────────────────────────────────────────────────
    print("Drawing Panel A (SCLC) ...")
    pc1_A = draw_panel(
        ax_hm        = ax_A_hm,
        ax_bar       = ax_A_bar,
        ax_sc        = ax_A_sc,
        mat          = mat_sclc,
        team1        = sclc_ne,
        team2        = sclc_nonne,
        team1_label  = "Neuroendocrine",
        team2_label  = "Non-Neuroendocrine",
        panel_title  = "CCLE SCLC cell lines — SCLC signature",
        colorbar_label = "Phenotypic Score\nmean(NE) − mean(non-NE)",
    )
    print(f"  Panel A PC1 variance: {pc1_A:.2f}%  (paper reports 56.97%)")

    # ── Draw Panel B ──────────────────────────────────────────────────────────
    print("Drawing Panel B (EMT) ...")
    pc1_B = draw_panel(
        ax_hm        = ax_B_hm,
        ax_bar       = ax_B_bar,
        ax_sc        = ax_B_sc,
        mat          = mat_emt,
        team1        = emt_e,
        team2        = emt_m,
        team1_label  = "Epithelial (E)",
        team2_label  = "Mesenchymal (M)",
        panel_title  = "CCLE non-SCLC cell lines — EMT signature",
        colorbar_label = "Phenotypic Score\nmean(E) − mean(M)",
    )
    print(f"  Panel B PC1 variance: {pc1_B:.2f}%  (paper reports 46.33%)")

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = BASE / "figure1AB_reproduction_final.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"\nSaved → {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
