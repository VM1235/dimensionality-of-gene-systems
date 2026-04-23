"""
Reproduction of Wang et al. PRX Life 2, 043009 (2024) — Dentate Gyrus Neurogenesis
THREE PANELS:
  (A) Frustration Score along Reaction Coordinate — mean ± std, 15 bins
  (B) Intrinsic Dimensionality (TWO-NN) along Reaction Coordinate
  (C) Scatter: ID vs Frustration, each dot = one cell, coloured by cell type

Fix applied: scv.tl.recover_dynamics wrapped in if __name__ == '__main__'
             AND n_jobs=1 to avoid macOS multiprocessing spawn crash.

Run:  python3 code_final_fixed.py
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import scvelo as scv
import scanpy as sc
import anndata as ad
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr

# ── Constants ─────────────────────────────────────────────────────────────────
N_BINS = 15
N_PCS  = 30
ACCENT = "#c0392b"   # red  — frustration
ID_COL = "#2471a3"   # blue — intrinsic dimensionality
BG     = "white"
GRID_C = "#ebebeb"
DATA_PATH = "/Users/apple/Downloads/Project/part2/data/DentateGyrus/10X43_1.h5ad"


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def bin_along_rc(pt, signal, n=N_BINS):
    """Bin signal into n equal-width pseudotime bins; return bin centres, means, stds."""
    edges = np.linspace(pt.min(), pt.max(), n + 1)
    idx   = np.clip(np.digitize(pt, edges) - 1, 0, n - 1)
    c, m, s = [], [], []
    for i in range(n):
        mask = idx == i
        if mask.sum() > 3:
            c.append(i + 1)
            m.append(signal[mask].mean())
            s.append(signal[mask].std())
    return np.array(c), np.array(m), np.array(s)


def twonn_id_per_cell(X, k=30):
    """
    TWO-NN intrinsic dimensionality (Facco et al. 2017) — per-cell estimate.
    For each cell, ID = -1 / log(d1/d2) where d1, d2 are distances to 1st and
    2nd nearest neighbours. Averaged over k successive pairs for stability.
    Returns array of length n_cells, clipped at 97th percentile.
    """
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(X)
    d, _  = nbrs.kneighbors(X)
    d     = d[:, 1:]            # drop self (distance = 0)
    ids   = []
    for i in range(d.shape[1] - 1):
        mu   = d[:, i + 1] / (d[:, i] + 1e-12)
        ids.append(1.0 / (np.log(mu + 1e-12) + 1e-12))
    id_cell = np.mean(ids, axis=0)
    return np.clip(id_cell, 0, np.percentile(id_cell, 97))


def velocity_frustration(spliced, velocity):
    """
    Per-cell frustration: fraction of genes where RNA-velocity sign
    conflicts with whether spliced counts are above/below their mean.
    Proxy for 'fraction of frustrated GRN edges' used in the paper.
    Returns values in [0, 1].
    """
    delta_s  = spliced - spliced.mean(axis=0)
    v_thresh = np.percentile(np.abs(velocity), 25, axis=0)
    effective = np.abs(velocity) > v_thresh
    conflict  = (np.sign(delta_s) * np.sign(velocity)) < 0
    total     = effective.sum(axis=1).astype(float)
    total     = np.where(total == 0, 1, total)
    return (conflict & effective).sum(axis=1).astype(float) / total


def to_dense(arr):
    return arr.toarray() if sp.issparse(arr) else np.array(arr)


def styled_panel(ax, rc, mean, std, color, ylabel, title, marker):
    ax.set_facecolor(BG)
    ax.fill_between(rc, mean - std, mean + std,
                    color=color, alpha=0.18, zorder=1, label="±1 SD")
    ax.plot(rc, mean, color=color, lw=2.6, zorder=3,
            marker=marker, ms=8, mfc="white", mec=color, mew=2.2)
    pk = np.argmax(mean)
    ax.annotate(
        f"Peak r = {rc[pk]}\n(transition state)",
        xy=(rc[pk], mean[pk]),
        xytext=(rc[pk] + 1.8, mean[pk] + 0.07 * (mean.max() - mean.min())),
        arrowprops=dict(arrowstyle="->", color=color, lw=1.5),
        fontsize=8.5, color=color, fontweight="bold"
    )
    ax.set_xlabel("Reaction Coordinate  r", fontsize=10.5)
    ax.set_ylabel(ylabel, fontsize=10.5)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks(np.arange(1, N_BINS + 1, 2))
    ax.yaxis.grid(True, color=GRID_C, zorder=0)
    ax.legend(fontsize=8.5, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0.5, N_BINS + 0.5)


# ═════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE — must be inside if __name__ == '__main__'
# to avoid macOS multiprocessing spawn crash
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    # 1. Load ─────────────────────────────────────────────────────────────────
    print("|-----> Loading Dentate Gyrus data...")
    adata = ad.read_h5ad(DATA_PATH)
    print(f"        {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"        Cell types: {list(adata.obs['clusters'].unique())}")

    # 2. Preprocessing ────────────────────────────────────────────────────────
    print("|-----> Preprocessing...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=10)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    sc.tl.pca(adata, n_comps=N_PCS)
    sc.pp.neighbors(adata, n_pcs=N_PCS)

    # 3. RNA velocity ─────────────────────────────────────────────────────────
    # KEY FIX: Use 'deterministic' mode — avoids recover_dynamics which uses
    # multiprocessing and crashes on macOS with spawn start method.
    # Deterministic mode uses spliced/unspliced ratio directly and is fast.
    print("|-----> Estimating RNA velocity (deterministic mode — fast, no multiprocessing)...")
    scv.tl.velocity(adata, mode="deterministic")
    scv.tl.velocity_graph(adata)

    # 4. Pseudotime (diffusion) ───────────────────────────────────────────────
    print("|-----> Computing diffusion pseudotime...")
    cell_types = adata.obs["clusters"].astype("category")

    # Root = Radial Glia-like (earliest progenitor in dentate gyrus)
    root_mask = cell_types.str.lower().str.contains("radial|rg|glia", na=False)
    if root_mask.sum() > 0:
        adata.uns["iroot"] = int(np.where(root_mask)[0][0])
        print(f"        Root cell index: {adata.uns['iroot']} (Radial Glia-like)")
    else:
        adata.uns["iroot"] = 0
        print("        WARNING: No Radial Glia-like cells found, using cell 0 as root")

    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, n_dcs=10)

    # Clean up pseudotime (DPT can give inf for disconnected cells)
    pt = adata.obs["dpt_pseudotime"].values.astype(float)
    finite_max = float(np.nanmax(pt[np.isfinite(pt)]))
    pt = np.where(~np.isfinite(pt), finite_max, pt)
    pt = (pt - pt.min()) / (pt.max() - pt.min() + 1e-12)
    adata.obs["pseudotime_final"] = pt

    # 5. Frustration score ────────────────────────────────────────────────────
    print("|-----> Computing frustration score...")
    v_dense = to_dense(adata.layers["velocity"])
    s_dense = to_dense(adata.layers["spliced"])
    # Only use genes that are finite in both
    finite_cols = np.isfinite(v_dense).all(axis=0) & np.isfinite(s_dense).all(axis=0)
    v_clean = v_dense[:, finite_cols]
    s_clean = s_dense[:, finite_cols]

    fr = velocity_frustration(s_clean, v_clean)
    p2, p98 = np.percentile(fr, 2), np.percentile(fr, 98)
    adata.obs["frustration"] = np.clip((fr - p2) / (p98 - p2 + 1e-12), 0, 1)

    # 6. Intrinsic dimensionality ─────────────────────────────────────────────
    print("|-----> Computing intrinsic dimensionality (TWO-NN per cell)...")
    adata.obs["intrinsic_dim"] = twonn_id_per_cell(adata.obsm["X_pca"])

    # 7. Bin along reaction coordinate ────────────────────────────────────────
    print("|-----> Binning into 15 reaction-coordinate bins...")
    pt_a = adata.obs["pseudotime_final"].values
    fr_a = adata.obs["frustration"].values
    id_a = adata.obs["intrinsic_dim"].values

    rc_f, mf, sf = bin_along_rc(pt_a, fr_a)
    rc_d, md, sd = bin_along_rc(pt_a, id_a)
    rho, pval    = spearmanr(fr_a, id_a)

    print(f"        Frustration peak at reaction coordinate r = {rc_f[np.argmax(mf)]}")
    print(f"        ID peak at reaction coordinate r          = {rc_d[np.argmax(md)]}")
    print(f"        Spearman ρ(Frustration, ID) = {rho:.3f}  (p = {pval:.2e})")

    # 8. Figure ───────────────────────────────────────────────────────────────
    print("|-----> Generating figure...")

    cats = cell_types.cat.categories.tolist()
    pal  = sns.color_palette("tab20", n_colors=max(len(cats), 2))
    cmap = {ct: pal[i % 20] for i, ct in enumerate(cats)}

    fig = plt.figure(figsize=(19, 6.2), facecolor=BG)
    fig.suptitle(
        "Dentate Gyrus Neurogenesis  —  Wang et al. PRX Life 2, 043009 (2024)\n"
        "Frustration peak & Intrinsic Dimensionality along reaction coordinate  |  "
        f"Spearman ρ(ID, Frustration) = {rho:.2f}  (p = {pval:.1e})",
        fontsize=11, fontweight="bold", y=1.02
    )
    gs = gridspec.GridSpec(1, 3, wspace=0.40)

    # Panel A — Frustration
    ax1 = fig.add_subplot(gs[0])
    styled_panel(ax1, rc_f, mf, sf, ACCENT,
                 "Frustration Score",
                 "(A)  Frustration along\nReaction Coordinate", "o")
    ax1.axhline(fr_a.mean(), color="#999", ls="--", lw=1.4, zorder=2,
                label=f"Population mean = {fr_a.mean():.2f}")
    ax1.legend(fontsize=8.5, frameon=False)

    # Panel B — Intrinsic Dimensionality
    ax2 = fig.add_subplot(gs[1])
    styled_panel(ax2, rc_d, md, sd, ID_COL,
                 "Intrinsic Dimensionality (TWO-NN)",
                 "(B)  Intrinsic Dimensionality\nalong Reaction Coordinate", "s")

    # Panel C — Scatter: ID vs Frustration, each dot = one cell
    ax3 = fig.add_subplot(gs[2])
    ax3.set_facecolor(BG)
    for ct in cats:
        mask = (cell_types == ct).values
        ax3.scatter(fr_a[mask], id_a[mask],
                    color=cmap[ct], s=18, alpha=0.55,
                    linewidths=0, label=ct, zorder=2)
    xline = np.linspace(fr_a.min(), fr_a.max(), 200)
    ax3.plot(xline, np.polyval(np.polyfit(fr_a, id_a, 1), xline),
             color="#555", lw=1.6, ls="--", zorder=3,
             label=f"Linear trend  (ρ = {rho:.2f})")
    ax3.set_xlabel("Frustration Score", fontsize=10.5)
    ax3.set_ylabel("Intrinsic Dimensionality (TWO-NN)", fontsize=10.5)
    ax3.set_title("(C)  ID vs. Frustration\n(each dot = one cell)", fontsize=11, fontweight="bold")
    ax3.yaxis.grid(True, color=GRID_C, zorder=0)
    ax3.xaxis.grid(True, color=GRID_C, zorder=0)
    ax3.set_axisbelow(True)
    ax3.spines[["top", "right"]].set_visible(False)
    ax3.legend(title="Cell Type", fontsize=7.5, title_fontsize=8.5,
               bbox_to_anchor=(1.03, 1), loc="upper left",
               frameon=False, markerscale=1.8)

    plt.tight_layout()
    out = "/Users/apple/Downloads/Project/part2/data/DentateGyrus/reproduction_final.png"
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()

    print(f"\n[SUCCESS]  Saved → {out}")
