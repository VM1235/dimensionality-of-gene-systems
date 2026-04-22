#!/usr/bin/env python3
"""
Compute Intrinsic Dimensionality (ID) across EMT time course for Cook et al. 2020
(Nature Communications, GSE147405).

Dataset: A549 cell line, TGFb1 treatment (and optionally all 12 combinations).

Methods (same as replicate_fig1ab.py in this repo):
- Normalization : counts / total UMIs per cell (relative abundances)
- Estimator     : TWO-NN (skdim.id.TwoNN) → one scalar ID per group of cells
- Subsampling   : floor(0.75 * min_group_size) cells, 10 draws → mean ± SD
- ID-score      : rescaled to [0, 1] within each panel (induction / reversal)

Time points:
  Induction : 0d → 8h → 1d → 3d → 7d
  Reversal  : 7d → 8h_rm → 1d_rm → 3d_rm   (_rm = stimulus removed)

Usage:
    python emt_id_cook2020.py \
        --matrix  data/GSE147405_A549_TGFB1_TimeCourse_UMI_matrix.csv \
        --meta    data/GSE147405_A549_TGFB1_TimeCourse_metadata.csv \
        --label   A549_TGFB1
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from skdim.id import TwoNN

# ── time-point ordering ───────────────────────────────────────────────────────
INDUCTION_ORDER = ["0d", "8h", "1d", "3d", "7d"]
INDUCTION_NUMERIC = [0, 0.33, 1, 3, 7]          # x-axis in days (8h ≈ 0.33 d)

REVERSAL_ORDER = ["7d", "8h_rm", "1d_rm", "3d_rm"]
REVERSAL_NUMERIC = [0, 0.33, 1, 3]               # x-axis = days since removal


# ── core functions (identical to replicate_fig1ab.py) ────────────────────────

def normalize_relative_abundance(X) -> np.ndarray:
    """Counts / total UMIs per cell → dense float32."""
    if sp.issparse(X):
        total = np.asarray(X.sum(axis=1)).ravel()
        total = np.maximum(total, 1e-12)
        Xn = X.multiply(1.0 / total[:, None])
        return np.asarray(Xn.todense(), dtype=np.float32)
    X = np.asarray(X, dtype=np.float64)
    total = X.sum(axis=1, keepdims=True)
    total = np.maximum(total, 1e-12)
    return (X / total).astype(np.float32)


def intrinsic_dimension_twonn(X: np.ndarray, seed: int) -> float:
    """TWO-NN global ID — single scalar for the point cloud."""
    rng = np.random.default_rng(seed)
    idx = rng.permutation(X.shape[0])
    Xs = np.ascontiguousarray(X[idx])
    return float(TwoNN().fit_transform(Xs))


def rescale_unit(y: np.ndarray) -> np.ndarray:
    lo, hi = float(np.min(y)), float(np.max(y))
    if hi - lo < 1e-15:
        return np.zeros_like(y)
    return (y - lo) / (hi - lo)


# ── ID computation per group of time points ───────────────────────────────────

def compute_id_over_timepoints(
    expr: np.ndarray,           # shape (n_genes, n_cells) — raw UMI matrix
    cell_ids: np.ndarray,       # cell barcodes (columns of expr)
    meta: pd.DataFrame,         # metadata with 'Time' column, indexed by barcode
    time_order: list[str],
    n_draws: int = 10,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each time point in time_order, subsample cells, normalize, run TWO-NN.
    Returns (means, stds) arrays of length len(time_order).
    """
    # build dict: time_label → row indices into expr (expr is genes×cells, so cols)
    idx_per_tp: dict[str, np.ndarray] = {}
    counts = []
    for tp in time_order:
        cells_at_tp = meta.index[meta["Time"] == tp]
        # find positions of these cells in the expression matrix columns
        cell_pos = np.where(np.isin(cell_ids, cells_at_tp))[0]
        idx_per_tp[tp] = cell_pos
        counts.append(len(cell_pos))
        print(f"  {tp}: {len(cell_pos)} cells")

    n_sub = int(0.75 * min(counts))
    print(f"  → subsampling {n_sub} cells per time point ({n_draws} draws)\n")

    means, stds = [], []
    for tp in time_order:
        pool = idx_per_tp[tp]
        ids = []
        for k in range(n_draws):
            rng = np.random.default_rng(seed + k * 7919)
            pick = rng.choice(pool, size=n_sub, replace=False)
            # expr is (genes × cells): take columns, then transpose → (cells × genes)
            X_sub = expr[:, pick].T
            X_norm = normalize_relative_abundance(X_sub)
            id_val = intrinsic_dimension_twonn(X_norm, seed=seed + k)
            ids.append(id_val)
        ids = np.array(ids)
        means.append(ids.mean())
        stds.append(ids.std(ddof=0))
        print(f"  {tp}: raw ID = {ids.mean():.4f} ± {ids.std(ddof=0):.4f}")

    return np.array(means), np.array(stds)


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_id(
    times_numeric: list[float],
    means: np.ndarray,
    stds: np.ndarray,
    time_labels: list[str],
    title: str,
    xlabel: str,
    outfile: Path,
    color: str = "C0",
):
    score = rescale_unit(means)
    lo, hi = float(np.min(means)), float(np.max(means))
    scale = 1.0 / (hi - lo) if (hi - lo) > 1e-15 else 0.0
    err = stds * scale

    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=150)
    ax.errorbar(
        times_numeric, score, yerr=err,
        fmt="o-", capsize=4, color=color,
        ecolor="gray", markersize=7, linewidth=1.8,
    )
    ax.set_xticks(times_numeric)
    ax.set_xticklabels(time_labels, rotation=20, ha="right")
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("ID-score (rescaled to [0,1])")
    ax.set_title(title)
    fig.tight_layout()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {outfile}\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True,
                   help="Path to UMI matrix CSV (genes × cells)")
    p.add_argument("--meta", required=True,
                   help="Path to metadata CSV (cells × features, must have 'Time' column)")
    p.add_argument("--label", default="A549_TGFB1",
                   help="Label used in plot titles and output filenames")
    p.add_argument("--outdir", default="figures",
                   help="Output directory for PNG files")
    p.add_argument("--n-draws", type=int, default=10)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    outdir = Path(args.outdir)

    # ── load data ─────────────────────────────────────────────────────────────
    print("Loading metadata...")
    meta = pd.read_csv(args.meta, index_col=0)
    # keep only Singlets (doublets already removed in this dataset, but safe to filter)
    if "Doublet" in meta.columns:
        meta = meta[meta["Doublet"] == "Singlet"]

    print(f"Loading UMI matrix (this may take a moment)...")
    expr_df = pd.read_csv(args.matrix, index_col=0)
    # expr_df: rows = genes, cols = cells
    cell_ids = np.array(expr_df.columns)
    expr = expr_df.values  # (n_genes, n_cells)

    print(f"\nMatrix shape: {expr.shape[0]} genes × {expr.shape[1]} cells")
    print(f"Metadata cells: {len(meta)}")
    print(f"Time points present: {sorted(meta['Time'].unique())}\n")

    # ── induction ─────────────────────────────────────────────────────────────
    print("=== EMT INDUCTION ===")
    means_ind, stds_ind = compute_id_over_timepoints(
        expr, cell_ids, meta,
        time_order=INDUCTION_ORDER,
        n_draws=args.n_draws,
        seed=args.seed,
    )
    plot_id(
        times_numeric=INDUCTION_NUMERIC,
        means=means_ind,
        stds=stds_ind,
        time_labels=INDUCTION_ORDER,
        title=f"{args.label} — EMT Induction",
        xlabel="Time",
        outfile=outdir / f"{args.label}_induction_ID.png",
        color="C0",
    )

    # ── reversal ──────────────────────────────────────────────────────────────
    print("=== EMT REVERSAL ===")
    means_rev, stds_rev = compute_id_over_timepoints(
        expr, cell_ids, meta,
        time_order=REVERSAL_ORDER,
        n_draws=args.n_draws,
        seed=args.seed,
    )
    plot_id(
        times_numeric=REVERSAL_NUMERIC,
        means=means_rev,
        stds=stds_rev,
        time_labels=REVERSAL_ORDER,
        title=f"{args.label} — EMT Reversal",
        xlabel="Days since stimulus removed",
        outfile=outdir / f"{args.label}_reversal_ID.png",
        color="C2",
    )


if __name__ == "__main__":
    main()
