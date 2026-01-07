from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
FIG_DIR = REPORTS_DIR / "figures"
TABLES_DIR = REPORTS_DIR / "tables"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--a", type=str, default="Oligodendroglioma", help="Class A (case)")
    p.add_argument("--b", type=str, default="Astrocytoma", help="Class B (control)")
    p.add_argument("--max-genes", type=int, default=8000, help="Top-variance genes to test (speed).")
    p.add_argument("--out", type=str, default=str(FIG_DIR / "fig6_volcano_oligo_vs_astro.png"))
    p.add_argument("--out-table", type=str, default=str(TABLES_DIR / "de_oligo_vs_astro.csv"))
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    X = pd.read_parquet(PROCESSED_DIR / "X_train.parquet")
    y = pd.read_csv(PROCESSED_DIR / "y_class_train.csv", index_col="sample_id")["class_label"]
    X = X.loc[y.index]

    mask_a = (y == args.a)
    mask_b = (y == args.b)

    Xa = X.loc[mask_a]
    Xb = X.loc[mask_b]

    if Xa.empty or Xb.empty:
        raise RuntimeError(f"Not enough samples for A={args.a} or B={args.b}")

    # variance filter for speed
    if args.max_genes and args.max_genes < X.shape[1]:
        keep = X.var(axis=0).sort_values(ascending=False).head(args.max_genes).index
        Xa = Xa[keep]
        Xb = Xb[keep]

    # log2 fold-change (mean difference)
    mean_a = Xa.mean(axis=0)
    mean_b = Xb.mean(axis=0)
    log2fc = np.log2((mean_a + 1e-6) / (mean_b + 1e-6))

    # t-test per gene
    pvals = []
    for col in Xa.columns:
        stat, pv = ttest_ind(Xa[col].values, Xb[col].values, equal_var=False, nan_policy="omit")
        pvals.append(pv)
    pvals = np.array(pvals, dtype=float)

    # FDR correction
    _, qvals, _, _ = multipletests(pvals, method="fdr_bh")

    de = pd.DataFrame({
        "gene": Xa.columns.astype(str),
        "log2FC": log2fc.values,
        "p_value": pvals,
        "q_value": qvals,
    }).sort_values("q_value")

    de.to_csv(args.out_table, index=False)
    LOG.info("✅ Saved DE table: %s", args.out_table)

    # Volcano plot
    x = de["log2FC"].values
    yv = -np.log10(de["q_value"].values + 1e-300)

    plt.figure(figsize=(8, 6))
    plt.scatter(x, yv, s=10, alpha=0.7)

    # threshold lines
    plt.axhline(-np.log10(0.05), linestyle="--")
    plt.axvline(1.0, linestyle="--")
    plt.axvline(-1.0, linestyle="--")

    plt.title(f"Volcano: {args.a} vs {args.b} (FDR BH)")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(q-value)")

    # label top 10 genes by q
    top = de.head(10)
    for _, r in top.iterrows():
        plt.text(r["log2FC"], -np.log10(r["q_value"] + 1e-300), str(r["gene"]), fontsize=8)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    LOG.info("✅ Saved: %s", args.out)


if __name__ == "__main__":
    main()
