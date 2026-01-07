from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
FIG_DIR = REPORTS_DIR / "figures"


def _load_train() -> tuple[pd.DataFrame, pd.Series]:
    X = pd.read_parquet(PROCESSED_DIR / "X_train.parquet")
    y = pd.read_csv(PROCESSED_DIR / "y_class_train.csv", index_col="sample_id")["class_label"]
    X = X.loc[y.index]
    return X, y


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--method", choices=["pca", "umap"], default="pca", help="UMAP requires umap-learn installed.")
    p.add_argument("--n-components", type=int, default=2)
    p.add_argument("--max-features", type=int, default=2000, help="Use top-variance features for speed.")
    p.add_argument("--out", type=str, default=str(FIG_DIR / "fig1_embedding_subtypes.png"))
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    X, y = _load_train()
    LOG.info("Loaded X=%s, y=%s", X.shape, y.shape)

    # feature downselect by variance
    if args.max_features and args.max_features < X.shape[1]:
        var = X.var(axis=0).sort_values(ascending=False)
        keep = var.head(args.max_features).index
        X = X[keep]
        LOG.info("Using top-variance features: %d", X.shape[1])

    Xv = X.values.astype(float)
    Xv = StandardScaler(with_mean=True, with_std=True).fit_transform(Xv)

    if args.method == "pca":
        emb = PCA(n_components=args.n_components, random_state=0).fit_transform(Xv)
        title = "PCA of RNA-seq (colored by histological subtype)"
        xlab, ylab = "PC1", "PC2"
    else:
        try:
            import umap  # type: ignore
        except Exception as e:
            raise RuntimeError("UMAP selected but umap-learn not installed. pip install umap-learn") from e
        emb = umap.UMAP(n_components=args.n_components, random_state=0).fit_transform(Xv)
        title = "UMAP of RNA-seq (colored by histological subtype)"
        xlab, ylab = "UMAP1", "UMAP2"

    classes = sorted(y.dropna().unique().tolist())
    class_to_int = {c: i for i, c in enumerate(classes)}
    colors = [class_to_int.get(v, -1) for v in y.tolist()]

    plt.figure(figsize=(7.5, 6))
    sc = plt.scatter(emb[:, 0], emb[:, 1], c=colors, s=18, alpha=0.8)

    # legend
    handles = []
    for c in classes:
        handles.append(plt.Line2D([0], [0], marker="o", linestyle="", label=c))
    plt.legend(handles=handles, title="Subtype", loc="best", frameon=True)

    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    LOG.info("✅ Saved: %s", args.out)


if __name__ == "__main__":
    main()
