from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"
FIG_DIR = REPORTS_DIR / "figures"


def _read_indexed_csv(path: str, preferred_index: str = "sample_id") -> pd.DataFrame:
    df = pd.read_csv(path)
    if preferred_index in df.columns:
        df = df.set_index(preferred_index)
    else:
        first = df.columns[0]
        if str(first).lower().startswith("unnamed") or str(first).lower() in {"index", "sample"}:
            df = df.set_index(first)
        else:
            df = df.set_index(first)
    df.index = df.index.astype(str)
    return df


def _zscore_rows(M: np.ndarray) -> np.ndarray:
    mu = np.nanmean(M, axis=1, keepdims=True)
    sd = np.nanstd(M, axis=1, keepdims=True) + 1e-8
    return (M - mu) / sd


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--biomarkers", type=str, default=str(TABLES_DIR / "biomarkers_survival_table.csv"))
    p.add_argument("--groups", type=str, default=str(TABLES_DIR / "km_groups_test.csv"))
    p.add_argument("--X-test", type=str, default=str(PROCESSED_DIR / "X_test.parquet"))
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--out", type=str, default=str(FIG_DIR / "fig5_biomarker_heatmap_risk.png"))
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    bio = pd.read_csv(args.biomarkers)
    if "gene" not in bio.columns:
        raise RuntimeError(f"'gene' column not found in {args.biomarkers}. Found: {bio.columns.tolist()}")
    genes = bio["gene"].dropna().astype(str).head(args.top_n).tolist()

    X = pd.read_parquet(args.X_test)
    g = _read_indexed_csv(args.groups, "sample_id")

    if "group" not in g.columns and "risk_group" in g.columns:
        g = g.rename(columns={"risk_group": "group"})
    if "group" not in g.columns:
        raise RuntimeError(f"'group' column not found in {args.groups}. Found: {g.columns.tolist()}")

    # align
    common = X.index.astype(str).intersection(g.index.astype(str))
    X = X.loc[common]
    g = g.loc[common]

    g["group"] = g["group"].astype(str).str.strip().str.lower()

    # order samples
    if "risk_score" in g.columns:
        g = g.sort_values(["group", "risk_score"])
    else:
        g = g.sort_values(["group"])

    genes_present = [gg for gg in genes if gg in X.columns]
    if not genes_present:
        raise RuntimeError(
            "None of the biomarker genes were found in X_test columns. "
            "Check if expression columns are gene symbols vs Ensembl IDs."
        )

    M = X.loc[g.index, genes_present].T.values.astype(float)  # genes x samples
    Mz = _zscore_rows(M)

    plt.figure(figsize=(10, max(5.5, 0.22 * Mz.shape[0] + 2)))
    plt.imshow(Mz, aspect="auto", interpolation="nearest")
    plt.colorbar(label="Row z-score")

    plt.yticks(np.arange(len(genes_present)), genes_present)
    plt.xticks([])

    # group boundaries
    groups = g["group"].astype(str).tolist()
    boundaries = []
    for i in range(1, len(groups)):
        if groups[i] != groups[i - 1]:
            boundaries.append(i)
    for b in boundaries:
        plt.axvline(b - 0.5, linewidth=2)

    plt.title("Top Biomarker Expression Heatmap (Test) by Risk Group")
    plt.xlabel("Samples (ordered by risk group)")
    plt.ylabel("Genes")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    LOG.info("✅ Saved: %s", args.out)


if __name__ == "__main__":
    main()
