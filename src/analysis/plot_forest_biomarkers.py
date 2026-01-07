from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"
FIG_DIR = REPORTS_DIR / "figures"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--biomarkers", type=str, default=str(TABLES_DIR / "biomarkers_survival_table.csv"))
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--out", type=str, default=str(FIG_DIR / "fig4_forest_cox_biomarkers.png"))
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.biomarkers)
    # Ensure columns exist
    need = {"gene", "hazard_ratio", "coef"}
    missing = need - set(df.columns)
    if missing:
        raise RuntimeError(f"biomarkers table missing columns: {sorted(missing)}")

    # strongest effects by |log(HR)|
    df = df.copy()
    df["hazard_ratio"] = pd.to_numeric(df["hazard_ratio"], errors="coerce")
    df = df.dropna(subset=["hazard_ratio"])
    df["abs_log_hr"] = np.abs(np.log(df["hazard_ratio"]))
    df = df.sort_values("abs_log_hr", ascending=False).head(args.top_n)

    # Plot (horizontal)
    genes = df["gene"].astype(str).tolist()
    hr = df["hazard_ratio"].values
    y = np.arange(len(genes))

    plt.figure(figsize=(8.5, max(5.0, 0.35 * len(genes) + 1.5)))
    plt.axvline(1.0, linestyle="--")

    plt.scatter(hr, y)
    plt.yticks(y, genes)
    plt.gca().invert_yaxis()

    plt.xscale("log")
    plt.xlabel("Hazard Ratio (log scale)")
    plt.title(f"Top Cox Biomarkers (n={len(genes)})")

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    LOG.info("✅ Saved: %s", args.out)


if __name__ == "__main__":
    main()
