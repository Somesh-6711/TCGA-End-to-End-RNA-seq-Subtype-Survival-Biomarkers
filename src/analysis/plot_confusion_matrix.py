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
    p.add_argument("--cm-csv", type=str, default=str(TABLES_DIR / "classification_confusion_matrix.csv"))
    p.add_argument("--out", type=str, default=str(FIG_DIR / "fig2_confusion_matrix.png"))
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    cm_df = pd.read_csv(args.cm_csv, index_col=0)
    cm = cm_df.values
    labels = cm_df.index.tolist()

    plt.figure(figsize=(7, 6))
    plt.imshow(cm, interpolation="nearest")
    plt.title("Confusion Matrix (Test)")
    plt.colorbar()

    tick_marks = np.arange(len(labels))
    plt.xticks(tick_marks, labels, rotation=35, ha="right")
    plt.yticks(tick_marks, labels)

    # annotate
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, str(cm[i, j]), ha="center", va="center")

    plt.ylabel("True label")
    plt.xlabel("Predicted label")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    LOG.info("✅ Saved: %s", args.out)


if __name__ == "__main__":
    main()
