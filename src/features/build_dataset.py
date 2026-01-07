from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"


def _ensure_dirs() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)


def _patient_split(patient_ids: pd.Series, test_size: float, seed: int) -> Tuple[pd.Index, pd.Index]:
    """
    Returns (train_patients, test_patients) as Index of unique patient IDs.
    """
    uniq = patient_ids.dropna().astype(str).unique()
    rng = np.random.default_rng(seed)
    rng.shuffle(uniq)

    n_test = int(round(len(uniq) * test_size))
    n_test = max(1, min(n_test, len(uniq) - 1))

    test_patients = pd.Index(uniq[:n_test], name="patient_id")
    train_patients = pd.Index(uniq[n_test:], name="patient_id")
    return train_patients, test_patients


def _top_variance_genes(X: pd.DataFrame, k: int) -> pd.DataFrame:
    """
    Keep top-k genes by variance (computed on TRAIN only).
    """
    if k <= 0 or k >= X.shape[1]:
        return X
    v = X.var(axis=0, skipna=True)
    keep = v.sort_values(ascending=False).head(k).index
    return X.loc[:, keep]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-genes", type=int, default=5000)
    parser.add_argument("--log1p", action="store_true", help="Apply log1p to expression features.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    _ensure_dirs()

    expr_path = INTERIM_DIR / "expression.parquet"
    labels_path = INTERIM_DIR / "labels.parquet"
    if not expr_path.exists() or not labels_path.exists():
        raise FileNotFoundError("Missing interim files. Run: python -m src.data.download_xena --cohort LGG --force")

    X = pd.read_parquet(expr_path)
    labels = pd.read_parquet(labels_path)

    # Align
    common = X.index.intersection(labels.index)
    X = X.loc[common].copy()
    labels = labels.loc[common].copy()

    # Define task masks
    mask_class = labels["class_label"].notna()
    mask_surv = labels["os_time"].notna() & labels["os_event"].notna()

    # We will create ONE split that works for both tasks:
    # Use patients that have at least one labeled sample for either task (union),
    # but later we subset per task.
    eligible = labels.loc[mask_class | mask_surv, "patient_id"].dropna().astype(str)
    if eligible.empty:
        raise RuntimeError("No eligible patients found with class_label or survival labels.")

    train_p, test_p = _patient_split(eligible, test_size=args.test_size, seed=args.seed)

    train_idx = labels.index[labels["patient_id"].astype(str).isin(train_p)]
    test_idx = labels.index[labels["patient_id"].astype(str).isin(test_p)]

    # Train/test expression
    X_train = X.loc[train_idx]
    X_test = X.loc[test_idx]

    if args.log1p:
        # Expression is already roughly log scale in many Xena matrices; still optional.
        X_train = np.log1p(X_train)
        X_test = np.log1p(X_test)

    # Feature selection on TRAIN only
    X_train = _top_variance_genes(X_train, args.top_genes)
    X_test = X_test.loc[:, X_train.columns]

    # Task 2: classification labels
    y_train = labels.loc[train_idx, ["patient_id", "class_label"]].copy()
    y_test = labels.loc[test_idx, ["patient_id", "class_label"]].copy()

    # Task 1: survival labels
    s_train = labels.loc[train_idx, ["patient_id", "os_time", "os_event"]].copy()
    s_test = labels.loc[test_idx, ["patient_id", "os_time", "os_event"]].copy()

    # Save
    X_train.to_parquet(PROCESSED_DIR / "X_train.parquet")
    X_test.to_parquet(PROCESSED_DIR / "X_test.parquet")

    y_train.to_csv(PROCESSED_DIR / "y_class_train.csv", index_label="sample_id")
    y_test.to_csv(PROCESSED_DIR / "y_class_test.csv", index_label="sample_id")

    s_train.to_csv(PROCESSED_DIR / "survival_train.csv", index_label="sample_id")
    s_test.to_csv(PROCESSED_DIR / "survival_test.csv", index_label="sample_id")

    # Report
    def _counts(df: pd.DataFrame) -> int:
        return int(df.shape[0])

    split_summary = pd.DataFrame([{
        "n_samples_total": int(labels.shape[0]),
        "n_patients_total": int(labels["patient_id"].astype(str).nunique()),
        "n_patients_train": int(len(train_p)),
        "n_patients_test": int(len(test_p)),
        "n_genes_out": int(X_train.shape[1]),
        "train_samples_total": int(len(train_idx)),
        "test_samples_total": int(len(test_idx)),
        "train_class_labeled": int((y_train["class_label"].notna()).sum()),
        "test_class_labeled": int((y_test["class_label"].notna()).sum()),
        "train_surv_labeled": int((s_train["os_time"].notna() & s_train["os_event"].notna()).sum()),
        "test_surv_labeled": int((s_test["os_time"].notna() & s_test["os_event"].notna()).sum()),
        "class_label_counts_train": dict(y_train["class_label"].value_counts(dropna=True)),
        "class_label_counts_test": dict(y_test["class_label"].value_counts(dropna=True)),
    }])
    split_summary.to_csv(TABLES_DIR / "split_summary.csv", index=False)

    LOG.info("✅ Saved: %s", PROCESSED_DIR / "X_train.parquet")
    LOG.info("✅ Saved: %s", PROCESSED_DIR / "X_test.parquet")
    LOG.info("✅ Saved: %s", PROCESSED_DIR / "y_class_train.csv")
    LOG.info("✅ Saved: %s", PROCESSED_DIR / "y_class_test.csv")
    LOG.info("✅ Saved: %s", PROCESSED_DIR / "survival_train.csv")
    LOG.info("✅ Saved: %s", PROCESSED_DIR / "survival_test.csv")
    LOG.info("✅ Saved: %s", TABLES_DIR / "split_summary.csv")
    LOG.info("Split summary:\n%s", split_summary.to_string(index=False))


if __name__ == "__main__":
    main()
