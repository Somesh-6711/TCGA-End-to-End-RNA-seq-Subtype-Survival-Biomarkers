from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Any, Tuple

import numpy as np
import pandas as pd

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"
METRICS_DIR = REPORTS_DIR / "metrics"


def _ensure_dirs() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    METRICS_DIR.mkdir(parents=True, exist_ok=True)


def _load() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    X_train = pd.read_parquet(PROCESSED_DIR / "X_train.parquet")
    X_test = pd.read_parquet(PROCESSED_DIR / "X_test.parquet")
    s_train = pd.read_csv(PROCESSED_DIR / "survival_train.csv", index_col="sample_id")
    s_test = pd.read_csv(PROCESSED_DIR / "survival_test.csv", index_col="sample_id")
    return X_train, X_test, s_train, s_test


def _filter_labeled(X: pd.DataFrame, s: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    s = s[s["os_time"].notna() & s["os_event"].notna()].copy()
    X = X.loc[s.index]
    # ensure numeric
    s["os_time"] = s["os_time"].astype(float)
    s["os_event"] = s["os_event"].astype(float)
    return X, s


def _aggregate_by_patient(X: pd.DataFrame, s: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate multiple samples per patient to one row:
    - expression: mean
    - labels: first (should be identical within patient for TCGA OS)
    """
    tmp = s[["patient_id", "os_time", "os_event"]].copy()
    tmp["patient_id"] = tmp["patient_id"].astype(str)

    # join patient_id onto X
    Xp = X.copy()
    Xp["patient_id"] = tmp["patient_id"].values

    X_agg = Xp.groupby("patient_id").mean(numeric_only=True)

    s_agg = tmp.groupby("patient_id", as_index=True).agg(
        os_time=("os_time", "first"),
        os_event=("os_event", "first"),
    )
    return X_agg, s_agg


def _select_top_variance(X_train: pd.DataFrame, X_test: pd.DataFrame, top_genes: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if top_genes <= 0 or top_genes >= X_train.shape[1]:
        return X_train, X_test
    var = X_train.var(axis=0)
    keep = var.sort_values(ascending=False).head(top_genes).index
    return X_train[keep], X_test[keep]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--penalizer", type=float, default=0.5, help="L2 penalizer for CoxPH (stability).")
    parser.add_argument("--top-genes", type=int, default=500, help="Train-only variance selection for survival.")
    parser.add_argument("--patient-agg", action="store_true", help="Aggregate samples to patient-level (recommended).")
    parser.add_argument("--top-genes-out", type=int, default=30)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    _ensure_dirs()

    X_train, X_test, s_train, s_test = _load()
    X_train, s_train = _filter_labeled(X_train, s_train)
    X_test, s_test = _filter_labeled(X_test, s_test)

    LOG.info("Raw survival rows: train=%d test=%d", X_train.shape[0], X_test.shape[0])

    if args.patient_agg:
        X_train, s_train = _aggregate_by_patient(X_train, s_train)
        X_test, s_test = _aggregate_by_patient(X_test, s_test)
        LOG.info("Patient-aggregated rows: train=%d test=%d", X_train.shape[0], X_test.shape[0])

    # Feature selection (train-only, leakage-safe)
    X_train, X_test = _select_top_variance(X_train, X_test, args.top_genes)
    LOG.info("Using features: %d (top-variance=%d)", X_train.shape[1], args.top_genes)

    # Standardize (train only)
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xtr = pd.DataFrame(scaler.fit_transform(X_train.values), index=X_train.index, columns=X_train.columns)
    Xte = pd.DataFrame(scaler.transform(X_test.values), index=X_test.index, columns=X_test.columns)

    df_tr = Xtr.copy()
    df_tr["os_time"] = s_train["os_time"].values
    df_tr["os_event"] = s_train["os_event"].values

    df_te = Xte.copy()
    df_te["os_time"] = s_test["os_time"].values
    df_te["os_event"] = s_test["os_event"].values

    cph = CoxPHFitter(penalizer=args.penalizer)
    cph.fit(df_tr, duration_col="os_time", event_col="os_event")

    risk_train = cph.predict_partial_hazard(df_tr).rename("risk_score")
    risk_test = cph.predict_partial_hazard(df_te).rename("risk_score")

    # Flip sign: higher risk -> shorter survival
    cindex_train = float(concordance_index(df_tr["os_time"], -risk_train.values, df_tr["os_event"]))
    cindex_test = float(concordance_index(df_te["os_time"], -risk_test.values, df_te["os_event"]))

    metrics: Dict[str, Any] = {
        "model": "CoxPH",
        "penalizer": float(args.penalizer),
        "top_genes": int(args.top_genes),
        "patient_agg": bool(args.patient_agg),
        "n_train": int(df_tr.shape[0]),
        "n_test": int(df_te.shape[0]),
        "n_features": int(X_train.shape[1]),
        "cindex_train": cindex_train,
        "cindex_test": cindex_test,
    }

    (METRICS_DIR / "survival_metrics.json").write_text(json.dumps(metrics, indent=2))
    LOG.info("✅ Saved: %s", METRICS_DIR / "survival_metrics.json")
    LOG.info("C-index: train=%.3f | test=%.3f", cindex_train, cindex_test)

    out = pd.DataFrame(
        {
            "os_time": df_te["os_time"].values,
            "os_event": df_te["os_event"].values,
            "risk_score": risk_test.values,
        },
        index=df_te.index,
    )
    out.index.name = "patient_id" if args.patient_agg else "sample_id"
    out.to_csv(TABLES_DIR / "survival_risk_scores_test.csv", index=True)
    LOG.info("✅ Saved: %s", TABLES_DIR / "survival_risk_scores_test.csv")

    coefs = cph.params_.sort_values(key=lambda s: s.abs(), ascending=False)
    top = coefs.head(args.top_genes_out).to_frame(name="coef")
    top.index.name = "gene"
    top.to_csv(TABLES_DIR / "top_survival_genes.csv", index=True)
    LOG.info("✅ Saved: %s", TABLES_DIR / "top_survival_genes.csv")


if __name__ == "__main__":
    main()
