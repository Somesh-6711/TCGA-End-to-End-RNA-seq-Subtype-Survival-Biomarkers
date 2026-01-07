from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Any

import numpy as np
import pandas as pd
from lifelines.statistics import logrank_test

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"
METRICS_DIR = REPORTS_DIR / "metrics"


def _ensure_dirs() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    METRICS_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--risk-file",
        type=str,
        default=str(TABLES_DIR / "survival_risk_scores_test.csv"),
        help="CSV produced by train_survival_cox.py",
    )
    parser.add_argument(
        "--split",
        choices=["median", "tertile", "quartile"],
        default="median",
        help="How to form risk groups.",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    _ensure_dirs()

    df = pd.read_csv(args.risk_file)

    # Determine ID column name (patient_agg -> patient_id, otherwise sample_id)
    id_col = "patient_id" if "patient_id" in df.columns else "sample_id"
    required = {id_col, "os_time", "os_event", "risk_score"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing required columns in {args.risk_file}: {sorted(missing)}")

    df = df[[id_col, "os_time", "os_event", "risk_score"]].copy()
    df = df[df["os_time"].notna() & df["os_event"].notna() & df["risk_score"].notna()].copy()

    # Build groups
    r = df["risk_score"].astype(float)

    if args.split == "median":
        thr = float(r.median())
        df["risk_group"] = np.where(r >= thr, "High", "Low")
        group_info = {"method": "median", "threshold": thr}

        high = df[df["risk_group"] == "High"]
        low = df[df["risk_group"] == "Low"]

        lr = logrank_test(
            low["os_time"], high["os_time"],
            event_observed_A=low["os_event"], event_observed_B=high["os_event"]
        )

        metrics: Dict[str, Any] = {
            "split": group_info,
            "n_total": int(df.shape[0]),
            "n_low": int(low.shape[0]),
            "n_high": int(high.shape[0]),
            "logrank_p_value": float(lr.p_value),
            "test_statistic": float(lr.test_statistic),
        }

    elif args.split == "tertile":
        q1, q2 = r.quantile([1/3, 2/3]).tolist()
        df["risk_group"] = pd.cut(r, bins=[-np.inf, q1, q2, np.inf], labels=["Low", "Mid", "High"])
        group_info = {"method": "tertile", "q1": float(q1), "q2": float(q2)}
        metrics = {"split": group_info, "n_total": int(df.shape[0]), "counts": df["risk_group"].value_counts().to_dict()}

    else:  # quartile
        q1, q2, q3 = r.quantile([0.25, 0.5, 0.75]).tolist()
        df["risk_group"] = pd.cut(r, bins=[-np.inf, q1, q2, q3, np.inf], labels=["Q1", "Q2", "Q3", "Q4"])
        group_info = {"method": "quartile", "q1": float(q1), "q2": float(q2), "q3": float(q3)}
        metrics = {"split": group_info, "n_total": int(df.shape[0]), "counts": df["risk_group"].value_counts().to_dict()}

    # Save grouped table
    out_path = TABLES_DIR / "km_groups_test.csv"
    df.to_csv(out_path, index=False)
    LOG.info("✅ Saved: %s", out_path)

    # Save metrics
    met_path = METRICS_DIR / "km_metrics.json"
    met_path.write_text(json.dumps(metrics, indent=2))
    LOG.info("✅ Saved: %s", met_path)

    if args.split == "median":
        LOG.info("Log-rank p-value (Low vs High): %.4g", metrics["logrank_p_value"])


if __name__ == "__main__":
    main()
