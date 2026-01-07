from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Any, Tuple

import numpy as np
import pandas as pd

from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    classification_report,
    confusion_matrix,
)

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
    y_train = pd.read_csv(PROCESSED_DIR / "y_class_train.csv", index_col="sample_id")
    y_test = pd.read_csv(PROCESSED_DIR / "y_class_test.csv", index_col="sample_id")
    return X_train, X_test, y_train, y_test


def _top_coeffs(model: LogisticRegression, feature_names: np.ndarray, class_names: np.ndarray, top_k: int) -> pd.DataFrame:
    """
    For multinomial LR: coef_ shape = (n_classes, n_features).
    We return top +/- genes per class by absolute weight.
    """
    coef = model.coef_
    rows = []
    for ci, cname in enumerate(class_names):
        w = coef[ci]
        # top positive and negative
        pos_idx = np.argsort(w)[-top_k:][::-1]
        neg_idx = np.argsort(w)[:top_k]

        for rank, j in enumerate(pos_idx, start=1):
            rows.append({"class": cname, "direction": "positive", "rank": rank, "gene": feature_names[j], "coef": float(w[j])})
        for rank, j in enumerate(neg_idx, start=1):
            rows.append({"class": cname, "direction": "negative", "rank": rank, "gene": feature_names[j], "coef": float(w[j])})

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cv", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-genes-out", type=int, default=25)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    _ensure_dirs()

    X_train, X_test, y_train_df, y_test_df = _load()

    # Keep labeled only (just in case)
    y_train_df = y_train_df[y_train_df["class_label"].notna()].copy()
    y_test_df = y_test_df[y_test_df["class_label"].notna()].copy()

    X_train = X_train.loc[y_train_df.index]
    X_test = X_test.loc[y_test_df.index]

    groups = y_train_df["patient_id"].astype(str).values

    le = LabelEncoder()
    y_train = le.fit_transform(y_train_df["class_label"].astype(str).values)
    y_test = le.transform(y_test_df["class_label"].astype(str).values)
    class_names = le.classes_

    LOG.info("Train samples: %d | Test samples: %d | Classes: %s", X_train.shape[0], X_test.shape[0], list(class_names))

    pipe = Pipeline(steps=[
    ("scaler", StandardScaler(with_mean=False)),
    ("clf", LogisticRegression(
        solver="lbfgs",
        l1_ratio=0.0,
        max_iter=5000,
        class_weight="balanced",
        random_state=args.seed,
        )),
    ])

    # Small grid (fast, strong baseline)
    param_grid = {
        "clf__C": [0.01, 0.1, 1.0, 3.0, 10.0],
    }

    cv = GroupKFold(n_splits=args.cv)

    gs = GridSearchCV(
        estimator=pipe,
        param_grid=param_grid,
        scoring="balanced_accuracy",
        cv=cv.split(X_train, y_train, groups=groups),
        n_jobs=-1,
        refit=True,
        verbose=0,
    )
    gs.fit(X_train, y_train)

    best = gs.best_estimator_
    LOG.info("Best params: %s | Best CV balanced_acc: %.4f", gs.best_params_, gs.best_score_)

    # Test evaluation
    y_pred = best.predict(X_test)

    metrics: Dict[str, Any] = {
        "cv_n_splits": args.cv,
        "best_params": gs.best_params_,
        "best_cv_balanced_accuracy": float(gs.best_score_),
        "test_accuracy": float(accuracy_score(y_test, y_pred)),
        "test_balanced_accuracy": float(balanced_accuracy_score(y_test, y_pred)),
        "test_f1_macro": float(f1_score(y_test, y_pred, average="macro")),
        "test_f1_weighted": float(f1_score(y_test, y_pred, average="weighted")),
        "classes": list(class_names),
        "n_train": int(X_train.shape[0]),
        "n_test": int(X_test.shape[0]),
        "n_features": int(X_train.shape[1]),
    }

    # Save metrics JSON
    (METRICS_DIR / "classification_metrics.json").write_text(json.dumps(metrics, indent=2))
    LOG.info("✅ Saved: %s", METRICS_DIR / "classification_metrics.json")

    # Save classification report CSV
    report = classification_report(y_test, y_pred, target_names=class_names, output_dict=True, zero_division=0)
    report_df = pd.DataFrame(report).T
    report_df.to_csv(TABLES_DIR / "classification_test_report.csv", index=True)
    LOG.info("✅ Saved: %s", TABLES_DIR / "classification_test_report.csv")

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=np.arange(len(class_names)))
    cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
    cm_df.to_csv(TABLES_DIR / "classification_confusion_matrix.csv", index=True)
    LOG.info("✅ Saved: %s", TABLES_DIR / "classification_confusion_matrix.csv")

    # Top genes (biomarkers)
    clf: LogisticRegression = best.named_steps["clf"]
    feature_names = np.array(X_train.columns)

    top_df = _top_coeffs(clf, feature_names, class_names, top_k=args.top_genes_out)
    top_df.to_csv(TABLES_DIR / "top_genes_coefficients.csv", index=False)
    LOG.info("✅ Saved: %s", TABLES_DIR / "top_genes_coefficients.csv")

    # Print quick view
    LOG.info(
        "Test metrics: acc=%.3f bal_acc=%.3f f1_macro=%.3f",
        metrics["test_accuracy"], metrics["test_balanced_accuracy"], metrics["test_f1_macro"]
    )


if __name__ == "__main__":
    main()
