from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"


def _load_processed() -> tuple[pd.DataFrame, pd.DataFrame]:
    X_train = pd.read_parquet(PROCESSED_DIR / "X_train.parquet")
    surv_train = pd.read_csv(PROCESSED_DIR / "survival_train.csv", index_col="sample_id")
    return X_train, surv_train


def _filter_survival_labeled(X: pd.DataFrame, s: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    s = s[s["os_time"].notna() & s["os_event"].notna()].copy()
    X = X.loc[s.index]
    s["os_time"] = s["os_time"].astype(float)
    s["os_event"] = s["os_event"].astype(float)
    return X, s


def _aggregate_patient(X: pd.DataFrame, s: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tmp = s[["patient_id", "os_time", "os_event"]].copy()
    tmp["patient_id"] = tmp["patient_id"].astype(str)

    Xp = X.copy()
    Xp["patient_id"] = tmp["patient_id"].values
    X_agg = Xp.groupby("patient_id").mean(numeric_only=True)

    s_agg = tmp.groupby("patient_id", as_index=True).agg(
        os_time=("os_time", "first"),
        os_event=("os_event", "first"),
    )
    return X_agg, s_agg


def _norm_gene(x) -> str | float:
    """Normalize gene identifiers for robust joining."""
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    # "1234.0" -> "1234"
    if s.endswith(".0") and s.replace(".0", "").isdigit():
        s = s.replace(".0", "")
    return s


def _ensure_gene_col(df: pd.DataFrame, desired: str = "gene") -> pd.DataFrame:
    """Ensure dataframe has a 'gene' column even if gene names were index/unnamed."""
    if desired in df.columns:
        return df
    df2 = df.copy().reset_index()
    if desired in df2.columns:
        return df2
    if "index" in df2.columns:
        return df2.rename(columns={"index": desired})
    return df2.rename(columns={df2.columns[0]: desired})


def _normalize_feature_columns(X: pd.DataFrame) -> pd.DataFrame:
    """Ensure feature columns are strings and normalized; collapse duplicates if created."""
    Xn = X.copy()
    Xn.columns = [_norm_gene(c) for c in Xn.columns]

    # If normalization created duplicate columns, collapse them by mean
    if len(set(Xn.columns)) != len(Xn.columns):
        LOG.warning("Duplicate feature columns detected after normalization; collapsing by mean.")
        Xn = Xn.T.groupby(level=0).mean().T

    return Xn


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=30, help="How many biomarkers to output.")
    parser.add_argument(
        "--survival-model-table",
        type=str,
        default=str(TABLES_DIR / "top_survival_genes.csv"),
        help="Output from train_survival_cox.py (genes + coef).",
    )
    parser.add_argument(
        "--cox-summary",
        type=str,
        default=str(TABLES_DIR / "cox_summary.csv"),
        help="Optional lifelines Cox summary saved from train_survival_cox.py.",
    )
    parser.add_argument("--out-csv", type=str, default=str(TABLES_DIR / "biomarkers_survival_table.csv"))
    parser.add_argument("--out-md", type=str, default=str(TABLES_DIR / "biomarkers_survival_table.md"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # ---------- Load top survival genes ----------
    top = pd.read_csv(args.survival_model_table)
    top = _ensure_gene_col(top, "gene")

    if "coef" not in top.columns:
        raise RuntimeError(
            f"Expected 'coef' in {args.survival_model_table}. Found columns: {top.columns.tolist()}"
        )

    top = top[["gene", "coef"]].copy()
    top["gene"] = top["gene"].map(_norm_gene)

    # Keep strongest effects (don’t assume file ordering)
    top = top.loc[top["coef"].abs().sort_values(ascending=False).index].head(args.top_n).reset_index(drop=True)

    # ---------- Load expression + survival labels ----------
    X_train, s_train = _load_processed()
    X_train, s_train = _filter_survival_labeled(X_train, s_train)

    # Patient aggregation (matches your survival training mode)
    Xp, _ = _aggregate_patient(X_train, s_train)
    Xp = _normalize_feature_columns(Xp)

    # ---------- Variance per gene (patient-level) ----------
    top_genes = top["gene"].dropna().astype(str).tolist()
    genes_present = [g for g in top_genes if g in Xp.columns]

    if not genes_present:
        raise RuntimeError(
            "None of the top genes were found in the expression matrix columns.\n"
            "This usually means you have an ID mismatch (gene symbols vs Entrez vs Ensembl).\n"
            "Quick check:\n"
            "  - print(list(X_train.columns[:10]))\n"
            "  - print(top['gene'].head(10).tolist())"
        )

    missing = sorted(set(top_genes) - set(genes_present))
    if missing:
        LOG.warning("Some top genes not found in expression matrix (first 10): %s", missing[:10])

    var = Xp[genes_present].var(axis=0).rename("variance").reset_index()
    var = var.rename(columns={"index": "gene"})
    var["gene"] = var["gene"].map(_norm_gene)

    # ---------- Build biomarker table ----------
    out = top.merge(var, on="gene", how="left")
    out["hazard_ratio"] = np.exp(out["coef"].astype(float))
    out["direction"] = np.where(out["coef"] > 0, "↑ risk (worse survival)", "↓ risk (better survival)")

    # Quality check: how many variances missing?
    miss_var = out["variance"].isna().mean() * 100
    LOG.info("Variance missing for %.1f%% of biomarkers", miss_var)

    # ---------- Optional: enrich with lifelines summary ----------
    cox_path = Path(args.cox_summary)
    if cox_path.exists():
        cox = pd.read_csv(cox_path)
        cox = _ensure_gene_col(cox, "gene")
        cox["gene"] = cox["gene"].map(_norm_gene)

        keep = ["gene"]
        for col in [
            "p",
            "exp(coef)",
            "exp(coef) lower 95%",
            "exp(coef) upper 95%",
            "coef lower 95%",
            "coef upper 95%",
        ]:
            if col in cox.columns:
                keep.append(col)

        out = out.merge(cox[keep], on="gene", how="left")
        LOG.info("✅ Enriched biomarkers with Cox summary: %s", cox_path)
    else:
        LOG.info("Cox summary not found (%s). Table will not include p-values/CI.", cox_path)

    # ---------- Sort by effect size ----------
    out["abs_log_hr"] = np.abs(np.log(out["hazard_ratio"].astype(float)))
    out = out.sort_values("abs_log_hr", ascending=False).drop(columns=["abs_log_hr"])

    # Save CSV
    out.to_csv(args.out_csv, index=False)
    LOG.info("✅ Saved: %s", args.out_csv)

    # Save Markdown (pretty formatting)
    md = out.copy()

    if "coef" in md.columns:
        md["coef"] = md["coef"].map(lambda x: f"{float(x):.4f}" if pd.notna(x) else "")
    if "hazard_ratio" in md.columns:
        md["hazard_ratio"] = md["hazard_ratio"].map(lambda x: f"{float(x):.3f}" if pd.notna(x) else "")
    if "variance" in md.columns:
        md["variance"] = md["variance"].map(lambda x: f"{float(x):.4f}" if pd.notna(x) else "")
    if "p" in md.columns:
        md["p"] = md["p"].map(lambda x: f"{float(x):.2e}" if pd.notna(x) else "")

    for ci_col in ["exp(coef) lower 95%", "exp(coef) upper 95%"]:
        if ci_col in md.columns:
            md[ci_col] = md[ci_col].map(lambda x: f"{float(x):.3f}" if pd.notna(x) else "")

    Path(args.out_md).write_text(md.to_markdown(index=False), encoding="utf-8")
    LOG.info("✅ Saved: %s", args.out_md)


if __name__ == "__main__":
    main()
