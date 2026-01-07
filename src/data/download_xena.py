from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Tuple, Optional, List

import numpy as np
import pandas as pd
import requests

LOG = logging.getLogger(__name__)

TCGA_HUB = "https://tcga.xenahubs.net"
PANCAN_HUB = "https://pancanatlas.xenahubs.net"

LGG_EXPR = "TCGA.LGG.sampleMap/HiSeqV2"
LGG_CLIN = "TCGA.LGG.sampleMap/LGG_clinicalMatrix"
LGG_SURV = "survival/LGG_survival.txt"
PANCAN_SUBTYPE = "TCGASubtype.20170308.tsv"

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw" / "xena"
INTERIM_DIR = DATA_DIR / "interim"
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"

_TCGA_PAT_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", re.IGNORECASE)
_TCGA_SAMP_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2})", re.IGNORECASE)

def _patient_id(x: str) -> str:
    s = str(x)
    m = _TCGA_PAT_RE.search(s)
    return m.group(1).upper() if m else s

def _sample_short(x: str) -> str:
    s = str(x)
    m = _TCGA_SAMP_RE.search(s)
    return m.group(1).upper() if m else s

def _ensure_dirs() -> None:
    for p in [RAW_DIR, INTERIM_DIR, TABLES_DIR]:
        p.mkdir(parents=True, exist_ok=True)

def _http_get(url: str, timeout: int = 120) -> requests.Response:
    headers = {
        "User-Agent": "tcga-xena-client/1.0",
        "Accept": "*/*",
        "Referer": "https://xenabrowser.net/",
    }
    return requests.get(url, headers=headers, timeout=timeout, stream=True)

def _download(url: str, out_path: Path, force: bool) -> None:
    if out_path.exists() and not force:
        LOG.info("Exists (skip): %s", out_path)
        return
    LOG.info("Downloading: %s", url)
    r = _http_get(url)
    if r.status_code != 200:
        raise RuntimeError(f"HTTP {r.status_code} for {url}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)
    LOG.info("Saved -> %s", out_path)

def _download_xena_dataset(hub: str, dataset_id: str, out_path: Path, force: bool) -> None:
    candidates = [
        f"{hub}/download/{dataset_id}.gz",
        f"{hub}/download/{dataset_id}.tsv.gz",
        f"{hub}/download/{dataset_id}",
    ]
    last_err = None
    for url in candidates:
        try:
            _download(url, out_path, force=force)
            return
        except Exception as e:
            last_err = e
            continue
    raise RuntimeError(f"Failed to download dataset {dataset_id}. Tried: {candidates}") from last_err

def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False

def _read_table(path: Path) -> pd.DataFrame:
    # auto-detect gzip content (Xena sometimes returns plain text)
    if _is_gzip(path):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)

def _clean_nan(df: pd.DataFrame) -> pd.DataFrame:
    return df.replace({"NaN": np.nan, "nan": np.nan, "NA": np.nan, "": np.nan}, regex=False)

def _read_expression_hi_seq_v2(path_gz: Path) -> pd.DataFrame:
    LOG.info("Reading expression matrix...")
    df = pd.read_csv(path_gz, sep="\t", compression="gzip", low_memory=False)
    feature_col = df.columns[0]
    df = df.set_index(feature_col)
    df = df.apply(pd.to_numeric, errors="coerce")
    X = df.T
    X.index.name = "sample_id"
    return X.astype("float32")

def _dedup_map(keys: pd.Series, values: pd.Series) -> pd.Series:
    df = pd.DataFrame({"k": keys.astype(str), "v": values})
    df["v"] = df["v"].replace({"NaN": np.nan, "nan": np.nan, "NA": np.nan, "": np.nan}, regex=False)

    def first_nonnull(x: pd.Series):
        x = x.dropna()
        return x.iloc[0] if len(x) else np.nan

    out = df.groupby("k", as_index=True)["v"].agg(first_nonnull)
    out.index = out.index.astype(str)
    return out

def _infer_os_time_and_event(samples: list[str], surv_df: pd.DataFrame, clin_df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    surv = surv_df.copy()
    clin = clin_df.copy()

    id_col = surv.columns[0]
    surv["_sid"] = surv[id_col].astype(str).apply(_sample_short)
    surv["_pid"] = surv[id_col].astype(str).apply(_patient_id)

    time_candidates = [c for c in surv.columns if c.lower() in ["os.time", "os_time"]]
    if not time_candidates:
        time_candidates = [c for c in surv.columns if ("os" in c.lower() and "time" in c.lower())]
    if not time_candidates:
        raise RuntimeError(f"Could not find OS.time in survival file. Columns: {list(surv.columns)[:30]}")
    os_time_col = time_candidates[0]

    event_candidates = [c for c in surv.columns if c.lower() in ["os", "os.event", "os_event", "event"]]
    event_candidates = [c for c in event_candidates if c != os_time_col]
    if not event_candidates:
        event_candidates = [c for c in surv.columns if "event" in c.lower() or ("os" in c.lower() and "status" in c.lower())]
    os_event_surv_col = event_candidates[0] if event_candidates else None

    sample_short = pd.Series([_sample_short(s) for s in samples], index=samples)
    patient_ids  = pd.Series([_patient_id(s) for s in samples], index=samples)

    map_sid = _dedup_map(surv["_sid"], surv[os_time_col])
    map_pid = _dedup_map(surv["_pid"], surv[os_time_col])

    os_time = sample_short.map(map_sid).fillna(patient_ids.map(map_pid))
    os_time = pd.to_numeric(os_time, errors="coerce")

    os_event = pd.Series(np.nan, index=samples, dtype="float")

    if os_event_surv_col is not None:
        emap_sid = _dedup_map(surv["_sid"], surv[os_event_surv_col])
        emap_pid = _dedup_map(surv["_pid"], surv[os_event_surv_col])
        tmp = sample_short.map(emap_sid).fillna(patient_ids.map(emap_pid))
        tmp = tmp.astype(str).str.strip().str.lower().replace({"dead": "1", "deceased": "1", "alive": "0"})
        os_event = pd.to_numeric(tmp, errors="coerce")

    # clinicalMatrix vital_status fallback
    clin_id_col = clin.columns[0]
    clin["_sid"] = clin[clin_id_col].astype(str).apply(_sample_short)
    clin["_pid"] = clin[clin_id_col].astype(str).apply(_patient_id)

    vital_cols = [c for c in clin.columns if c.lower() == "vital_status" or "vital" in c.lower()]
    if vital_cols:
        vcol = vital_cols[0]
        vmap_sid = _dedup_map(clin["_sid"], clin[vcol])
        vmap_pid = _dedup_map(clin["_pid"], clin[vcol])

        vital = sample_short.map(vmap_sid).fillna(patient_ids.map(vmap_pid))
        vital = vital.astype(str).str.strip().str.lower().replace({"dead": "1", "deceased": "1", "alive": "0"})
        vital_event = pd.to_numeric(vital, errors="coerce")

        os_event = os_event.fillna(vital_event)

    return os_time, os_event

def _choose_class_label(
    clin_df: pd.DataFrame,
    min_coverage: float = 0.60,
    max_classes: int = 12,
) -> Tuple[Optional[str], Optional[pd.Series]]:
    """
    Pick a clinically meaningful categorical label for classification.

    Strategy:
    1) Try high-value targets first (histology-related).
    2) Only fall back to easy demographics like gender if nothing else works.
    """
    high_value = [
        "histological_type",
        "icd_o_3_histology",
        "icd_10",
        "icd_o_3_site",
        "laterality",
    ]
    fallback = [
        "eastern_cancer_oncology_group",
        "karnofsky_performance_score",  # might be numeric-like; if too many unique, rejected
        "inherited_genetic_syndrome_found",
        "family_history_of_primary_brain_tumor",
        "gender",  # LAST resort
    ]

    def prep(col: str) -> pd.Series:
        s = clin_df[col].copy()
        s = s.replace({"NaN": np.nan, "nan": np.nan, "NA": np.nan, "": np.nan}, regex=False)
        s = s.astype(str).str.strip()
        s = s.replace({"nan": np.nan, "NaN": np.nan, "NA": np.nan, "": np.nan})
        # restore NaN where original was NaN
        s = s.where(clin_df[col].notna(), np.nan)
        return s

    def acceptable(s: pd.Series) -> bool:
        coverage = float(s.notna().mean())
        if coverage < min_coverage:
            return False
        n_classes = int(s.dropna().nunique())
        if n_classes < 2 or n_classes > max_classes:
            return False
        return True

    # Pass 1: high-value columns
    for c in high_value:
        if c in clin_df.columns:
            s = prep(c)
            if acceptable(s):
                return c, s

    # Pass 2: fallback columns
    for c in fallback:
        if c in clin_df.columns:
            s = prep(c)
            if acceptable(s):
                return c, s

    return None, None

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True, choices=["LGG"])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    _ensure_dirs()

    lgg_dir = RAW_DIR / "LGG"
    pancan_dir = RAW_DIR / "pancan"
    lgg_dir.mkdir(parents=True, exist_ok=True)
    pancan_dir.mkdir(parents=True, exist_ok=True)

    expr_path   = lgg_dir / "expression.gz"
    clin_path   = lgg_dir / "clinicalMatrix.maybe_gz"
    surv_path   = lgg_dir / "survival.maybe_gz"
    pancan_path = pancan_dir / "TCGASubtype.20170308.tsv.maybe_gz"

    _download_xena_dataset(TCGA_HUB, LGG_EXPR, expr_path, force=args.force)
    _download_xena_dataset(TCGA_HUB, LGG_CLIN, clin_path, force=args.force)
    _download_xena_dataset(TCGA_HUB, LGG_SURV, surv_path, force=args.force)
    _download_xena_dataset(PANCAN_HUB, PANCAN_SUBTYPE, pancan_path, force=args.force)

    X = _read_expression_hi_seq_v2(expr_path)
    samples = X.index.tolist()
    LOG.info("Expression samples: %d | Cohort: LGG", len(samples))

    clin_df = _clean_nan(_read_table(clin_path))
    surv_df = _clean_nan(_read_table(surv_path))
    pancan_df = _clean_nan(_read_table(pancan_path))

    # ---- survival labels ----
    os_time, os_event = _infer_os_time_and_event(samples, surv_df, clin_df)

    # ---- optional pancan subtype ----
    subtype_primary = pd.Series(np.nan, index=samples, name="subtype_primary")
    if "sampleID" in pancan_df.columns:
        primary_col = "Subtype_Integrative" if "Subtype_Integrative" in pancan_df.columns else (
            "Subtype_Selected" if "Subtype_Selected" in pancan_df.columns else (
                "Subtype_mRNA" if "Subtype_mRNA" in pancan_df.columns else None
            )
        )
        if primary_col:
            pancan_df["patient_id"] = pancan_df["sampleID"].astype(str).apply(_patient_id)
            pm = pancan_df.groupby("patient_id")[primary_col].apply(
                lambda s: s.dropna().iloc[0] if s.dropna().shape[0] else np.nan
            )
            subtype_primary = pd.Series([pm.get(_patient_id(s), np.nan) for s in samples],
                                        index=samples, name="subtype_primary")

    # ---- choose classification label from clinicalMatrix ----
    class_col, class_series = _choose_class_label(clin_df, min_coverage=0.70, max_classes=10)
    if class_col is None:
        LOG.warning("No good classification label found in clinicalMatrix with current thresholds.")
        class_label = pd.Series(np.nan, index=samples, name="class_label")
        n_classes = 0
    else:
        LOG.info("Classification label selected: %s", class_col)
        # map clinical rows to cohort samples by sample_short/patient_id
        # clinicalMatrix rows are per sampleID (often TCGA-..-..-01)
        clin_id_col = "sampleID" if "sampleID" in clin_df.columns else clin_df.columns[0]
        clin_df["_sid"] = clin_df[clin_id_col].astype(str).apply(_sample_short)
        clin_df["_pid"] = clin_df[clin_id_col].astype(str).apply(_patient_id)

        cohort_sid = pd.Series([_sample_short(s) for s in samples], index=samples)
        cohort_pid = pd.Series([_patient_id(s) for s in samples], index=samples)

        m_sid = _dedup_map(clin_df["_sid"], class_series)
        m_pid = _dedup_map(clin_df["_pid"], class_series)

        class_label = cohort_sid.map(m_sid).fillna(cohort_pid.map(m_pid))
        class_label.name = "class_label"
        n_classes = int(class_label.dropna().nunique())

    # ---- build outputs ----
    clinical = pd.DataFrame(index=pd.Index(samples, name="sample_id"))
    clinical["patient_id"] = [_patient_id(s) for s in samples]
    clinical["os_time"] = os_time.values
    clinical["os_event"] = os_event.values
    clinical["subtype_primary"] = subtype_primary.values
    clinical["class_label"] = class_label.values
    clinical["class_label_name"] = class_col if class_col else np.nan

    labels = pd.DataFrame(index=pd.Index(samples, name="sample_id"))
    labels["patient_id"] = clinical["patient_id"].values
    labels["os_time"] = os_time.values
    labels["os_event"] = os_event.values
    labels["subtype_primary"] = subtype_primary.values
    labels["class_label"] = class_label.values

    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    X.to_parquet(INTERIM_DIR / "expression.parquet")
    clinical.to_parquet(INTERIM_DIR / "clinical.parquet")
    labels.to_parquet(INTERIM_DIR / "labels.parquet")

    os_labeled = int(((labels["os_time"].notna()) & (labels["os_event"].notna())).sum())
    subtype_labeled = int(labels["subtype_primary"].notna().sum())
    class_labeled = int(labels["class_label"].notna().sum())

    summary = pd.DataFrame([{
        "n_samples": len(labels),
        "os_labeled": os_labeled,
        "class_label_name": class_col if class_col else "",
        "class_labeled": class_labeled,
        "n_classes": n_classes,
        "subtype_labeled": subtype_labeled,
    }])
    summary.to_csv(TABLES_DIR / "label_summary.csv", index=False)

    LOG.info("✅ Saved: %s", INTERIM_DIR / "expression.parquet")
    LOG.info("✅ Saved: %s", INTERIM_DIR / "clinical.parquet")
    LOG.info("✅ Saved: %s", INTERIM_DIR / "labels.parquet")
    LOG.info("✅ Saved: %s", TABLES_DIR / "label_summary.csv")
    LOG.info("Label summary:\n%s", summary.to_string(index=False))

if __name__ == "__main__":
    main()
