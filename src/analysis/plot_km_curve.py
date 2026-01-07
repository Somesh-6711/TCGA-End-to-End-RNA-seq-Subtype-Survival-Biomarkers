from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"
FIG_DIR = REPORTS_DIR / "figures"


def _read_indexed_csv(path: str, preferred_index: str = "sample_id") -> pd.DataFrame:
    """Read CSV and set index robustly (sample_id / Unnamed: 0 / first column)."""
    df = pd.read_csv(path)
    if preferred_index in df.columns:
        df = df.set_index(preferred_index)
    else:
        first = df.columns[0]
        if str(first).lower().startswith("unnamed") or str(first).lower() in {"index", "sample"}:
            df = df.set_index(first)
        else:
            # fallback: assume first column is an index-like column
            df = df.set_index(first)
    df.index = df.index.astype(str)
    return df


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--groups", type=str, default=str(TABLES_DIR / "km_groups_test.csv"))
    p.add_argument("--survival", type=str, default=str(PROCESSED_DIR / "survival_test.csv"))
    p.add_argument("--out", type=str, default=str(FIG_DIR / "fig3_km_risk_groups.png"))
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    g = _read_indexed_csv(args.groups, "sample_id")
    s = _read_indexed_csv(args.survival, "sample_id")

    # Normalize expected columns
    if "group" not in g.columns and "risk_group" in g.columns:
        g = g.rename(columns={"risk_group": "group"})
    if "group" not in g.columns:
        raise RuntimeError(f"'group' column not found in {args.groups}. Found: {g.columns.tolist()}")

    # If km_groups already has survival columns, don't join survival_test again
    if ("os_time" in g.columns) and ("os_event" in g.columns):
        df = g.copy()
    else:
    # join survival from file
        df = g.join(s[["os_time", "os_event"]], how="inner")

    df = df[df["os_time"].notna() & df["os_event"].notna() & df["group"].notna()].copy()

    df["os_event"] = df["os_event"].astype(int)

    # normalize group names
    df["group"] = df["group"].astype(str).str.strip().str.lower()

    low = df[df["group"].isin(["low", "low risk", "l"])]
    high = df[df["group"].isin(["high", "high risk", "h"])]

    if low.empty or high.empty:
        # maybe they are exactly "Low"/"High" already but not lowercased?
        # or different naming -> show what exists
        raise RuntimeError(
            f"Could not find both low/high groups. Unique groups: {sorted(df['group'].unique().tolist())}"
        )

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7.5, 6))

    kmf.fit(low["os_time"], event_observed=low["os_event"], label=f"Low risk (n={len(low)})")
    ax = kmf.plot(ci_show=True)

    kmf.fit(high["os_time"], event_observed=high["os_event"], label=f"High risk (n={len(high)})")
    kmf.plot(ax=ax, ci_show=True)

    res = logrank_test(
        low["os_time"], high["os_time"],
        event_observed_A=low["os_event"], event_observed_B=high["os_event"]
    )

    plt.title(f"Kaplan–Meier Survival: Risk Stratification (log-rank p={res.p_value:.3g})")
    plt.xlabel("Time (days)")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    LOG.info("✅ Saved: %s", args.out)


if __name__ == "__main__":
    main()
