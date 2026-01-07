from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

LOG = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORTS_DIR = REPO_ROOT / "reports"
TABLES_DIR = REPORTS_DIR / "tables"
FIGURES_DIR = REPORTS_DIR / "figures"
METRICS_DIR = REPORTS_DIR / "metrics"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--groups-file",
        type=str,
        default=str(TABLES_DIR / "km_groups_test.csv"),
        help="CSV produced by km_risk_stratification.py",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=str(FIGURES_DIR / "km_test.png"),
        help="Output image path",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.groups_file)
    required = {"os_time", "os_event", "risk_group"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing columns in {args.groups_file}: {sorted(missing)}")

    # Only Low/High for the median split
    df = df[df["risk_group"].isin(["Low", "High"])].copy()

    low = df[df["risk_group"] == "Low"]
    high = df[df["risk_group"] == "High"]

    lr = logrank_test(
        low["os_time"], high["os_time"],
        event_observed_A=low["os_event"], event_observed_B=high["os_event"]
    )
    p = lr.p_value

    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()

    plt.figure()
    ax = plt.gca()

    kmf_low.fit(low["os_time"], event_observed=low["os_event"], label=f"Low risk (n={len(low)})")
    kmf_low.plot_survival_function(ax=ax)

    kmf_high.fit(high["os_time"], event_observed=high["os_event"], label=f"High risk (n={len(high)})")
    kmf_high.plot_survival_function(ax=ax)

    ax.set_title(f"Kaplan–Meier Survival Curves (Test) | Log-rank p={p:.3g}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    ax.grid(True, alpha=0.25)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    LOG.info("✅ Saved: %s", out_path)


if __name__ == "__main__":
    main()
