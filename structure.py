

from __future__ import annotations

from pathlib import Path

PROJECT_NAME = "tcga-translational-genomics"
PYTHON_PACKAGE = "src"


FOLDERS = [
    "data/raw",
    "data/interim",
    "data/processed",
    "data/external",
    "notebooks",
    "reports/figures",
    "reports/tables",
    "app",
    f"{PYTHON_PACKAGE}/data",
    f"{PYTHON_PACKAGE}/features",
    f"{PYTHON_PACKAGE}/models",
    f"{PYTHON_PACKAGE}/viz",
    f"{PYTHON_PACKAGE}/utils",
    "tests",
]


FILES = {
    "README.md": f"# {PROJECT_NAME}\n\nEnd-to-end TCGA translational genomics pipeline.\n",
    ".gitignore": "\n".join(
        [
            "# Python",
            "__pycache__/",
            "*.pyc",
            ".pytest_cache/",
            ".mypy_cache/",
            ".ruff_cache/",
            "",
            "# Environments",
            ".venv/",
            "venv/",
            "env/",
            "",
            "# Data (keep folders, ignore contents)",
            "data/raw/*",
            "data/interim/*",
            "data/processed/*",
            "data/external/*",
            "!data/**/.gitkeep",
            "",
            "# Notebooks",
            ".ipynb_checkpoints/",
            "",
            "# OS",
            ".DS_Store",
            "Thumbs.db",
            "",
            "# App secrets",
            ".streamlit/secrets.toml",
            "",
            "# Reports build artifacts",
            "reports/**/*.html",
        ]
    )
    + "\n",
    "requirements.txt": "# filled later by you\n",
    "app/streamlit_app.py": """\
import streamlit as st

st.set_page_config(page_title="TCGA Translational Genomics", layout="wide")

st.title("TCGA Translational Genomics Dashboard")
st.write("Coming soon: risk scoring + explainability + pathway views.")
""",
    f"{PYTHON_PACKAGE}/__init__.py": "",
    f"{PYTHON_PACKAGE}/config.py": """\
from dataclasses import dataclass

@dataclass(frozen=True)
class Config:
    random_state: int = 42
    test_size: float = 0.2
""",
    f"{PYTHON_PACKAGE}/paths.py": """\
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"

REPORTS_DIR = PROJECT_ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
TABLES_DIR = REPORTS_DIR / "tables"
""",
    f"{PYTHON_PACKAGE}/utils/logger.py": """\
import logging

def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        fmt = logging.Formatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    return logger
""",
    f"{PYTHON_PACKAGE}/data/download_xena.py": """\
\"\"\"Download TCGA data from UCSC Xena (placeholder).

Next step: we will fill this with a robust downloader and dataset IDs.
\"\"\"

def main() -> None:
    raise NotImplementedError("We will implement this next.")

if __name__ == "__main__":
    main()
""",
    "tests/test_imports.py": f"""\
def test_imports():
    import {PYTHON_PACKAGE}  # noqa: F401
""",
}

GITKEEP_FOLDERS = [
    "data/raw",
    "data/interim",
    "data/processed",
    "data/external",
]


def write_if_missing(path: Path, content: str) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def main() -> None:
    root = Path.cwd()

    # Create folders
    for folder in FOLDERS:
        (root / folder).mkdir(parents=True, exist_ok=True)

    # Add .gitkeep so empty folders survive git
    for folder in GITKEEP_FOLDERS:
        keep = root / folder / ".gitkeep"
        write_if_missing(keep, "")

    # Create files
    for rel_path, content in FILES.items():
        write_if_missing(root / rel_path, content)

    print("✅ Project scaffold created.")
    print("Next: edit requirements.txt, create environment, implement download script.")


if __name__ == "__main__":
    main()
