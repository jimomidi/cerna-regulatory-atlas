"""
Shared helper functions used across the analysis and figure scripts.

Centralizing these small utilities avoids repeating the same plot-styling
and I/O boilerplate in every script.
"""

from pathlib import Path
from typing import Union

import matplotlib as mpl

from config import FONT_FAMILY, RESULTS_DIR


def apply_default_style(font_size: int = 12) -> None:
    """Apply the project-wide matplotlib style (font family and base size).

    Falls back gracefully if the requested font is not installed, since
    matplotlib will simply substitute its default sans/serif font.
    """
    mpl.rcParams["font.family"] = FONT_FAMILY
    mpl.rcParams["font.size"] = font_size
    mpl.rcParams["axes.labelweight"] = "bold"


def ensure_dir(path: Union[str, Path]) -> Path:
    """Create ``path`` (and parents) if it does not already exist, returning it."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def results_path(*parts: str) -> Path:
    """Build a path under the project's results directory, creating it if needed."""
    ensure_dir(RESULTS_DIR)
    return RESULTS_DIR.joinpath(*parts)
