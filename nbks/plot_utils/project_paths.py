import pathlib
import time

from plotly.io import write_image

ROOT_DIR = pathlib.Path(__file__)
while not (ROOT_DIR / "PuningAnalysis").exists():
    ROOT_DIR = ROOT_DIR.parent

ROOT_DIR = ROOT_DIR / "PuningAnalysis"

OUTPUT_ROOT = ROOT_DIR / "outputs"

FIG_DIR = OUTPUT_ROOT / "figs"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SUPP_FIG_DIR = OUTPUT_ROOT / "supp_figs"
SUPP_FIG_DIR.mkdir(parents=True, exist_ok=True)

TABLE_DIR = OUTPUT_ROOT / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

# SUPP_TABLE_DIR = OUTPUT_ROOT / "tables_supp"
# SUPP_TABLE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = ROOT_DIR / "data"


class pdf_writer:
    """class that handles super annoying mathjax warning box in plotly pdf's"""

    def __init__(self, width=None, height=None) -> None:
        self._done_once = False

    def __call__(self, fig, path, width=None, height=None):
        # the sleep, plus successive write, is ESSENTIAL to avoid the super annoying
        # "[MathJax]/extensions/MathMenu.js" text box error
        # but we only need to do this once
        if not self._done_once:
            write_image(fig, path)
            time.sleep(2)
            self._done_once = True

        write_image(fig, path, width=width, height=height)


write_pdf = pdf_writer(width=600, height=800)
