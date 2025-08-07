import pathlib
import time

from plotly.io import write_image

ROOT_DIR = pathlib.Path(__file__).parent.parent.parent

OUTPUT_ROOT = ROOT_DIR / "outputs"

FIG_DIR = OUTPUT_ROOT / "figs"
FIG_DIR.mkdir(parents=True, exist_ok=True)

TABLE_DIR = OUTPUT_ROOT / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

# SUPP_TABLE_DIR = OUTPUT_ROOT / "tables_supp"
# SUPP_TABLE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = "/Users/gulugulu/clock/mammal_orthologs_hsap_1"
# RESULT_DIR = ROOT_DIR / "MutationDiseqAnalysis/results"


class pdf_writer:
    """class that handles super annoying mathjax warning box in plotly pdf's"""

    def __init__(self) -> None:
        self._done_once = False

    def __call__(self, fig, path):
        # the sleep, plus successive write, is ESSENTIAL to avoid the super annoying
        # "[MathJax]/extensions/MathMenu.js" text box error
        # but we only need to do this once
        if not self._done_once:
            write_image(fig, path)
            time.sleep(2)
            self._done_once = True

        write_image(fig, path)

write_pdf = pdf_writer()