"""Shared configuration for the Python parts of the pipeline."""
from pathlib import Path
import os

def find_repo_root(start: Path = None) -> Path:
    p = (start or Path.cwd()).resolve()
    for cand in [p, *p.parents]:
        if (cand / ".repo-root").exists():
            return cand
    raise RuntimeError("Could not locate the repository root (.repo-root marker).")

ROOT = Path(os.environ.get("PBRAS_ROOT") or find_repo_root(Path(__file__).parent))
DATA, COUNTS, ANNOT = ROOT / "data", ROOT / "data/counts", ROOT / "data/annotation"
RESULTS = ROOT / "results"
DE_DIR, GS_DIR, FIG_DIR = RESULTS / "de_tables", RESULTS / "geneset", RESULTS / "figures"
for d in (RESULTS, DE_DIR, GS_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    "KO1": ("control",   "KO-1",   19.28),
    "KO3": ("control",   "KO-3",   18.58),
    "WT1": ("treatment", "3610-1", 57.80),
    "WT2": ("treatment", "3610-2", 77.75),
}
DESIGNS = {
    "d22a": dict(samples=["KO1", "KO3", "WT1", "WT2"],
                 label="2 WT vs 2 KO", short="2WTx2KO"),
}
DESIGN = DESIGNS["d22a"]
FTM = {"PMG11_03146":"FtmA","PMG11_03147":"FtmC","PMG11_03148":"FtmD",
       "PMG11_03149":"FtmB","PMG11_03150":"FtmE","PMG11_03151":"FtmF",
       "PMG11_03152":"FtmG","PMG11_03153":"FtmH"}
FTM_ORDER = ["FtmA","FtmB","FtmC","FtmD","FtmE","FtmF","FtmG","FtmH"]
FC_CUTOFF, FDR_CUTOFF = 2, 0.05
