"""Input/output paths for the Python figure scripts (mirror of R/paths.R).

REPO_ROOT is set by run_all.sh; without it the repo root is taken as two levels above this file.
"""
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT") or Path(__file__).resolve().parents[1])
RESULTS = REPO_ROOT / "results"
RNA_NS_DIR = RESULTS / "rnaseq_de" / "gene_NS"        # = Haydar_RNAseq/analysis/deseq2/results/gene_NS
DMGSEA_DIR = RESULTS / "methylation" / "dmgsea"        # = Haydar_Methylation_Array/results/dmgsea
INT_DIR = RESULTS / "integration"                      # = integration/results
FIG_ROOT = REPO_ROOT / "figures"                       # = integration/paper_figures
FONT_DIR = REPO_ROOT / "workflow" / "05_figures" / "fonts"
CANDIDATE_GENE_SETS_R = REPO_ROOT / "workflow" / "04_integration" / "candidate_gene_sets.R"
