#!/usr/bin/env bash
# Reproduce every number, table and figure in the manuscript and Additional file 1.
#
#   ./run_all.sh                         full run (about 1 h on two cores)
#   DEPBENCH_NPERM=50 DEPBENCH_NBOOT=200 DEPBENCH_SIMREPS=10 DEPBENCH_DESIGNREPS=30 ./run_all.sh
#                                        smoke test (about 10 min); Monte Carlo
#                                        quantities will differ from the paper
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
for f in data/ClinVar_All_DCS_A.csv data/VUS_Clinical_Reclassification.csv \
         data/case_study_NARS1_PGK1.csv data/loo_predictions.csv; do
  [[ -f "$f" ]] || { echo "ERROR: $f not found (see data/README.md)" >&2; exit 1; }
done
run () { echo ">>> $1"; python3 "src/$1"; }
run 00_build_cohort.py      # cohort, reconstruction of score components, encodings
run 01_ablation.py          # 2^4 logistic and 2^3 three-family ablation
run 02_diagnostics.py       # association, pair composition, baselines, permutation nulls
run 03_normalisation.py     # provenance, gene identity, transfer, increment over AlphaMissense
run 03b_rank_convention.py  # percentile ranking convention
run 05_ensemble_audit.py    # exploratory ensemble: weights, nested blending
run 06_within_gene.py       # within-gene tests, conditional logistic regression
run 04_simulation.py        # simulation with known estimands
run 07_design.py            # genes required for a stable estimate
run 08_tables.py            # Tables S1-S2
run 11_calibration.py       # calibration, fold corrections, leave-two-genes-out
run 12_robustness.py        # gene-drop influence, power, paired within-gene tests
run 09_figures.py           # Figures 1-5
run 10_supp_figures.py      # Figures S1-S4
run 13_figure6.py           # Figure 6
python3 manuscript/build_manuscript.py
python3 manuscript/build_supplement.py
node manuscript/build_docx.js manuscript/manuscript.json manuscript/Manuscript.docx
node manuscript/build_docx.js manuscript/supplement.json manuscript/AdditionalFile1.docx --supplement
python3 src/audit.py
