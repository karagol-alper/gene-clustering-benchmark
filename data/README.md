# Input data

Everything in `results/` is regenerated from these five files by `run_all.sh`.
No analysis script reads labels from anything except `ClinVar_All_DCS_A.csv`.

| File | Content | Used by |
|---|---|---|
| `ClinVar_All_DCS_A.csv` | 768 ClinVar missense variants in 28 ATLAS proteins, with the published DCS, RCS and normalised DCS, and both AlphaMissense scores | `00_build_cohort.py` |
| `VUS_Clinical_Reclassification.csv` | residues of uncertain significance with the raw ingredients of the scores (site conservation, per-replica RMSF and N_eq, normalised mobility, raw DCS) | `00_build_cohort.py` |
| `case_study_NARS1_PGK1.csv` | the same ingredients for two further proteins | `00_build_cohort.py` |
| `loo_predictions.csv` | leave-one-variant-out predictions of the five base learners and the blend of the exploratory ensemble | `05_ensemble_audit.py` |
| `AdditionalFile2_cohort_v1.csv` | an earlier distribution of the labelled cohort, used only for a consistency check; the pipeline runs without it | `00_build_cohort.py` (optional) |

The raw ingredients are used only to recover each protein's rescaling map and
to reconstruct conservation, normalised mobility and mean RMSF for every
scored residue. Labels are never taken from them.

## Provenance and terms

These files are redistributed here so that the analysis is reproducible from a
single clone. They are derived from the following public resources, each of
which carries its own terms. The MIT licence in the repository root covers the
source code, not these files.

| Source | What it contributes | Terms |
|---|---|---|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (NCBI) | variant identifiers and clinical classifications | public domain (US Government work); NCBI asks that the resource be cited |
| [ATLAS](https://www.dsimb.inserm.fr/ATLAS) | the molecular dynamics trajectories behind the RMSF values | CC BY 4.0 |
| AlphaMissense, via [HegeLab](https://alphamissense.hegelab.org) | `AM_Specific_Score`, `AM_Mean_Residue_Score`, `AM_Status` | **CC BY-NC-SA 4.0** (Google DeepMind): attribution, non-commercial use, share-alike |
| ADEPT (Karagöl & Karagöl, *Funct Integr Genomics* 2026) | DCS, RCS, normalised DCS and their ingredients | see <https://www.karagolresearch.com/adept> |
| UniProt | protein accessions | CC BY 4.0 |

The AlphaMissense terms are the most restrictive of these and apply to the
columns named above wherever they appear, including the copies in
`results/cohort_full.csv` and `results/AdditionalFile2_cohort.csv`. Anyone
redistributing this repository or those columns should keep the attribution
and the non-commercial, share-alike conditions.
