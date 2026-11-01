# Gene-level clustering in variant-classifier benchmarks

Analysis code for *Gene-level clustering confounds cross-validated benchmarks of
structure-informed missense variant classifiers* (Karagöl & Karagöl). It
reproduces every number, table and figure in the paper and its Additional
file 1, and provides the diagnostics as a small library you can apply to your
own cohort.

## What the code is about

Clinically labelled missense variants are distributed very unevenly across
genes, and many genes contribute variants of one class only. In the cohort
studied here, 15 of 18 genes are class-pure, so 84% of the pathogenic–benign
pairs that the AUC averages over compare variants from *different* genes. A
pooled AUC on such a cohort is therefore largely a statement about how a score
ranks genes, not variants, and the design of the cross-validation decides which
of the two it measures.

The core tool is an exact decomposition of any AUC into the concordance among
within-gene pairs and the concordance among between-gene pairs,

    AUC = w · AUC_W + (1 − w) · AUC_B,

where `w` is the fraction of pairs that lie within a gene. Around it the
repository implements gene-prevalence baselines, structure-preserving
permutation nulls, clustered inference, calibration corrections for
leave-one-gene-out designs, and an exact power calculation for the within-gene
test.

Nothing in the library is specific to genes: `groups` can be any clustering of
the rows (protein, family, study, patient).

## Install and run

```bash
git clone https://github.com/<owner>/gene-clustering-benchmark.git
cd gene-clustering-benchmark
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./run_all.sh                    # full run, roughly one hour on two cores
```

A fast end-to-end check that exercises every script in about ten minutes (Monte
Carlo quantities will differ from the paper):

```bash
DEPBENCH_NPERM=50 DEPBENCH_NBOOT=200 DEPBENCH_SIMREPS=8 DEPBENCH_DESIGNREPS=30 ./run_all.sh
```

Rebuilding the manuscript documents additionally needs Node with the `docx`
package; the analysis itself does not.

Results land in `results/` (tables, JSON summaries, out-of-fold predictions,
permutation distributions, simulation output) and `results/figures/`
(600-dpi PNG, vector PDF, 300-dpi TIFF). `src/audit.py` then recomputes the key
quantities straight from `results/` and checks them against the generated
manuscript text, failing if anything disagrees:

```bash
PYTHONPATH=src python3 src/audit.py     # 115 checks
```

## Using the diagnostics on your own cohort

```python
from depbench.metrics import (auc_decompose, within_gene_perm_p, obuchowski_test,
                              cluster_bootstrap, percentile_ci, fast_auc)
from depbench.common import gene_prevalence, gene_prevalence_logo

d = auc_decompose(y, score, gene)        # auc, auc_within, auc_between, pair counts
obs, p_one, p_two, _ = within_gene_perm_p(y, score, gene)
a1, a2, stat, p, se = obuchowski_test(y, score1, score2, gene)   # t on (groups − 1) df
boot = cluster_bootstrap(y, gene, lambda yb, idx, gb: fast_auc(yb, score[idx]))
lo, hi = percentile_ci(boot)
```

`gene_prevalence(df, leave_one_out=True)` is the variant-level lookup baseline:
if your model cannot beat it, it may be recognising genes. `gene_prevalence_logo(df)`
is the corresponding reference under leave-one-gene-out.

**Is a null result informative?** Because the within-gene permutation null
depends only on each group's class counts, power can be computed exactly for
your cohort before you interpret anything:

```python
from depbench.power import within_null, power_curve, min_detectable_within_auc

null  = within_null(y, gene)                              # computed once, reused
curve = power_curve(y, gene, [0.6, 0.65, 0.7, 0.75, 0.8], null)
min_detectable_within_auc(curve)      # smallest AUC_W reaching 80% power
```

**Is a below-chance leave-one-gene-out AUC real?** Under that design each
between-group pair is ranked by two different models, trained on different
class balances. Both corrections below act within a fold and are strictly
increasing, so neither can change `AUC_W`:

```python
from depbench.calibration import calibration_fit, prior_shift, leave_two_groups_out_auc

calibration_fit(y, p)                 # recalibration intercept, slope, Brier
q = prior_shift(p, folds, y)          # remove each fold's training prior
auc_b, n = leave_two_groups_out_auc(fit_predict, X, y, gene)   # no cross-fold comparison
```

## Pipeline

| Script | Produces |
|---|---|
| `00_build_cohort.py` | labelled cohort; recovery of each protein's rescaling map and reconstruction of conservation, normalised mobility and mean RMSF for every residue; label-free within-protein encodings |
| `01_ablation.py` | 2⁴ logistic ablation (split × scaling × selection × tuning) and 2³ ablation for three model families, each with the within/between decomposition (Table 2, Tables S5–S6) |
| `02_diagnostics.py` | gene–class association, random-intercept heterogeneity, pair composition, gene-prevalence baselines, permutation nulls (Table 1, Table S4) |
| `03_normalisation.py` | rescaling provenance, gene-identity recovery, transfer of alternative encodings, increment over AlphaMissense (Tables S7–S10) |
| `03b_rank_convention.py` | mid-rank versus rank/n percentile conventions |
| `05_ensemble_audit.py` | exploratory ensemble: recovered weights, nested blending (Table S12) |
| `06_within_gene.py` | within-gene discrimination, stratified permutation and bootstrap, conditional logistic regression (Tables S3, S11) |
| `04_simulation.py` | simulation with known estimands (Table S13) |
| `07_design.py` | genes required for a stable estimate (Table S14) |
| `08_tables.py` | Tables S1–S2 |
| `11_calibration.py` | calibration of fold-specific models, prior-shift and fold-rank corrections, leave-two-genes-out between-gene AUC (Tables S15–S16) |
| `12_robustness.py` | gene-drop influence, power of the within-gene test, paired within-gene comparisons (Tables S17–S18) |
| `09_figures.py`, `10_supp_figures.py`, `13_figure6.py` | Figures 1–6 and S1–S4 |
| `audit.py` | recomputes key quantities from `results/` and checks them against the manuscript; checks citation order |

The library under `src/depbench/` holds the reusable parts: `metrics` (AUC
decomposition, clustered and permutation inference), `common` (cohort
construction, cross-validation, baselines), `calibration`, `power`, `glmm`
(random-intercept model by dense-grid integration), `fastlr` (Newton solver
used inside permutation loops) and `style` (figure theme).

## Design notes

* `common.cross_validate` implements every configuration in the ablation;
  `FoldCache` caches label-independent preprocessing per split, and
  `engine="fast"` solves the identical penalised logistic objective by Newton's
  method inside permutation loops.
* Monte Carlo sizes come from `DEPBENCH_NPERM` (default 2000), `DEPBENCH_NBOOT`
  (5000), `DEPBENCH_SIMREPS` (150) and `DEPBENCH_DESIGNREPS` (300). Seeds are
  fixed in `depbench.common.SEED`.
* Paired AUC comparisons use Obuchowski's clustered variance with the statistic
  referred to a *t* distribution on (number of groups − 1) degrees of freedom;
  with few clusters the normal reference is anticonservative.
* The random-intercept model integrates each gene's random effect on a dense
  uniform grid. Gauss–Hermite rules, fixed or adaptive, are inaccurate when the
  between-group SD is large and most groups are class-pure.
* The reconstruction of score components uses residues of uncertain
  significance only through their unlabelled ingredients.

## Data

See [`data/README.md`](data/README.md) for the contents of each input file and
for provenance and licensing. ClinVar, ATLAS and UniProt are permissively
licensed; the AlphaMissense columns are CC BY-NC-SA 4.0 and that condition
travels with them. Dynamics-aware scores can be regenerated with
[ADEPT](https://www.karagolresearch.com/adept) as described in Karagöl &
Karagöl, *Funct Integr Genomics* 2026;26:228.

## Citing

See `CITATION.cff`. Please cite the article; the journal and DOI will be added
once it is published.

## Licence

MIT for the code (see `LICENSE`). The input data files keep the terms of their
original sources, listed in `data/README.md`.
