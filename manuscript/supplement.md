Additional file 1

Supplementary methods, tables and figures for: *Gene-level clustering confounds cross-validated benchmarks of structure-informed missense variant classifiers*

Taner Karagöl (McGovern Institute for Brain Research, MIT; karagolt@mit.edu) and Alper Karagöl (Department of Mechanical Engineering, MIT; akaragol@mit.edu)

All analyses use the same cohort: 88 ClinVar-classified missense variants in 18 genes (36 pathogenic or likely pathogenic, 52 benign or likely benign; base rate 0.409). Every value is produced by the scripts in Additional file 3; file names in brackets refer to its results directory.

## Contents

Supplementary Methods S1–S10

Table S1. Discrimination with gene-cluster and naive bootstrap intervals

Table S2. Paired comparisons: DeLong and clustered tests

Table S3. Per-gene composition, multipliers and within-gene discrimination

Table S4. Permutation null distributions

Table S5. Validation ablation for three model families

Table S6. Factorial effects for logistic regression

Table S7. Per-protein rescaling parameters

Table S8. Gene-identity recovery

Table S9. Transfer of alternative encodings

Table S10. Incremental value over AlphaMissense

Table S11. Within-gene discrimination of raw scores

Table S12. Audit of the exploratory ensemble

Table S13. Simulation summary

Table S14. Design requirements

Table S15. Calibration of out-of-fold predictions

Table S16. Fold-wise corrections and leave-two-genes-out

Table S17. Influence of each gene

Table S18. What the cohort can detect within genes

Table S19. DOME checklist

Figures S1–S4

## Supplementary Methods

### S1. Within-gene and between-gene decomposition of the AUC

Let *P* and *B* be the sets of pathogenic and benign variants, *s* a score and ψ(*a*, *b*) = 1 if *s*_{a} > *s*_{b}, 1/2 if *s*_{a} = *s*_{b} and 0 otherwise. The Mann–Whitney form of the AUC is the mean of ψ over all |*P*|·|*B*| pairs. Partition the pairs into *W*, in which both variants belong to the same gene, and its complement. Then

AUC = (1/|*P*||*B*|) Σ_{(a,b)} ψ(a,b) = *w*·AUC_{W} + (1 − *w*)·AUC_{B},

with *w* = |*W*|/(|*P*||*B*|), AUC_{W} the mean of ψ over *W* and AUC_{B} the mean over the complement. The decomposition is exact for any score, including out-of-fold predictions. A score that is constant within genes has AUC_{W} = 1/2, so its AUC is 1/2 + (1 − *w*)(AUC_{B} − 1/2), and *w* bounds how far within-gene information can move the pooled AUC. AUC_{W} equals the pair-weighted mean of the per-gene AUCs, weighted by *n*_{P,g}·*n*_{B,g}, and this equality is exact. It is therefore a within-covariate (stratified) concordance in the spirit of the covariate-adjusted ROC of Janes and Pepe (Am J Epidemiol 2008) with gene as the covariate, differing in that genes enter in proportion to their pairs rather than their cases. Under leave-one-gene-out validation, both members of a within-gene pair are predicted by the model fitted without their gene, so AUC_{W} is a fold-internal statistic; every between-gene pair compares predictions of two different fitted models.

Within-gene inference holds every gene's class counts fixed. The null of AUC_{W} is obtained by permuting labels within each mixed gene (20,000 Monte Carlo replicates, which approximate the exact stratified permutation null); intervals resample variants within gene and class. Because only genes with both classes contribute, AUC_{W} is undefined when no such gene is present.

### S2. Reconstruction of the position-level scores

For residue *i* of protein *k*, with conservation *c*_{i} and mean backbone RMSF *r*_{i} (three ATLAS replicas), the published quantities are

*m*_{i} = (*r*_{i} − min_{k} *r*)/(max_{k} *r* − min_{k} *r*),   DCS_{i} = 10·(*c*_{i}*m*_{i} − min_{k} *cm*)/(max_{k} *cm* − min_{k} *cm*),

RCS_{i} = 10·(*c*_{i}(1 − *m*_{i}) − min_{k})/(max_{k} − min_{k}),   nDCS_{i} = log(1 + DCS_{i})/log 11,

where minima and maxima run over all residues of protein *k*. Each published score is thus an affine function of a raw product with protein-specific slope and intercept. For residues with published ingredients (variants of uncertain significance at other positions of the same proteins) the identities Raw DCS = *c*·*m* and nDCS = log(1 + DCS)/log 11 hold to machine precision. Linear regression of each published score on its raw product recovers the slope and intercept for each of 26 proteins with a maximum residual of 6.0 × 10^{−10} score units. Inverting both maps gives *c*·*m* and *c*(1 − *m*) for every scored residue; their sum is *c* and their ratio gives *m*, and the per-protein linear map between *m* and *r* gives mean RMSF. On the 40 labelled residues whose ingredients are available directly, the maximum absolute error was 8.8 × 10^{−10}. The slope of the DCS map (score units per unit of *c*·*m*) is reported as the protein's DCS multiplier (Table S7).

### S3. Clustered paired comparison of AUCs

Obuchowski's estimator treats genes *j* = 1, …, *I* as independent units. For each score, placement values *V*_{10}(*a*) = mean_{b∈B} ψ(*a*, *b*) and *V*_{01}(*b*) = mean_{a∈P} ψ(*a*, *b*) are summed within genes to *V*_{10,j} and *V*_{01,j}. With *m*_{j} and *n*_{j} the numbers of positives and negatives in gene *j*, *M* = Σ*m*_{j}, *N* = Σ*n*_{j} and θ the AUC,

S_{10} = I/((I − 1)M) Σ_{j}(V_{10,j} − m_{j}θ)^{2},  S_{01} = I/((I − 1)N) Σ_{j}(V_{01,j} − n_{j}θ)^{2},  S_{11} = I/(I − 1) Σ_{j}(V_{10,j} − m_{j}θ)(V_{01,j} − n_{j}θ),

and Var(θ) = S_{10}/M + S_{01}/N + 2S_{11}/(MN). Covariances between two scores use the corresponding cross-products. Because the cohort has only *I* = 18 genes, the standardised difference in AUC is referred to a *t* distribution on *I* − 1 degrees of freedom rather than to the standard normal; with few clusters the variance is itself estimated imprecisely and the normal reference is anticonservative. When every gene contains one variant the estimator reduces to the DeLong variance up to the factor I/(I − 1). Intervals for differences also come from a paired gene-cluster bootstrap that resamples genes with replacement (5,000 replicates), holding out-of-fold predictions fixed.

### S4. Random-intercept model

logit P(*y*_{ij} = 1) = α + *u*_{j}, *u*_{j} ~ N(0, σ_{g}^{2}), was fitted by maximum likelihood. The random effect of each gene was integrated out by the trapezoidal rule on a uniform grid of 96,001 points spanning ±12 standard deviations; a fixed Gauss–Hermite rule, adaptive or not, is inaccurate here, because when σ_{g} is large the integrand of a class-pure gene is a normal density truncated sharply where that gene's log-odds change sign, and the resulting quadrature error inflates the maximised likelihood. The grid result was confirmed against a grid twice as fine and twice as wide (log-likelihoods agreeing to eight decimals). A profile-likelihood 95% interval for σ_{g} was obtained by maximising over α at fixed σ_{g} and locating where twice the drop in log-likelihood reaches 3.84. The maximum-likelihood estimate is 4.3 (95% interval 1.6–16.2); the latent-scale intraclass correlation σ_{g}^{2}/(σ_{g}^{2} + π^{2}/3) is 0.85 (0.44–0.99). The upper limit is wide because fifteen of the eighteen genes are class-pure.

### S5. Simulation and design calculations

For each cell of the simulation grid, the intercept *a* was set by bisection so that the marginal pathogenic rate equalled 0.41. Each replicate drew a cohort, computed the cross-validated estimators (leave-one-variant-out and leave-one-gene-out with pooled predictions) for the three predictors, refitted each predictor on the complete simulated cohort and applied it to (i) 40 new variants in each existing gene, with that gene's *u*_{j} and *f*_{j}, for the seen-gene AUC and (ii) 400 new genes with 10 variants each for the new-gene AUC and its within-gene component. Logistic models were fitted by Newton's method with the scikit-learn objective (L2 penalty, C = 1, balanced class weights) after standardisation. The nearest-neighbour learner used the seven nearest training variants in standardised (*x̃*, *h̃*) space with inverse-distance weights; leave-one-out excluded the variant itself, and leave-one-gene-out excluded all variants of its gene. Cells were replicated 150 times.

For the design calculations, the within-gene effect β was found by bisection so that the population within-gene AUC of the score *x* equalled 0.65 or 0.75. For each number of genes, 300 cohorts were drawn and the sampling standard deviations of the pooled and within-gene AUCs of *x* were computed. Half-widths 1.96 × SD were regressed on the number of genes on log–log axes (slopes close to −1/2) and solved for ±0.05 and ±0.03.

### S6. Exploratory ensemble

Input features (twelve): DCS, RCS and mean RMSF (fixed), plus RMSF of each of the three replicas, residue-level N_{eq} of the second replica, |Δhydrophobicity|, |Δcharge|, |Δvolume|, alternate-residue hydrophobicity and an alternate-residue glycine/proline indicator, selected by ANOVA F on the complete labelled set after imputation and Yeo–Johnson transformation fitted on the complete labelled set. Learners and tuning: XGBoost, LightGBM, CatBoost, random forest and logistic regression, each tuned by 75 Optuna trials against repeated stratified 5 × 5-fold AUC on the complete labelled set. Evaluation: leave-one-variant-out predictions for each tuned learner. Blending: weights chosen by 300 Optuna trials maximising the leave-one-variant-out AUC of the weighted average. Four stages therefore had access to the held-out variant: preprocessing, selection, tuning and blending; the variant-level split additionally exposed the held-out variant's gene. The nested re-estimation of the blend is described in the main text.

### S7. Calibration and cross-fold estimation

Let p_{ij} be the out-of-fold prediction for variant *i* of gene *j*.  The calibration intercept and slope are the coefficients of logit P(y = 1) = a + b·logit(p), fitted by unpenalised logistic regression; b = 1 with a = 0 denotes perfect calibration, b < 1 an overconfident predictor and b < 0 a reversed one.  Under leave-one-gene-out validation the model that predicts gene *j* is trained on a set whose pathogenic fraction π_{−j} excludes gene *j*, so π_{−j} moves away from gene *j*'s own fraction whenever that gene is class-pure.  The prior-shift correction removes this by setting

logit(p′_{ij}) = logit(p_{ij}) − logit(π_{−j}) + logit(π),

with π the cohort fraction.  The fold-rank correction instead maps each prediction to its percentile among the predictions the same model makes on its own training set.  Both act within a fold and are strictly increasing, so AUC_{W}, which compares only variants of the same gene, is numerically unchanged; this is verified in Table S16.  A predictor that scores every variant by π_{−j} is constant within a fold, so its AUC is produced entirely by differences between folds, and the prior shift must return it to exactly 0.5.  Leave-two-genes-out fits a model on the complement of each unordered pair of genes and scores both held-out genes with it, so that each between-gene pair is ranked by one model that saw neither gene; summing over all pairs covers every between-gene pair exactly once.

### S8. Power of the within-gene permutation test

Within-gene permutation holds each gene's class counts fixed, so for a score with no ties the permutation distribution of AUC_{W} is determined by the counts alone and not by the score: the test is distribution-free given the cohort structure.  The null may therefore be computed once and reused.  To obtain power at a population within-gene AUC of τ, scores were drawn as *s* ~ N(0, 1) with √2·Φ^{−1}(τ) added for pathogenic variants, which gives P(*s*_{P} > *s*_{B}) = τ within every gene, and AUC_{W} was compared with the two-sided 5% critical values of the stored null.  Differences between two scores were summarised by a stratified bootstrap resampling variants within gene and class, which preserves the pairing that makes the two within-gene AUCs comparable.

### S9. Software and reproducibility

Python 3.13; scikit-learn 1.9.1; SciPy 1.18.1; statsmodels 0.15.0; NumPy 2.5.3; pandas 3.0.5; Matplotlib 3.11.2. Random seeds are fixed. `run_all.sh` regenerates every result, table and figure, and `audit.py` checks every number in the manuscript against the results directory and the citation list against the reference list. Logistic regression in permutation loops uses an exact Newton solver for the scikit-learn objective, with label-independent preprocessing cached per split; observed values computed both ways agree to three decimals (C1 0.680, C8 0.594, position-level scores 0.152). Reporting follows the DOME recommendations (Table S19).

### S10. Contents of Additional file 2

**Columns of Additional file 2.**

Column(s) | Content

gene, uniprot_accession, atlas_protein_id | gene symbol, UniProt accession and ATLAS identifier

clinvar_variation_id | ClinVar variation (VCV) accession

protein_change, residue_position, ref_aa, alt_aa | HGVS protein change and parsed residues

clinvar_classification, pathogenic | ClinVar classification and binary label (1 = P/LP)

n_variants_in_gene, gene_pathogenic_fraction | gene size and pathogenic fraction in the cohort

DCS_Score, RCS_Score, Normalized_DCS | published position-level scores (per-protein 0-10)

conservation_reconstructed, norm_dynamics_reconstructed, mean_rmsf_reconstructed_A | components recovered by inverting each protein's rescaling map (Methods S2)

dcs_protein_scale_factor | DCS multiplier of the protein (score units per unit of conservation x mobility)

DCS_within_protein_percentile, RCS_within_protein_percentile | mid-rank percentile among all ClinVar-observed residues of the protein (label-free)

AM_variant_score, AM_residue_mean_score | AlphaMissense variant score and residue mean

dHydro ... hinge_disr | 16 substitution descriptors (Methods)

**Table S1. Discrimination with gene-cluster and naive bootstrap intervals.** 5,000 resamples each. The cluster bootstrap resamples genes; the naive bootstrap resamples variants. Ratio: width of the cluster interval over the naive interval. Obuchowski and DeLong SE: analytic standard errors treating genes or variants as independent. [TableS1_discrimination.csv]

Score or model | AUC | AP | AUC_{W} | AUC_{B} | Cluster 95% CI | Naive 95% CI | Ratio | SE clustered | SE DeLong

AlphaMissense (variant score) | 0.919 | 0.870 | 0.896 | 0.924 | 0.858–0.989 | 0.855–0.970 | 1.13 | 0.029 | 0.029

AlphaMissense (residue mean) | 0.872 | 0.810 | 0.772 | 0.891 | 0.770–0.975 | 0.787–0.940 | 1.34 | 0.049 | 0.039

DCS | 0.538 | 0.444 | 0.681 | 0.511 | 0.244–0.785 | 0.408–0.664 | 2.12 | 0.143 | 0.066

RCS | 0.518 | 0.488 | 0.668 | 0.490 | 0.260–0.734 | 0.389–0.646 | 1.84 | 0.111 | 0.066

Normalised DCS | 0.538 | 0.444 | 0.681 | 0.511 | 0.244–0.785 | 0.408–0.664 | 2.12 | 0.143 | 0.066

Conservation (reconstructed) | 0.525 | 0.461 | 0.668 | 0.498 | 0.263–0.763 | 0.398–0.653 | 1.96 | 0.118 | 0.065

Normalised mobility (reconstructed) | 0.483 | 0.409 | 0.483 | 0.483 | 0.315–0.639 | 0.364–0.607 | 1.33 | 0.070 | 0.063

|Δ hydrophobicity| | 0.688 | 0.570 | 0.673 | 0.691 | 0.554–0.798 | 0.567–0.803 | 1.03 | 0.042 | 0.059

|Δ volume| | 0.629 | 0.557 | 0.565 | 0.642 | 0.513–0.789 | 0.506–0.746 | 1.15 | 0.054 | 0.062

|Δ charge| | 0.663 | 0.531 | 0.619 | 0.672 | 0.573–0.782 | 0.564–0.764 | 1.05 | 0.045 | 0.052

Gene prevalence (variant excluded) | 0.739 | 0.722 | 0.000 | 0.879 | 0.340–0.980 | 0.629–0.841 | 3.02 | 0.217 | 0.054

Exploratory ensemble (variant-level, as reported) | 0.756 | 0.658 | 0.718 | 0.763 | 0.547–0.865 | 0.648–0.856 | 1.53 | 0.062 | 0.053

Full 19-feature model (C8) † | 0.594 | 0.527 | 0.550 | 0.602 | 0.501–0.730 | 0.463–0.719 | 0.90 | 0.042 | 0.065

Position-level scores only † | 0.152 | 0.267 | 0.289 | 0.126 | 0.028–0.267 | 0.075–0.242 | 1.43 | 0.067 | 0.042

Substitution descriptors only † | 0.635 | 0.546 | 0.607 | 0.640 | 0.520–0.759 | 0.506–0.753 | 0.97 | 0.045 | 0.062

AlphaMissense only † | 0.874 | 0.816 | 0.896 | 0.870 | 0.776–0.978 | 0.787–0.946 | 1.27 | 0.036 | 0.040

AlphaMissense + position-level scores † | 0.854 | 0.807 | 0.889 | 0.848 | 0.703–0.968 | 0.757–0.935 | 1.49 | 0.044 | 0.045

AlphaMissense + full feature set † | 0.873 | 0.832 | 0.893 | 0.869 | 0.768–0.962 | 0.789–0.944 | 1.25 | 0.037 | 0.039

† Gene-held-out logistic model (out-of-fold predictions). Ratios exceed 1 for 16 of 18 rows; the largest belong to scores whose pooled AUC rests almost entirely on between-gene pairs, where the clustering is, so that treating variants as independent understates the uncertainty most.

**Table S2. Paired comparisons of AUCs.** DeLong treats variants as independent; the clustered test (Obuchowski) treats genes as independent. [TableS2_paired_tests.csv]

Comparison | AUC 1 | AUC 2 | Difference | DeLong z | DeLong p | Clustered z | Clustered p

AlphaMissense (variant score) vs Exploratory ensemble (variant-level, as reported) | 0.919 | 0.756 | +0.163 | 3.44 | 5.9 × 10^{−4} | 3.40 | 0.00344

Exploratory ensemble (variant-level, as reported) vs Gene prevalence (variant excluded) | 0.756 | 0.739 | +0.017 | 0.25 | 0.803 | 0.08 | 0.934

AlphaMissense (variant score) vs Full 19-feature model (C8) † | 0.919 | 0.594 | +0.325 | 4.65 | 3.3 × 10^{−6} | 8.19 | 2.7 × 10^{−7}

AlphaMissense (variant score) vs AlphaMissense only † | 0.919 | 0.874 | +0.045 | 2.44 | 0.0148 | 2.84 | 0.0113

AlphaMissense only † vs AlphaMissense + position-level scores † | 0.874 | 0.854 | +0.020 | 1.98 | 0.0474 | 1.59 | 0.129

AlphaMissense only † vs AlphaMissense + full feature set † | 0.874 | 0.873 | +0.001 | 0.08 | 0.94 | 0.07 | 0.944

Full 19-feature model (C8) † vs Substitution descriptors only † | 0.594 | 0.635 | -0.041 | -1.79 | 0.073 | -3.10 | 0.00646

Full 19-feature model (C8) † vs Position-level scores only † | 0.594 | 0.152 | +0.442 | 6.05 | 1.5 × 10^{−9} | 4.46 | 3.4 × 10^{−4}

Substitution descriptors only † vs Position-level scores only † | 0.635 | 0.152 | +0.483 | 6.33 | 2.4 × 10^{−10} | 4.88 | 1.4 × 10^{−4}

DCS vs Conservation (reconstructed) | 0.538 | 0.525 | +0.013 | 0.20 | 0.838 | 0.16 | 0.871

DCS vs RCS | 0.538 | 0.518 | +0.020 | 0.26 | 0.795 | 0.22 | 0.826

AlphaMissense (residue mean) vs DCS | 0.872 | 0.538 | +0.334 | 4.07 | 4.7 × 10^{−5} | 1.89 | 0.0763

AlphaMissense (variant score) vs AlphaMissense (residue mean) | 0.919 | 0.872 | +0.048 | 2.23 | 0.0259 | 1.89 | 0.0763

† Gene-held-out logistic model.

**Table S3. Per-gene composition, DCS multiplier and within-gene discrimination.** Within-gene AUCs are defined only for genes containing both classes. [TableS_per_gene.csv]

Gene | UniProt | n | P/LP | B/LB | Pathogenic fraction | Within pairs | DCS multiplier | AUC AM | AUC DCS | AUC RCS | AUC conservation

PGK1 | P00558 | 34 | 16 | 18 | 0.47 | 288 | 12.1 | 0.90 | 0.70 | 0.69 | 0.69

GBP1 | P32455 | 13 | 0 | 13 | 0.00 | 0 | 10.8 | — | — | — | —

HMGCS2 | P54868 | 9 | 9 | 0 | 1.00 | 0 | 11.2 | — | — | — | —

GHRHR | Q02643 | 5 | 2 | 3 | 0.40 | 6 | 31.7 | 1.00 | 0.00 | 0.00 | 0.00

KAT8 | Q9H7Z6 | 5 | 4 | 1 | 0.80 | 4 | 28.8 | 0.50 | 0.50 | 0.00 | 0.00

PARP3 | Q9Y6F1 | 4 | 0 | 4 | 0.00 | 0 | 26.5 | — | — | — | —

LGALS3BP | Q08380 | 3 | 0 | 3 | 0.00 | 0 | 22.5 | — | — | — | —

BLMH | Q13867 | 3 | 0 | 3 | 0.00 | 0 | 10.0 | — | — | — | —

ANGPTL6 | Q8NI99 | 2 | 0 | 2 | 0.00 | 0 | 24.6 | — | — | — | —

NARS1 | O43776 | 2 | 2 | 0 | 1.00 | 0 | 30.7 | — | — | — | —

NT5C2 | P49902 | 1 | 0 | 1 | 0.00 | 0 | 14.2 | — | — | — | —

MSN | P26038 | 1 | 1 | 0 | 1.00 | 0 | 22.1 | — | — | — | —

PPARA | Q07869 | 1 | 0 | 1 | 0.00 | 0 | 10.3 | — | — | — | —

SPRTN | Q9H040 | 1 | 1 | 0 | 1.00 | 0 | 12.8 | — | — | — | —

STAMBPL1 | Q96FJ0 | 1 | 0 | 1 | 0.00 | 0 | 21.8 | — | — | — | —

TBRG1 | Q3YBR2 | 1 | 1 | 0 | 1.00 | 0 | 11.8 | — | — | — | —

TPH1 | P17752 | 1 | 0 | 1 | 0.00 | 0 | 12.4 | — | — | — | —

TTC5 | Q8N0Z6 | 1 | 0 | 1 | 0.00 | 0 | 14.1 | — | — | — | —

Gene–class association: χ² = 44.7, df = 17, Monte Carlo p ≤ 1 × 10^{−5} (86% of expected counts below five; asymptotic p not valid); Cramér's V = 0.71.

**Table S4. Permutation null distributions.** Complete pipeline refitted in every replicate. p (upper): proportion of null AUCs at or above the observed value; p (lower): at or below; p (2-sided): twice the smaller tail; all (b + 1)/(B + 1). [permutation_summary.csv, permutation_distributions.csv]

Pipeline | Null | Replicates | Observed | Null mean | SD | 2.5% | 97.5% | p (upper) | p (lower) | p (2-sided) | Null AUC_{W} | Null AUC_{B}

Variant-level, leaked (C1) | global | 2,000 | 0.680 | 0.554 | 0.070 | 0.410 | 0.686 | 0.030 | 0.971 | 0.061 | 0.550 | 0.555

Variant-level, leaked (C1) | within gene | 2,000 | 0.680 | 0.619 | 0.064 | 0.487 | 0.736 | 0.170 | 0.833 | 0.341 | 0.535 | 0.635

Gene-level, controlled (C8) | global | 2,000 | 0.594 | 0.492 | 0.085 | 0.325 | 0.653 | 0.112 | 0.889 | 0.225 | 0.501 | 0.490

Gene-level, controlled (C8) | within gene | 2,000 | 0.594 | 0.466 | 0.096 | 0.274 | 0.650 | 0.086 | 0.915 | 0.173 | 0.502 | 0.459

Position-level scores, gene-held-out | global | 2,000 | 0.152 | 0.477 | 0.089 | 0.291 | 0.627 | 1.000 | 0.001 | 0.001 | 0.496 | 0.473

Position-level scores, gene-held-out | within gene | 2,000 | 0.152 | 0.291 | 0.086 | 0.171 | 0.480 | 1.000 | 0.001 | 0.002 | 0.498 | 0.252

Position-level scores, gene-held-out | donor protein | 1,000 | 0.152 | 0.471 | 0.121 | 0.227 | 0.682 | 1.000 | 0.001 | 0.002 | 0.499 | 0.466

Position-level scores, variant-level | global | 2,000 | 0.278 | 0.466 | 0.118 | 0.129 | 0.629 | 0.932 | 0.068 | 0.136 | 0.443 | 0.471

Position-level scores, variant-level | within gene | 2,000 | 0.278 | 0.387 | 0.124 | 0.036 | 0.554 | 0.839 | 0.162 | 0.324 | 0.437 | 0.377

**Table S5. Validation ablation for three model families (2^{3} design, fixed hyperparameters).** Logistic regression C = 1; random forest 400 trees, minimum leaf 3; gradient boosting 60 trees of depth 2. Each cell: AUC (AP) and, below, AUC_{W}/AUC_{B}. [ablation_full.csv]

Config. | Split | Scaling | Selection | Logistic regression | Random forest | Gradient boosting

C1 | variant | full set | full set | 0.680 (0.597); 0.604/0.694 | 0.648 (0.558); 0.644/0.649 | 0.615 (0.500); 0.641/0.610

C2 | variant | full set | in fold | 0.689 (0.586); 0.671/0.693 | 0.643 (0.548); 0.627/0.646 | 0.642 (0.523); 0.685/0.633

C3 | variant | in fold | full set | 0.679 (0.595); 0.604/0.694 | 0.648 (0.558); 0.644/0.649 | 0.615 (0.500); 0.641/0.610

C4 | variant | in fold | in fold | 0.667 (0.575); 0.601/0.680 | 0.648 (0.551); 0.634/0.651 | 0.628 (0.504); 0.641/0.625

C5 | gene | full set | full set | 0.652 (0.586); 0.574/0.666 | 0.623 (0.555); 0.678/0.613 | 0.540 (0.458); 0.658/0.518

C6 | gene | full set | in fold | 0.591 (0.526); 0.540/0.600 | 0.616 (0.516); 0.698/0.600 | 0.594 (0.482); 0.691/0.576

C7 | gene | in fold | full set | 0.648 (0.586); 0.571/0.663 | 0.623 (0.555); 0.678/0.613 | 0.539 (0.456); 0.658/0.516

C8 | gene | in fold | in fold | 0.594 (0.527); 0.550/0.602 | 0.607 (0.507); 0.695/0.590 | 0.586 (0.475); 0.685/0.567

**Table S6. Factorial effects on AUC for logistic regression (2^{4} design of Table 2 in the main text).** Effects are standard contrasts: mean AUC with the factor at 'sees held-out data' (variant split; full-set scaling, selection or tuning) minus mean AUC at the other level; interactions are the corresponding products. [ablation_effects.csv]

Effect | Order | AUC effect

split | 1 | +0.0435

selection | 1 | +0.0180

tuning | 1 | -0.0124

scaling | 1 | -0.0049

split x selection | 2 | -0.0276

selection x tuning | 2 | +0.0175

split x tuning | 2 | +0.0107

scaling x tuning | 2 | +0.0083

split x scaling | 2 | +0.0078

scaling x selection | 2 | +0.0039

split x scaling x selection | 3 | -0.0085

split x selection x tuning | 3 | -0.0069

scaling x selection x tuning | 3 | -0.0053

split x scaling x tuning | 3 | -0.0033

split x scaling x selection x tuning | 4 | +0.0028

**Table S7. Per-protein rescaling parameters recovered from the published scores.** DCS and RCS multipliers: score units per unit of raw product. RMSF at normalised mobility 0 and 1: the protein's minimum and maximum mean RMSF. Max m (ClinVar): the largest normalised mobility among ClinVar-observed residues. [TableS_provenance.csv]

Protein | Gene | Cohort | Reference residues | DCS multiplier | RCS multiplier | RMSF min (Å) | RMSF max (Å) | Max m (ClinVar) | Max residual

proteinh1 | LGALS3BP | yes | 14 | 22.53 | 12.70 | 0.60 | 12.17 | 0.47 | 6.0e-10

proteinh26 | TBRG1 | yes | 21 | 11.80 | 10.35 | 0.51 | 6.87 | 0.42 | 5.6e-10

proteinh56 | ANGPTL6 | yes | 26 | 24.62 | 10.42 | 0.34 | 7.85 | 0.12 | 5.7e-10

proteinh52 | SPRTN | yes | 7 | 12.77 | 11.27 | 1.18 | 12.20 | 0.28 | 5.0e-10

proteinh48 | NT5C2 | yes | 4 | 14.21 | 10.36 | 0.50 | 7.51 | 0.65 | 2.0e-10

proteinh46 | NARS1 | yes | 13 | 30.71 | 10.47 | 0.49 | 4.30 | 0.33 | 4.7e-10

proteinh38 | PARP3 | yes | 43 | 26.46 | 10.98 | 0.44 | 6.59 | 0.58 | 5.6e-10

proteinh3 | MSN | yes | 13 | 22.12 | 12.11 | 3.74 | 18.69 | 0.37 | 4.6e-10

proteinh28 | TTC5 | yes | 7 | 14.06 | 12.72 | 0.39 | 7.74 | 0.21 | 3.8e-10

proteinh27 | GHRHR | yes | 19 | 31.68 | 11.02 | 0.75 | 8.06 | 0.42 | 5.2e-10

proteinh29 | STAMBPL1 | yes | 14 | 21.83 | 11.64 | 0.42 | 10.42 | 0.16 | 6.0e-10

proteinh25 | PGK1 | yes | 170 | 12.05 | 10.86 | 0.48 | 3.14 | 0.97 | 5.4e-10

proteinh24 | HMGCS2 | yes | 35 | 11.24 | 10.90 | 0.41 | 4.58 | 0.93 | 5.4e-10

proteinh22 | KAT8 | yes | 41 | 28.82 | 10.37 | 0.51 | 9.57 | 0.28 | 5.6e-10

proteinh21 | PPARA | yes | 4 | 10.31 | 13.11 | 0.45 | 4.08 | 0.52 | 4.5e-10

proteinh18 | BLMH | yes | 38 | 10.01 | 10.44 | 0.56 | 17.69 | 0.61 | 5.0e-10

proteinh14 | TPH1 | yes | 32 | 12.44 | 10.23 | 0.46 | 9.03 | 0.18 | 5.3e-10

proteinh8 | GBP1 | yes | 68 | 10.81 | 11.37 | 0.56 | 8.02 | 0.75 | 5.2e-10

proteinh12 | SERPINE1 | no | 6 | 19.61 | 11.47 | 0.36 | 3.80 | 0.20 | 4.0e-10

proteinh33 | MTERF3 | no | 4 | 19.23 | 12.52 | 0.91 | 5.62 | 0.28 | 4.6e-10

proteinh34 | P4HA1 | no | 21 | 10.01 | 14.47 | 1.77 | 10.66 | 0.80 | 5.1e-10

proteinh36 | BECN1 | no | 4 | 25.36 | 11.34 | 0.68 | 11.04 | 0.16 | 2.8e-10

proteinh2 | CPA2 | no | 5 | 22.71 | 12.96 | 0.39 | 2.64 | 0.64 | 4.1e-10

proteinh5 | MEF2A | no | 2 | 12.20 | 11.88 | 2.95 | 19.33 | 0.20 | 8.9e-16

proteinh16 | SRC | no | 3 | 29.18 | 13.61 | 0.41 | 3.93 | 0.19 | 2.5e-10

proteinh55 | NFS1 | no | 7 | 10.64 | 12.01 | 1.30 | 18.85 | 0.23 | 4.9e-10

**Table S8. Gene-identity recovery from each encoding.** Random forest, five-fold cross-validation stratified by gene and grouped by residue, five repeats; null from 50 label permutations. [gene_identity_recovery.csv]

Set | Encoding | Genes | Balanced accuracy | Accuracy | Majority rate | Null mean | Null SD | p

Labelled cohort | DCS, RCS, normalised DCS as published (per-protein min-max, 0-10) | 8 | 0.264 | 0.560 | 0.447 | 0.127 | 0.027 | 0.020

Labelled cohort | conservation x mobility products before per-protein rescaling | 8 | 0.246 | 0.513 | 0.447 | 0.128 | 0.036 | 0.020

Labelled cohort | conservation and normalised mobility | 8 | 0.211 | 0.508 | 0.447 | 0.125 | 0.041 | 0.059

Labelled cohort | conservation alone | 8 | 0.197 | 0.300 | 0.447 | 0.120 | 0.041 | 0.078

Labelled cohort | normalised mobility alone | 8 | 0.205 | 0.340 | 0.447 | 0.117 | 0.038 | 0.020

Labelled cohort | mean backbone RMSF (absolute, Angstrom) | 8 | 0.201 | 0.371 | 0.447 | 0.117 | 0.043 | 0.059

Labelled cohort | within-protein percentile, label-free reference | 8 | 0.067 | 0.226 | 0.447 | 0.120 | 0.035 | 1.000

Labelled cohort | AlphaMissense residue mean (absolute scale) | 8 | 0.190 | 0.303 | 0.447 | 0.107 | 0.041 | 0.078

Labelled cohort | AlphaMissense residue mean, per-protein min-max | 8 | 0.182 | 0.303 | 0.447 | 0.112 | 0.038 | 0.078

Labelled cohort | AlphaMissense residue mean x published DCS multiplier | 8 | 0.259 | 0.421 | 0.447 | 0.121 | 0.040 | 0.020

Labelled cohort | AlphaMissense residue mean, within-protein percentile | 8 | 0.067 | 0.179 | 0.447 | 0.118 | 0.048 | 0.882

Labelled cohort | 16 substitution descriptors | 8 | 0.174 | 0.384 | 0.447 | 0.133 | 0.050 | 0.176

All ClinVar residues | DCS, RCS, normalised DCS as published (per-protein min-max, 0-10) | 15 | 0.144 | 0.256 | 0.273 | 0.066 | 0.010 | 0.020

All ClinVar residues | conservation x mobility products before per-protein rescaling | 15 | 0.127 | 0.268 | 0.273 | 0.067 | 0.012 | 0.020

All ClinVar residues | conservation and normalised mobility | 15 | 0.156 | 0.293 | 0.273 | 0.069 | 0.010 | 0.020

All ClinVar residues | conservation alone | 15 | 0.116 | 0.159 | 0.273 | 0.068 | 0.014 | 0.020

All ClinVar residues | normalised mobility alone | 15 | 0.105 | 0.208 | 0.273 | 0.068 | 0.013 | 0.020

All ClinVar residues | mean backbone RMSF (absolute, Angstrom) | 15 | 0.125 | 0.183 | 0.273 | 0.065 | 0.013 | 0.020

All ClinVar residues | within-protein percentile, label-free reference | 15 | 0.015 | 0.054 | 0.273 | 0.065 | 0.010 | 1.000

All ClinVar residues | AlphaMissense residue mean (absolute scale) | 15 | 0.104 | 0.187 | 0.273 | 0.066 | 0.013 | 0.020

All ClinVar residues | AlphaMissense residue mean, per-protein min-max | 15 | 0.071 | 0.145 | 0.273 | 0.066 | 0.011 | 0.333

All ClinVar residues | AlphaMissense residue mean x published DCS multiplier | 15 | 0.162 | 0.255 | 0.273 | 0.069 | 0.016 | 0.020

All ClinVar residues | AlphaMissense residue mean, within-protein percentile | 15 | 0.001 | 0.002 | 0.273 | 0.064 | 0.011 | 1.000

**Table S9. Gene-held-out and variant-level AUC of alternative encodings.** Logistic regression; preprocessing in fold. [transfer_by_encoding.csv]

Encoding | Features | Gene-held-out AUC | AP | AUC_{W} | AUC_{B} | Variant-level AUC

DCS, RCS, normalised DCS as published (per-protein min-max, 0-10) | 3 | 0.152 | 0.267 | 0.289 | 0.126 | 0.278

conservation x mobility products before per-protein rescaling | 2 | 0.209 | 0.286 | 0.272 | 0.197 | 0.376

conservation and normalised mobility | 2 | 0.222 | 0.288 | 0.302 | 0.206 | 0.381

conservation alone | 1 | 0.230 | 0.304 | 0.299 | 0.217 | 0.373

normalised mobility alone | 1 | 0.360 | 0.329 | 0.517 | 0.330 | 0.434

mean backbone RMSF (absolute, Angstrom) | 1 | 0.583 | 0.439 | 0.517 | 0.595 | 0.607

within-protein percentile, label-free reference | 3 | 0.254 | 0.293 | 0.262 | 0.253 | 0.340

within-gene percentile among labelled variants | 3 | 0.336 | 0.324 | 0.248 | 0.353 | 0.442

within-gene z-score among labelled variants | 3 | 0.342 | 0.352 | 0.275 | 0.355 | 0.442

within-gene centring among labelled variants | 3 | 0.289 | 0.332 | 0.262 | 0.294 | 0.437

AlphaMissense residue mean (absolute scale) | 1 | 0.845 | 0.765 | 0.772 | 0.858 | 0.855

AlphaMissense residue mean, per-protein min-max | 1 | 0.847 | 0.768 | 0.772 | 0.862 | 0.855

AlphaMissense residue mean x published DCS multiplier | 1 | 0.858 | 0.765 | 0.772 | 0.875 | 0.855

AlphaMissense residue mean, within-protein percentile | 1 | 0.792 | 0.758 | 0.772 | 0.795 | 0.791

16 substitution descriptors | 16 | 0.635 | 0.546 | 0.607 | 0.640 | 0.667

**Table S10. Incremental value over AlphaMissense (gene-held-out logistic regression).** Δ: AUC minus that of AlphaMissense alone. CI: paired gene-cluster bootstrap. [incremental_over_am.csv]

Model | Selection | AUC | AP | AUC_{W} | AUC_{B} | Δ | Δ 95% CI | DeLong p | Clustered p

AM alone | none | 0.874 | 0.816 | 0.896 | 0.870 | — | — | — | —

AM + position scores (published) | none | 0.854 | 0.807 | 0.889 | 0.848 | -0.020 | -0.100 to -0.002 | 0.047 | 0.129

AM + position scores (unrescaled) | none | 0.857 | 0.816 | 0.909 | 0.848 | -0.017 | -0.096 to +0.004 | 0.089 | 0.259

AM + conservation | none | 0.856 | 0.800 | 0.866 | 0.855 | -0.018 | -0.100 to -0.004 | 0.054 | 0.093

AM + position scores (within-protein percentile) | none | 0.855 | 0.736 | 0.778 | 0.869 | -0.019 | -0.073 to +0.022 | 0.220 | 0.464

AM + position scores (within-gene percentile, labelled) | none | 0.852 | 0.740 | 0.772 | 0.867 | -0.022 | -0.081 to +0.007 | 0.129 | 0.371

AM + position scores (within-gene z, labelled) | none | 0.865 | 0.769 | 0.819 | 0.874 | -0.009 | -0.049 to +0.013 | 0.525 | 0.588

AM + position scores (within-gene centred, labelled) | none | 0.846 | 0.708 | 0.775 | 0.860 | -0.028 | -0.084 to +0.012 | 0.122 | 0.298

AM + substitution descriptors | k=9 | 0.886 | 0.843 | 0.886 | 0.886 | +0.012 | -0.052 to +0.043 | 0.410 | 0.405

AM + full feature set | k=9 | 0.873 | 0.832 | 0.893 | 0.869 | -0.001 | -0.061 to +0.033 | 0.940 | 0.944

**Table S11. Within-gene discrimination of raw scores.** AUC_{W}: concordance over the 298 within-gene pairs; CI from a stratified bootstrap; p from stratified permutation, two-sided on |AUC_{W} − 0.5| (20,000 Monte Carlo replicates). CLR: conditional logistic regression stratified by gene, odds ratio per cohort SD. [within_gene_discrimination.csv]

Score | AUC | AUC_{B} | AUC_{W} | AUC_{W} 95% CI | p | PGK1 | GHRHR | KAT8 | CLR OR (95% CI) | CLR p

AlphaMissense (variant) | 0.919 | 0.924 | 0.896 | 0.772–0.983 | < 0.0001 | 0.90 | 1.00 | 0.50 | 6.47 (2.16–19.39) | 0.001

AlphaMissense (residue mean) | 0.872 | 0.891 | 0.772 | 0.601–0.906 | 0.0042 | 0.78 | 1.00 | 0.00 | 3.90 (1.53–9.91) | 0.004

DCS | 0.538 | 0.511 | 0.681 | 0.497–0.849 | 0.0652 | 0.70 | 0.00 | 0.50 | 0.84 (0.37–1.93) | 0.690

RCS | 0.518 | 0.490 | 0.668 | 0.480–0.839 | 0.0883 | 0.69 | 0.00 | 0.00 | 1.15 (0.64–2.06) | 0.647

Normalised DCS | 0.538 | 0.511 | 0.681 | 0.490–0.849 | 0.0652 | 0.70 | 0.00 | 0.50 | 1.04 (0.51–2.12) | 0.923

Raw DCS (before rescaling) | 0.527 | 0.498 | 0.681 | 0.490–0.849 | 0.0652 | 0.70 | 0.00 | 0.50 | 1.00 (0.49–2.06) | 0.998

Raw RCS (before rescaling) | 0.536 | 0.511 | 0.668 | 0.483–0.842 | 0.0883 | 0.69 | 0.00 | 0.00 | 1.13 (0.63–2.04) | 0.675

Conservation | 0.525 | 0.498 | 0.668 | 0.483–0.832 | 0.0880 | 0.69 | 0.00 | 0.00 | 1.13 (0.60–2.13) | 0.699

Normalised mobility | 0.483 | 0.483 | 0.483 | 0.292–0.674 | 0.8744 | 0.49 | 0.00 | 1.00 | 0.93 (0.50–1.71) | 0.807

Mean RMSF | 0.376 | 0.355 | 0.483 | 0.285–0.678 | 0.8744 | 0.49 | 0.00 | 1.00 | 0.91 (0.33–2.46) | 0.847

|Δ hydrophobicity| | 0.688 | 0.691 | 0.673 | 0.488–0.842 | 0.0745 | 0.69 | 0.50 | 0.00 | 1.58 (0.85–2.93) | 0.150

|Δ volume| | 0.629 | 0.642 | 0.565 | 0.378–0.745 | 0.4944 | 0.58 | 0.00 | 0.00 | 1.02 (0.57–1.81) | 0.945

|Δ charge| | 0.663 | 0.672 | 0.619 | 0.475–0.763 | 0.1298 | 0.61 | 1.00 | 0.88 | 2.28 (1.09–4.76) | 0.028

**Table S12. Audit of the exploratory ensemble.** [ensemble_audit.json]

Quantity | Value

Reported (in-sample blend) | 0.756

Best of 3,006 candidate weightings, in sample | 0.753

Equal weights | 0.723

Blend with weights nested in the leave-one-out loop | 0.732

Single best learner, nested choice | 0.685

Base learner: XGBoost | 0.703

Base learner: LightGBM | 0.731

Base learner: CatBoost | 0.701

Base learner: RandomForest | 0.725

Base learner: Logistic | 0.724

Recovered weights | XGBoost 5.8%, LightGBM 30.2%, CatBoost 1.2%, RandomForest 5.7%, Logistic 57.0%

Within-gene / between-gene AUC (reported) | 0.718 / 0.763

Gene-cluster 95% CI (reported) | 0.547–0.865

Versus gene-prevalence lookup: difference, DeLong p, clustered p | +0.017, 0.80, 0.93

Spearman correlation with gene prevalence (in sample) | 0.40

**Table S13. Simulation summary (within-gene effect β = 0.75, feature offset γ = 1).** Means over 150 replicates. kNN: seven-nearest-neighbour learner on two descriptors. Bias columns: estimator minus estimand. [simulation.csv; full grid including β = 0 and γ = 0 in the file]

Design | Gene SD | kNN variant-LOO | kNN seen-gene | kNN new-gene | Variant-LOO − new | Gene-LOO − new | Gene-LOO (W) − within | LR variant-LOO − new | Lookup variant-LOO | Lookup gene-LOO

12 genes × 10 | 0.0 | 0.560 | 0.572 | 0.565 | -0.005 | -0.009 | +0.001 | -0.003 | 0.462 | 0.335

12 genes × 10 | 0.5 | 0.567 | 0.582 | 0.554 | +0.013 | -0.019 | -0.004 | -0.006 | 0.527 | 0.302

12 genes × 10 | 1.0 | 0.593 | 0.611 | 0.550 | +0.043 | -0.012 | +0.007 | +0.002 | 0.621 | 0.243

12 genes × 10 | 1.5 | 0.614 | 0.626 | 0.533 | +0.081 | -0.024 | +0.007 | +0.003 | 0.718 | 0.180

12 genes × 10 | 2.0 | 0.627 | 0.643 | 0.524 | +0.103 | -0.042 | -0.006 | +0.004 | 0.773 | 0.144

12 genes × 10 | 3.0 | 0.666 | 0.676 | 0.515 | +0.151 | -0.045 | +0.015 | +0.015 | 0.847 | 0.093

12 genes × 10 | 4.0 | 0.684 | 0.700 | 0.510 | +0.174 | -0.059 | +0.008 | +0.024 | 0.899 | 0.058

12 genes × 2 | 0.0 | 0.481 | 0.553 | 0.549 | -0.068 | -0.076 | +0.044 | -0.052 | 0.489 | 0.130

12 genes × 2 | 0.5 | 0.480 | 0.565 | 0.549 | -0.069 | -0.068 | +0.008 | -0.041 | 0.489 | 0.129

12 genes × 2 | 1.0 | 0.486 | 0.578 | 0.533 | -0.047 | -0.075 | +0.013 | -0.060 | 0.557 | 0.102

12 genes × 2 | 1.5 | 0.488 | 0.590 | 0.524 | -0.035 | -0.082 | +0.026 | -0.078 | 0.601 | 0.082

12 genes × 2 | 2.0 | 0.455 | 0.600 | 0.520 | -0.065 | -0.132 | -0.009 | -0.067 | 0.659 | 0.060

12 genes × 2 | 3.0 | 0.518 | 0.626 | 0.513 | +0.004 | -0.097 | +0.016 | -0.033 | 0.731 | 0.041

12 genes × 2 | 4.0 | 0.513 | 0.647 | 0.508 | +0.005 | -0.112 | -0.006 | -0.042 | 0.783 | 0.027

12 genes × 5 | 0.0 | 0.547 | 0.565 | 0.558 | -0.012 | -0.010 | +0.003 | -0.016 | 0.466 | 0.266

12 genes × 5 | 0.5 | 0.526 | 0.572 | 0.551 | -0.025 | -0.044 | -0.004 | -0.038 | 0.508 | 0.244

12 genes × 5 | 1.0 | 0.547 | 0.594 | 0.543 | +0.004 | -0.057 | +0.018 | -0.026 | 0.616 | 0.190

12 genes × 5 | 1.5 | 0.570 | 0.612 | 0.530 | +0.040 | -0.051 | -0.009 | -0.005 | 0.676 | 0.155

12 genes × 5 | 2.0 | 0.592 | 0.644 | 0.527 | +0.065 | -0.052 | -0.006 | +0.014 | 0.738 | 0.124

12 genes × 5 | 3.0 | 0.633 | 0.670 | 0.512 | +0.121 | -0.057 | -0.008 | +0.017 | 0.845 | 0.065

12 genes × 5 | 4.0 | 0.622 | 0.674 | 0.507 | +0.115 | -0.085 | -0.019 | +0.029 | 0.870 | 0.052

48 genes × 5 | 0.0 | 0.558 | 0.563 | 0.561 | -0.003 | -0.005 | +0.002 | -0.011 | 0.487 | 0.256

48 genes × 5 | 0.5 | 0.561 | 0.567 | 0.558 | +0.003 | -0.007 | +0.006 | -0.005 | 0.535 | 0.233

48 genes × 5 | 1.0 | 0.555 | 0.569 | 0.545 | +0.010 | -0.010 | -0.003 | -0.006 | 0.623 | 0.187

48 genes × 5 | 1.5 | 0.556 | 0.574 | 0.535 | +0.021 | -0.014 | -0.004 | -0.007 | 0.710 | 0.140

48 genes × 5 | 2.0 | 0.559 | 0.577 | 0.528 | +0.031 | -0.016 | -0.002 | -0.007 | 0.770 | 0.107

48 genes × 5 | 3.0 | 0.557 | 0.582 | 0.517 | +0.040 | -0.030 | +0.001 | -0.024 | 0.849 | 0.064

48 genes × 5 | 4.0 | 0.576 | 0.593 | 0.509 | +0.067 | -0.016 | -0.003 | -0.011 | 0.896 | 0.040

cohort-matched (18 genes) | 0.0 | 0.556 | 0.567 | 0.566 | -0.010 | -0.016 | -0.003 | -0.002 | 0.382 | 0.356

cohort-matched (18 genes) | 0.5 | 0.556 | 0.564 | 0.558 | -0.002 | -0.032 | -0.006 | -0.017 | 0.439 | 0.336

cohort-matched (18 genes) | 1.0 | 0.582 | 0.574 | 0.542 | +0.040 | -0.028 | -0.003 | -0.011 | 0.536 | 0.287

cohort-matched (18 genes) | 1.5 | 0.608 | 0.587 | 0.536 | +0.072 | -0.052 | -0.018 | +0.023 | 0.615 | 0.243

cohort-matched (18 genes) | 2.0 | 0.632 | 0.597 | 0.525 | +0.107 | -0.042 | +0.010 | -0.004 | 0.669 | 0.204

cohort-matched (18 genes) | 3.0 | 0.668 | 0.609 | 0.514 | +0.154 | -0.075 | -0.014 | +0.029 | 0.783 | 0.138

cohort-matched (18 genes) | 4.0 | 0.668 | 0.603 | 0.507 | +0.161 | -0.095 | -0.007 | +0.041 | 0.806 | 0.136

**Table S14. Genes required for a stable estimate (five variants per gene).** [design_requirements.csv, design_simulation.csv]

Gene SD | Within-gene AUC | Estimator | Half-width | Genes | Variants

1 | 0.65 | pooled | ±0.05 | 104 | 520

1 | 0.65 | pooled | ±0.03 | 291 | 1,455

1 | 0.65 | within | ±0.05 | 172 | 860

1 | 0.65 | within | ±0.03 | 492 | 2,460

1 | 0.75 | pooled | ±0.05 | 86 | 430

1 | 0.75 | pooled | ±0.03 | 243 | 1,215

1 | 0.75 | within | ±0.05 | 133 | 665

1 | 0.75 | within | ±0.03 | 373 | 1,865

2 | 0.65 | pooled | ±0.05 | 103 | 515

2 | 0.65 | pooled | ±0.03 | 281 | 1,405

2 | 0.65 | within | ±0.05 | 220 | 1,100

2 | 0.65 | within | ±0.03 | 607 | 3,035

2 | 0.75 | pooled | ±0.05 | 100 | 500

2 | 0.75 | pooled | ±0.03 | 276 | 1,380

2 | 0.75 | within | ±0.05 | 167 | 835

2 | 0.75 | within | ±0.03 | 460 | 2,300

4 | 0.65 | pooled | ±0.05 | 110 | 550

4 | 0.65 | pooled | ±0.03 | 304 | 1,520

4 | 0.65 | within | ±0.05 | 378 | 1,890

4 | 0.65 | within | ±0.03 | 1,012 | 5,060

4 | 0.75 | pooled | ±0.05 | 113 | 565

4 | 0.75 | pooled | ±0.03 | 314 | 1,570

4 | 0.75 | within | ±0.05 | 270 | 1,350

4 | 0.75 | within | ±0.03 | 730 | 3,650

**Table S15. Calibration of out-of-fold predictions under each design.** Intercept and slope of the logistic recalibration of the labels on the predicted log-odds; a well-calibrated set of predictions has intercept 0 and slope 1, and a negative slope indicates a reversed ordering. [calibration_by_design.csv]

Model | Resampling | AUC | AUC_{W} | AUC_{B} | Cal. intercept | Cal. slope | Brier | Mean prediction

Position-level scores | variant | 0.278 | 0.446 | 0.246 | -0.48 | -10.47 | 0.261 | 0.498

Position-level scores | gene | 0.152 | 0.289 | 0.126 | -1.27 | -5.53 | 0.298 | 0.461

Full 19-feature model | variant | 0.667 | 0.601 | 0.680 | -0.35 | +0.50 | 0.231 | 0.473

Full 19-feature model | gene | 0.594 | 0.550 | 0.602 | -0.26 | +0.28 | 0.254 | 0.411

AlphaMissense only | variant | 0.884 | 0.822 | 0.896 | -0.37 | +1.05 | 0.116 | 0.453

AlphaMissense only | gene | 0.874 | 0.896 | 0.870 | -0.62 | +0.91 | 0.131 | 0.487

Across the 18 gene-held-out folds the mean prediction of the position-level model correlates with the training prevalence at r = 0.74 and with the held-out gene's own prevalence at r = −0.45; training and held-out prevalence correlate at r = −0.63 [fold_calibration.csv].

**Table S16. Fold-wise corrections and leave-two-genes-out estimation.** Each correction is applied separately within every fold and is strictly increasing, so none of them can change AUC_{W}. Prior shift subtracts the fold's training prior from the predicted log-odds; fold rank replaces each prediction by its percentile among the training predictions of the same fold. [fold_corrections.csv, leave_two_genes_out.csv]

Model | Correction | AUC | AUC_{W} | AUC_{B} | Cluster 95% CI

Position-level scores | none | 0.152 | 0.289 | 0.126 | 0.028–0.267

Position-level scores | prior shift | 0.297 | 0.289 | 0.299 | 0.118–0.570

Position-level scores | fold rank | 0.150 | 0.289 | 0.123 | 0.032–0.272

Full 19-feature model | none | 0.594 | 0.550 | 0.602 | 0.501–0.730

Full 19-feature model | prior shift | 0.630 | 0.550 | 0.645 | 0.531–0.788

Full 19-feature model | fold rank | 0.617 | 0.550 | 0.630 | 0.532–0.758

AlphaMissense only | none | 0.874 | 0.896 | 0.870 | 0.776–0.978

AlphaMissense only | prior shift | 0.888 | 0.896 | 0.886 | 0.812–0.989

AlphaMissense only | fold rank | 0.877 | 0.896 | 0.873 | 0.799–0.991

Training prevalence only | none | 0.189 | 0.500 | 0.130 | —

Training prevalence only | prior shift | 0.500 | 0.500 | 0.500 | —

Leave-two-genes-out between-gene AUC, in which every between-gene pair is ranked by one model fitted without either of its genes (153 gene pairs, 1,574 variant pairs): Position-level scores, 0.216 against 0.126 pooled leave-one-gene-out; Full 19-feature model, 0.554 against 0.602 pooled leave-one-gene-out; AlphaMissense only, 0.924 against 0.870 pooled leave-one-gene-out.

**Table S17. Influence of each gene.** The cohort is rebuilt without the named gene and every estimate recomputed. [gene_drop_influence.csv]

Gene removed | Variants | P/LP | Within-gene pairs left | Position-level AUC | AUC_{W} | AUC_{B} | Within-gene AUC: AM | DCS | Conservation

PGK1 | 34 | 16 | 10 | 0.216 | 0.600 | 0.210 | 0.800 | 0.200 | 0.000

GHRHR | 5 | 2 | 292 | 0.223 | 0.284 | 0.210 | 0.894 | 0.695 | 0.681

KAT8 | 5 | 4 | 294 | 0.249 | 0.459 | 0.203 | 0.901 | 0.684 | 0.677

BLMH | 3 | 0 | 298 | 0.141 | 0.279 | 0.113 | 0.896 | 0.681 | 0.668

GBP1 | 13 | 0 | 298 | 0.184 | 0.295 | 0.154 | 0.896 | 0.681 | 0.668

HMGCS2 | 9 | 9 | 298 | 0.548 | 0.513 | 0.558 | 0.896 | 0.681 | 0.668

LGALS3BP | 3 | 0 | 298 | 0.181 | 0.289 | 0.160 | 0.896 | 0.681 | 0.668

ANGPTL6 | 2 | 0 | 298 | 0.147 | 0.285 | 0.119 | 0.896 | 0.681 | 0.668

MSN | 1 | 1 | 298 | 0.141 | 0.285 | 0.112 | 0.896 | 0.681 | 0.668

NARS1 | 2 | 2 | 298 | 0.161 | 0.285 | 0.136 | 0.896 | 0.681 | 0.668

NT5C2 | 1 | 0 | 298 | 0.168 | 0.279 | 0.146 | 0.896 | 0.681 | 0.668

PARP3 | 4 | 0 | 298 | 0.209 | 0.292 | 0.192 | 0.896 | 0.681 | 0.668

PPARA | 1 | 0 | 298 | 0.143 | 0.285 | 0.115 | 0.896 | 0.681 | 0.668

SPRTN | 1 | 1 | 298 | 0.147 | 0.279 | 0.122 | 0.896 | 0.681 | 0.668

STAMBPL1 | 1 | 0 | 298 | 0.141 | 0.292 | 0.112 | 0.896 | 0.681 | 0.668

TBRG1 | 1 | 1 | 298 | 0.143 | 0.292 | 0.114 | 0.896 | 0.681 | 0.668

TPH1 | 1 | 0 | 298 | 0.141 | 0.275 | 0.114 | 0.896 | 0.681 | 0.668

TTC5 | 1 | 0 | 298 | 0.146 | 0.282 | 0.120 | 0.896 | 0.681 | 0.668

**Table S18. What the cohort can detect within genes.** Left, power of the two-sided stratified permutation test at level 0.05 against a score of known within-gene AUC in this cohort's gene and class structure (4,000 replicates per point). Right, observed differences in within-gene AUC between scores, with stratified bootstrap 95% intervals. [within_gene_power.csv, within_gene_paired.csv]

Target AUC_{W} | Power |  | Comparison | AUC_{W} 1 | AUC_{W} 2 | Difference | 95% interval

0.50 | 0.054 |  | DCS vs Conservation | 0.681 | 0.668 | +0.013 | -0.164 to +0.201

0.55 | 0.072 |  | RCS vs Conservation | 0.668 | 0.668 | +0.000 | -0.111 to +0.114

0.60 | 0.176 |  | DCS vs RCS | 0.681 | 0.668 | +0.013 | -0.232 to +0.269

0.65 | 0.353 |  | Normalised mobility vs Conservation | 0.483 | 0.668 | -0.185 | -0.480 to +0.141

0.70 | 0.553 |  | AlphaMissense (variant) vs DCS | 0.896 | 0.681 | +0.215 | +0.007 to +0.419

0.75 | 0.775 |  | AlphaMissense (variant) vs Conservation | 0.896 | 0.668 | +0.228 | +0.081 to +0.406

0.80 | 0.921 |  |  |  |  |  | 

0.85 | 0.989 |  |  |  |  |  | 

0.90 | 1.000 |  |  |  |  |  | 

The permutation null is distribution-free given the class counts: two independent 20,000-replicate nulls agree to 0.007 at every quantile examined. Power reaches 80% at AUC_{W} = 0.76.

**Table S19. DOME checklist.**

Section | Item | Response

Data | Provenance | ClinVar classifications for missense variants in proteins with ATLAS trajectories; 88 labelled variants, 18 genes; identifiers in Additional file 2

Data | Dataset splits | Leave-one-variant-out and leave-one-gene-out; nested inner folds grouped by gene for tuning

Data | Redundancy between splits | 3 genes contain both classes; 84% of pathogenic–benign pairs between genes; reported throughout

Data | Availability | Additional files 2 and 3; repository listed in the main text

Optimisation | Algorithms | Logistic regression (L2), random forest, gradient boosting; exploratory ensemble of five learners

Optimisation | Data encoding | Published per-protein scores; reconstructed components; within-protein encodings; 16 substitution descriptors

Optimisation | Parameters | Nine selected features of 19; C tuned from 11 values; fixed tree hyperparameters

Optimisation | Features selected in fold | Yes in controlled configurations; leakage quantified by factorial ablation

Optimisation | Fitting and regularisation | Class-balanced weights; L2 penalty; nested tuning

Model | Interpretability | Linear coefficients; AUC decomposition into within-gene and between-gene components

Model | Output | Probabilities; out-of-fold predictions provided

Evaluation | Metrics | AUC, average precision, within/between components, balanced accuracy for gene identity

Evaluation | Comparison | Gene-prevalence baseline; AlphaMissense; permutation nulls

Evaluation | Confidence | Gene-cluster bootstrap; Obuchowski clustered tests; Monte Carlo SE for permutation p

Evaluation | Availability | All code, predictions and audit in Additional file 3

## Supplementary Figures

**Figure S1. Mobility is normalised by each protein's most mobile residue.** Normalised mobility of every ClinVar-observed residue (blue) and of labelled variants (orange) in the 18 cohort proteins, ordered by the largest value observed. Numbers above give the maximum mean RMSF of the protein, which fixes the top of its scale.

**Figure S2. Ablation for three model families with the AUC decomposition.** Pooled AUC (filled), within-gene component (open circle) and between-gene component (cross) for the eight configurations of Table S5. Labels: configuration, split (V, variant; G, gene) and scaling and selection on the full labelled set (L) or in fold (F).

**Figure S3. Further simulation results.** Left, with no individual signal (β = 0), variant-level validation of the nearest-neighbour learner (solid) rises above 0.5 and gene-level validation (dashed) falls below it as between-gene heterogeneity grows. Centre, a logistic model on a single shifted feature cannot recognise genes and shows little bias. Right, the gene-prevalence lookup without individual signal.

**Figure S4. Width of gene-cluster relative to naive bootstrap intervals** for the scores and models of Table S1 († gene-held-out models).