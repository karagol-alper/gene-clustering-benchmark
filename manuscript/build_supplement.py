#!/usr/bin/env python3
"""Additional file 1: supplementary methods, tables and figures (JSON for build_docx.js)."""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mstools import Citer, Doc          # noqa: E402
from refs import REFS                  # noqa: E402
from values import load                # noqa: E402

V = load()
RES = os.path.join(os.path.dirname(HERE), "results")
FIG = os.path.join(RES, "figures")
rd = lambda n: pd.read_csv(os.path.join(RES, n))


def sci(x, digits=1):
    """0.00059 -> '5.9 × 10^{−4}'; values >= 0.001 as plain decimals."""
    x = float(x)
    if x == 0:
        return "0"
    if abs(x) >= 0.001:
        return f"{x:.3g}"
    e = int(np.floor(np.log10(abs(x))))
    return f"{x / 10 ** e:.{digits}f} × 10^{{−{-e}}}"


nm = lambda t: str(t).replace("|delta ", "|Δ ")
d = Doc(Citer(REFS), title="Additional file 1", running_head="Additional file 1 · Karagöl & Karagöl")

d.p("Additional file 1", "title")
d.p("Supplementary methods, tables and figures for: *Gene-level clustering confounds "
    "cross-validated benchmarks of structure-informed missense variant classifiers*", "authors")
d.p("Taner Karagöl (McGovern Institute for Brain Research, MIT; karagolt@mit.edu) and Alper "
    "Karagöl (Department of Mechanical Engineering, MIT; akaragol@mit.edu)", "affil")
d.p(f"All analyses use the same cohort: {V['n']} ClinVar-classified missense variants in "
    f"{V['n_genes']} genes ({V['n_p']} pathogenic or likely pathogenic, {V['n_b']} benign or "
    f"likely benign; base rate {V['base_rate']}). Every value is produced by the scripts in "
    "Additional file 3; file names in brackets refer to its results directory.", "first")

d.h(2, "Contents")
for line in ["Supplementary Methods S1–S10",
             "Table S1. Discrimination with gene-cluster and naive bootstrap intervals",
             "Table S2. Paired comparisons: DeLong and clustered tests",
             "Table S3. Per-gene composition, multipliers and within-gene discrimination",
             "Table S4. Permutation null distributions",
             "Table S5. Validation ablation for three model families",
             "Table S6. Factorial effects for logistic regression",
             "Table S7. Per-protein rescaling parameters",
             "Table S8. Gene-identity recovery",
             "Table S9. Transfer of alternative encodings",
             "Table S10. Incremental value over AlphaMissense",
             "Table S11. Within-gene discrimination of raw scores",
             "Table S12. Audit of the exploratory ensemble",
             "Table S13. Simulation summary",
             "Table S14. Design requirements",
             "Table S15. Calibration of out-of-fold predictions",
             "Table S16. Fold-wise corrections and leave-two-genes-out",
             "Table S17. Influence of each gene",
             "Table S18. What the cohort can detect within genes",
             "Table S19. DOME checklist",
             "Figures S1–S4"]:
    d.p(line, "first")

# ══════════════════════════════════════════════════════ METHODS
d.pagebreak()
d.h(2, "Supplementary Methods")
d.h(3, "S1. Within-gene and between-gene decomposition of the AUC")
d.p("Let *P* and *B* be the sets of pathogenic and benign variants, *s* a score and "
    "ψ(*a*, *b*) = 1 if *s*_{a} > *s*_{b}, 1/2 if *s*_{a} = *s*_{b} and 0 otherwise. The "
    "Mann–Whitney form of the AUC is the mean of ψ over all |*P*|·|*B*| pairs. Partition the "
    "pairs into *W*, in which both variants belong to the same gene, and its complement. Then", "first")
d.p("AUC = (1/|*P*||*B*|) Σ_{(a,b)} ψ(a,b) = *w*·AUC_{W} + (1 − *w*)·AUC_{B},", "eq")
d.p("with *w* = |*W*|/(|*P*||*B*|), AUC_{W} the mean of ψ over *W* and AUC_{B} the mean over the "
    "complement. The decomposition is exact for any score, including out-of-fold predictions. A score "
    "that is constant within genes has AUC_{W} = 1/2, so its AUC is 1/2 + (1 − *w*)(AUC_{B} − 1/2), "
    "and *w* bounds how far within-gene information can move the pooled AUC. AUC_{W} equals the "
    "pair-weighted mean of the per-gene AUCs, weighted by *n*_{P,g}·*n*_{B,g}, and this equality is "
    "exact. It is therefore a "
    "within-covariate (stratified) concordance in the spirit of the covariate-adjusted ROC of "
    "Janes and Pepe (Am J Epidemiol 2008) with gene as the covariate, differing in that genes "
    "enter in proportion to their pairs rather than their cases. Under leave-one-gene-out "
    "validation, both members of a within-gene pair are predicted by the model fitted without "
    "their gene, so AUC_{W} is a fold-internal statistic; every between-gene pair compares "
    "predictions of two different fitted models.")
d.p("Within-gene inference holds every gene's class counts fixed. The null of AUC_{W} is "
    "obtained by permuting labels within each mixed gene (20,000 Monte Carlo replicates, which "
    "approximate the exact stratified permutation null); intervals resample "
    "variants within gene and class. Because only genes with both classes contribute, AUC_{W} is "
    "undefined when no such gene is present.")

d.h(3, "S2. Reconstruction of the position-level scores")
d.p("For residue *i* of protein *k*, with conservation *c*_{i} and mean backbone RMSF *r*_{i} "
    "(three ATLAS replicas), the published quantities are", "first")
d.p("*m*_{i} = (*r*_{i} − min_{k} *r*)/(max_{k} *r* − min_{k} *r*),   "
    "DCS_{i} = 10·(*c*_{i}*m*_{i} − min_{k} *cm*)/(max_{k} *cm* − min_{k} *cm*),", "eq")
d.p("RCS_{i} = 10·(*c*_{i}(1 − *m*_{i}) − min_{k})/(max_{k} − min_{k}),   "
    "nDCS_{i} = log(1 + DCS_{i})/log 11,", "eq")
d.p("where minima and maxima run over all residues of protein *k*. Each published score is thus "
    "an affine function of a raw product with protein-specific slope and intercept. For residues "
    "with published ingredients (variants of uncertain significance at other positions of the "
    "same proteins) the identities Raw DCS = *c*·*m* and nDCS = log(1 + DCS)/log 11 hold to "
    f"machine precision. Linear regression of each published score on its raw product recovers the "
    f"slope and intercept for each of {V['n_maps']} proteins with a maximum residual of "
    f"{sci(V['map_resid'])} score units. Inverting both maps gives *c*·*m* and *c*(1 − *m*) for every "
    "scored residue; their sum is *c* and their ratio gives *m*, and the per-protein linear map "
    "between *m* and *r* gives mean RMSF. On the "
    f"{V['rec_n']} labelled residues whose ingredients are available directly, the maximum "
    f"absolute error was {sci(V['rec_err'])}. The slope of the DCS map (score units per unit of "
    "*c*·*m*) is reported as the protein's DCS multiplier (Table S7).")

d.h(3, "S3. Clustered paired comparison of AUCs")
d.p("Obuchowski's estimator treats genes *j* = 1, …, *I* as independent units. For each score, "
    "placement values *V*_{10}(*a*) = mean_{b∈B} ψ(*a*, *b*) and *V*_{01}(*b*) = "
    "mean_{a∈P} ψ(*a*, *b*) are summed within genes to *V*_{10,j} and *V*_{01,j}. With *m*_{j} and "
    "*n*_{j} the numbers of positives and negatives in gene *j*, *M* = Σ*m*_{j}, *N* = Σ*n*_{j} "
    "and θ the AUC,", "first")
d.p("S_{10} = I/((I − 1)M) Σ_{j}(V_{10,j} − m_{j}θ)^{2},  S_{01} = I/((I − 1)N) "
    "Σ_{j}(V_{01,j} − n_{j}θ)^{2},  S_{11} = I/(I − 1) Σ_{j}(V_{10,j} − m_{j}θ)(V_{01,j} − n_{j}θ),", "eq")
d.p("and Var(θ) = S_{10}/M + S_{01}/N + 2S_{11}/(MN). Covariances between two scores use the "
    "corresponding cross-products. Because the cohort has only *I* = 18 genes, the standardised "
    "difference in AUC is referred to a *t* distribution on *I* − 1 degrees of freedom rather "
    "than to the standard normal; with few clusters the variance is itself estimated imprecisely "
    "and the normal reference is anticonservative. When every gene contains one variant the "
    "estimator reduces to the DeLong "
    "variance up to the factor I/(I − 1). Intervals for differences also come from a paired "
    "gene-cluster bootstrap that resamples genes with replacement (5,000 replicates), holding "
    "out-of-fold predictions fixed.")

d.h(3, "S4. Random-intercept model")
d.p("logit P(*y*_{ij} = 1) = α + *u*_{j}, *u*_{j} ~ N(0, σ_{g}^{2}), was fitted by maximum "
    "likelihood. The random effect of each gene was integrated out by the trapezoidal rule on a "
    "uniform grid of 96,001 points spanning ±12 standard deviations; a fixed Gauss–Hermite rule, "
    "adaptive or not, is inaccurate here, because when σ_{g} is large the integrand of a "
    "class-pure gene is a normal density truncated sharply where that gene's log-odds change "
    "sign, and the resulting quadrature error inflates the maximised likelihood. The grid result "
    "was confirmed against a grid twice as fine and twice as wide (log-likelihoods agreeing to "
    "eight decimals). A profile-likelihood 95% interval for "
    "σ_{g} was obtained by maximising over α at fixed σ_{g} and locating where twice the drop in "
    f"log-likelihood reaches 3.84. The maximum-likelihood estimate is {V['sigma_g']} (95% "
    f"interval {V['sigma_ci']}); the latent-scale intraclass correlation "
    f"σ_{{g}}^{{2}}/(σ_{{g}}^{{2}} + π^{{2}}/3) is {V['icc']} ({V['icc_ci']}). The upper limit is "
    "wide because fifteen of the eighteen genes are class-pure.", "first")

d.h(3, "S5. Simulation and design calculations")
d.p("For each cell of the simulation grid, the intercept *a* was set by bisection so that the "
    "marginal pathogenic rate equalled 0.41. Each replicate drew a cohort, computed the "
    "cross-validated estimators (leave-one-variant-out and leave-one-gene-out with pooled "
    "predictions) for the three predictors, refitted each predictor on the complete simulated "
    "cohort and applied it to (i) 40 new variants in each existing gene, with that gene's *u*_{j} "
    "and *f*_{j}, for the seen-gene AUC and (ii) 400 new genes with 10 variants each for the "
    "new-gene AUC and its within-gene component. Logistic models were fitted by Newton's method "
    "with the scikit-learn objective (L2 penalty, C = 1, balanced class weights) after "
    "standardisation. The nearest-neighbour learner used the seven nearest training variants "
    "in standardised (*x̃*, *h̃*) space with inverse-distance weights; leave-one-out excluded the "
    "variant itself, and leave-one-gene-out excluded all variants of its gene. Cells were "
    "replicated 150 times.", "first")
d.p("For the design calculations, the within-gene effect β was found by bisection so that the "
    "population within-gene AUC of the score *x* equalled 0.65 or 0.75. For each number of "
    "genes, 300 cohorts were drawn and the sampling standard deviations of the pooled and "
    "within-gene AUCs of *x* were computed. Half-widths 1.96 × SD were regressed on the number "
    "of genes on log–log axes (slopes close to −1/2) and solved for ±0.05 and ±0.03.")

d.h(3, "S6. Exploratory ensemble")
d.p("Input features (twelve): DCS, RCS and mean RMSF (fixed), plus RMSF of each of the three "
    "replicas, residue-level N_{eq} of the second replica, |Δhydrophobicity|, |Δcharge|, "
    "|Δvolume|, alternate-residue hydrophobicity and an alternate-residue glycine/proline "
    "indicator, selected by ANOVA F on the complete labelled set after imputation and "
    "Yeo–Johnson transformation fitted on the complete labelled set. Learners and tuning: "
    "XGBoost, LightGBM, CatBoost, random forest and logistic regression, each tuned by 75 "
    "Optuna trials against repeated stratified 5 × 5-fold AUC on the complete labelled set. "
    "Evaluation: leave-one-variant-out predictions for each tuned learner. Blending: weights "
    "chosen by 300 Optuna trials maximising the leave-one-variant-out AUC of the weighted "
    "average. Four stages therefore had access to the held-out variant: preprocessing, "
    "selection, tuning and blending; the variant-level split additionally exposed the held-out "
    "variant's gene. The nested re-estimation of the blend is described in the main text.", "first")

d.h(3, "S7. Calibration and cross-fold estimation")
d.p("Let p_{ij} be the out-of-fold prediction for variant *i* of gene *j*.  The calibration "
    "intercept and slope are the coefficients of logit P(y = 1) = a + b·logit(p), fitted by "
    "unpenalised logistic regression; b = 1 with a = 0 denotes perfect calibration, b < 1 an "
    "overconfident predictor and b < 0 a reversed one.  Under leave-one-gene-out validation the "
    "model that predicts gene *j* is trained on a set whose pathogenic fraction π_{−j} excludes "
    "gene *j*, so π_{−j} moves away from gene *j*'s own fraction whenever that gene is "
    "class-pure.  The prior-shift correction removes this by setting", "first")
d.p("logit(p′_{ij}) = logit(p_{ij}) − logit(π_{−j}) + logit(π),", "eq")
d.p("with π the cohort fraction.  The fold-rank correction instead maps each prediction to its "
    "percentile among the predictions the same model makes on its own training set.  Both act "
    "within a fold and are strictly increasing, so AUC_{W}, which compares only variants of the "
    "same gene, is numerically unchanged; this is verified in Table S16.  A predictor that "
    "scores every variant by π_{−j} is constant within a fold, so its AUC is produced entirely "
    "by differences between folds, and the prior shift must return it to exactly 0.5.  "
    "Leave-two-genes-out fits a model on the complement of each unordered pair of genes and "
    "scores both held-out genes with it, so that each between-gene pair is ranked by one model "
    "that saw neither gene; summing over all pairs covers every between-gene pair exactly once.")

d.h(3, "S8. Power of the within-gene permutation test")
d.p("Within-gene permutation holds each gene's class counts fixed, so for a score with no ties "
    "the permutation distribution of AUC_{W} is determined by the counts alone and not by the "
    "score: the test is distribution-free given the cohort structure.  The null may therefore be "
    "computed once and reused.  To obtain power at a population within-gene AUC of τ, scores "
    "were drawn as *s* ~ N(0, 1) with √2·Φ^{−1}(τ) added for pathogenic variants, which gives "
    "P(*s*_{P} > *s*_{B}) = τ within every gene, and AUC_{W} was compared with the two-sided 5% "
    "critical values of the stored null.  Differences between two scores were summarised by a "
    "stratified bootstrap resampling variants within gene and class, which preserves the "
    "pairing that makes the two within-gene AUCs comparable.", "first")

d.h(3, "S9. Software and reproducibility")
d.p("Python 3.13; scikit-learn 1.9.1; SciPy 1.18.1; statsmodels 0.15.0; NumPy 2.5.3; pandas "
    "3.0.5; Matplotlib 3.11.2. Random seeds are fixed. `run_all.sh` regenerates every result, "
    "table and figure, and `audit.py` checks every number in the manuscript against the results "
    "directory and the citation list against the reference list. Logistic regression in "
    "permutation loops uses an exact Newton solver for the scikit-learn objective, with "
    "label-independent preprocessing cached per split; observed values computed both ways agree "
    "to three decimals (C1 0.680, C8 0.594, position-level scores 0.152). Reporting follows the "
    "DOME recommendations (Table S19).", "first")

d.h(3, "S10. Contents of Additional file 2")
cols = [
    ("gene, uniprot_accession, atlas_protein_id", "gene symbol, UniProt accession and ATLAS identifier"),
    ("clinvar_variation_id", "ClinVar variation (VCV) accession"),
    ("protein_change, residue_position, ref_aa, alt_aa", "HGVS protein change and parsed residues"),
    ("clinvar_classification, pathogenic", "ClinVar classification and binary label (1 = P/LP)"),
    ("n_variants_in_gene, gene_pathogenic_fraction", "gene size and pathogenic fraction in the cohort"),
    ("DCS_Score, RCS_Score, Normalized_DCS", "published position-level scores (per-protein 0-10)"),
    ("conservation_reconstructed, norm_dynamics_reconstructed, mean_rmsf_reconstructed_A",
     "components recovered by inverting each protein's rescaling map (Methods S2)"),
    ("dcs_protein_scale_factor", "DCS multiplier of the protein (score units per unit of conservation x mobility)"),
    ("DCS_within_protein_percentile, RCS_within_protein_percentile",
     "mid-rank percentile among all ClinVar-observed residues of the protein (label-free)"),
    ("AM_variant_score, AM_residue_mean_score", "AlphaMissense variant score and residue mean"),
    ("dHydro ... hinge_disr", "16 substitution descriptors (Methods)"),
]
d.table("**Columns of Additional file 2.**", [["Column(s)", "Content"]], [[a, b] for a, b in cols],
        widths=[3.2, 4.0], align=["l", "l"], font_size=8.5)

# ══════════════════════════════════════════════════════ TABLES
d.section(landscape=True)
t1 = rd("TableS1_discrimination.csv")
d.table("**Table S1. Discrimination with gene-cluster and naive bootstrap intervals.** 5,000 "
        "resamples each. The cluster bootstrap resamples genes; the naive bootstrap resamples "
        "variants. Ratio: width of the cluster interval over the naive interval. Obuchowski and "
        "DeLong SE: analytic standard errors treating genes or variants as independent. "
        "[TableS1_discrimination.csv]",
        [["Score or model", "AUC", "AP", "AUC_{W}", "AUC_{B}", "Cluster 95% CI", "Naive 95% CI",
          "Ratio", "SE clustered", "SE DeLong"]],
        [[nm(r.score), f"{r.auc:.3f}", f"{r.ap:.3f}", f"{r.auc_within:.3f}", f"{r.auc_between:.3f}",
          r.cluster_ci.replace("-", "–"), r.naive_ci.replace("-", "–"), f"{r.width_ratio:.2f}",
          f"{r.obuchowski_se:.3f}", f"{r.delong_se:.3f}"] for _, r in t1.iterrows()],
        widths=[3.4, 0.7, 0.7, 0.7, 0.7, 1.2, 1.2, 0.6, 0.8, 0.8],
        align=["l"] + ["c"] * 9, font_size=8,
        notes=f"† Gene-held-out logistic model (out-of-fold predictions). Ratios exceed 1 for "
              f"{int((t1.width_ratio > 1).sum())} of {len(t1)} rows; the largest belong to scores "
              f"whose pooled AUC rests almost entirely on between-gene pairs, where the "
              f"clustering is, so that treating variants as independent understates the "
              f"uncertainty most.")

t2 = rd("TableS2_paired_tests.csv")
d.table("**Table S2. Paired comparisons of AUCs.** DeLong treats variants as independent; the "
        "clustered test (Obuchowski) treats genes as independent. [TableS2_paired_tests.csv]",
        [["Comparison", "AUC 1", "AUC 2", "Difference", "DeLong z", "DeLong p", "Clustered z",
          "Clustered p"]],
        [[nm(r.comparison), f"{r.auc_1:.3f}", f"{r.auc_2:.3f}", f"{r.difference:+.3f}", f"{r.delong_z:.2f}",
          sci(r.delong_p), f"{r.clustered_z:.2f}", sci(r.clustered_p)] for _, r in t2.iterrows()],
        widths=[5.0, 0.7, 0.7, 0.8, 0.7, 0.9, 0.8, 0.9], align=["l"] + ["c"] * 7, font_size=8,
        notes="† Gene-held-out logistic model.")

t3 = rd("TableS_per_gene.csv")
fmt = lambda v: "—" if pd.isna(v) else f"{v:.2f}"
d.table("**Table S3. Per-gene composition, DCS multiplier and within-gene discrimination.** "
        "Within-gene AUCs are defined only for genes containing both classes. [TableS_per_gene.csv]",
        [["Gene", "UniProt", "n", "P/LP", "B/LB", "Pathogenic fraction", "Within pairs",
          "DCS multiplier", "AUC AM", "AUC DCS", "AUC RCS", "AUC conservation"]],
        [[r.Gene_Name, r.uniprot, int(r.n), int(r.n_pathogenic), int(r.n_benign),
          f"{r.pathogenic_fraction:.2f}", int(r.within_pairs), f"{r.dcs_multiplier:.1f}",
          fmt(r.auc_AM), fmt(r.auc_DCS), fmt(r.auc_RCS), fmt(r.auc_conservation)]
         for _, r in t3.iterrows()],
        widths=[1.1, 0.9, 0.5, 0.5, 0.5, 0.9, 0.7, 0.8, 0.7, 0.7, 0.7, 0.9], font_size=8,
        notes=f"Gene–class association: χ² = {V['chi2']}, df = {V['chi2_df']}, Monte Carlo "
              f"p {V['chi2_p_mc']} ({V['frac_exp_lt5']}% of expected counts below five; "
              f"asymptotic p not valid); Cramér's V = {V['cramer']}.")

t4 = rd("permutation_summary.csv")
names = {"C1": "Variant-level, leaked (C1)", "C8": "Gene-level, controlled (C8)",
         "POS": "Position-level scores, gene-held-out",
         "POSV": "Position-level scores, variant-level"}
nulls = {"global": "global", "within_gene": "within gene", "donor": "donor protein"}
d.table("**Table S4. Permutation null distributions.** Complete pipeline refitted in every "
        "replicate. p (upper): proportion of null AUCs at or above the observed value; p (lower): "
        "at or below; p (2-sided): twice the smaller tail; all (b + 1)/(B + 1). "
        "[permutation_summary.csv, permutation_distributions.csv]",
        [["Pipeline", "Null", "Replicates", "Observed", "Null mean", "SD", "2.5%", "97.5%",
          "p (upper)", "p (lower)", "p (2-sided)", "Null AUC_{W}", "Null AUC_{B}"]],
        [[names[r.design], nulls[r.null], f"{int(r.n):,}", f"{r.observed:.3f}", f"{r['mean']:.3f}",
          f"{r.sd:.3f}", f"{r.pct2_5:.3f}", f"{r.pct97_5:.3f}", f"{r.p_upper:.3f}",
          f"{r.p_lower:.3f}", f"{min(1.0, 2 * min(r.p_upper, r.p_lower)):.3f}",
          f"{r.mean_within:.3f}", f"{r.mean_between:.3f}"]
         for _, r in t4.iterrows()],
        widths=[2.2, 1.0, 0.7, 0.65, 0.7, 0.55, 0.55, 0.55, 0.7, 0.7, 0.75, 0.75, 0.75],
        align=["l", "l"] + ["c"] * 11, font_size=8)

t5 = rd("ablation_full.csv")
d.table("**Table S5. Validation ablation for three model families (2^{3} design, fixed "
        "hyperparameters).** Logistic regression C = 1; random forest 400 trees, minimum leaf 3; "
        "gradient boosting 60 trees of depth 2. Each cell: AUC (AP) and, below, AUC_{W}/AUC_{B}. "
        "[ablation_full.csv]",
        [["Config.", "Split", "Scaling", "Selection", "Logistic regression", "Random forest",
          "Gradient boosting"]],
        [[r.config, r.split, r.scaling, r.selection] +
         [f"{r[m]:.3f} ({r[m + '_ap']:.3f}); {r[m + '_within']:.3f}/{r[m + '_between']:.3f}"
          for m in ("lr", "rf", "gb")] for _, r in t5.iterrows()],
        widths=[0.7, 0.8, 0.8, 0.8, 2.2, 2.2, 2.2], font_size=8)

t6 = rd("ablation_effects.csv")
d.table("**Table S6. Factorial effects on AUC for logistic regression (2^{4} design of Table "
        "2 in the main text).** Effects are standard contrasts: mean AUC with the factor at 'sees held-out data' "
        "(variant split; full-set scaling, selection or tuning) minus mean AUC at the other level; "
        "interactions are the corresponding products. [ablation_effects.csv]",
        [["Effect", "Order", "AUC effect"]],
        [[r.effect, int(r.order), f"{r.auc_effect:+.4f}"] for _, r in
         t6.sort_values(["order", "auc_effect"], key=lambda s: s if s.name == "order" else -s.abs()).iterrows()],
        widths=[3.0, 0.8, 1.0], align=["l", "c", "c"], font_size=8)

t7 = rd("TableS_provenance.csv")
t7 = t7.sort_values("in_cohort", ascending=False)
d.table("**Table S7. Per-protein rescaling parameters recovered from the published scores.** "
        "DCS and RCS multipliers: score units per unit of raw product. RMSF at normalised "
        "mobility 0 and 1: the protein's minimum and maximum mean RMSF. Max m (ClinVar): the "
        "largest normalised mobility among ClinVar-observed residues. [TableS_provenance.csv]",
        [["Protein", "Gene", "Cohort", "Reference residues", "DCS multiplier", "RCS multiplier",
          "RMSF min (Å)", "RMSF max (Å)", "Max m (ClinVar)", "Max residual"]],
        [[r.Protein_ID, r.Gene_Name, "yes" if r.in_cohort else "no", int(r.n_ref),
          f"{r.dcs_slope:.2f}", f"{r.rcs_slope:.2f}", f"{r.rmsf_at_nd0:.2f}", f"{r.rmsf_at_nd1:.2f}",
          f"{r.max_nd_clinvar_residues:.2f}", f"{r.max_abs_residual:.1e}"] for _, r in t7.iterrows()],
        widths=[1.1, 1.0, 0.6, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9], font_size=7.5)

t8 = rd("gene_identity_recovery.csv")
setn = {"labelled": "Labelled cohort", "all_clinvar": "All ClinVar residues"}
d.table("**Table S8. Gene-identity recovery from each encoding.** Random forest, five-fold "
        "cross-validation stratified by gene and grouped by residue, five repeats; null from 50 "
        "label permutations. [gene_identity_recovery.csv]",
        [["Set", "Encoding", "Genes", "Balanced accuracy", "Accuracy", "Majority rate",
          "Null mean", "Null SD", "p"]],
        [[setn[r.set], r.description, int(r.k), f"{r.bacc:.3f}", f"{r.acc:.3f}", f"{r.majority:.3f}",
          f"{r.null_mean:.3f}", f"{r.null_sd:.3f}", f"{r.p:.3f}"] for _, r in t8.iterrows()],
        widths=[1.4, 4.0, 0.6, 1.0, 0.8, 0.9, 0.8, 0.7, 0.6], align=["l", "l"] + ["c"] * 7, font_size=7.5)

t9 = rd("transfer_by_encoding.csv")
d.table("**Table S9. Gene-held-out and variant-level AUC of alternative encodings.** Logistic "
        "regression; preprocessing in fold. [transfer_by_encoding.csv]",
        [["Encoding", "Features", "Gene-held-out AUC", "AP", "AUC_{W}", "AUC_{B}",
          "Variant-level AUC"]],
        [[r.description, int(r.n_features), f"{r.auc_gene:.3f}", f"{r.ap_gene:.3f}",
          f"{r.within_gene:.3f}", f"{r.between_gene:.3f}", f"{r.auc_variant:.3f}"]
         for _, r in t9.iterrows()],
        widths=[4.6, 0.8, 1.1, 0.7, 0.8, 0.8, 1.1], align=["l"] + ["c"] * 6, font_size=8)

t10 = rd("incremental_over_am.csv")
d.table("**Table S10. Incremental value over AlphaMissense (gene-held-out logistic regression).** "
        "Δ: AUC minus that of AlphaMissense alone. CI: paired gene-cluster bootstrap. "
        "[incremental_over_am.csv]",
        [["Model", "Selection", "AUC", "AP", "AUC_{W}", "AUC_{B}", "Δ", "Δ 95% CI", "DeLong p",
          "Clustered p"]],
        [[r.model, r.feature_selection, f"{r.auc:.3f}", f"{r.ap:.3f}", f"{r['within']:.3f}",
          f"{r['between']:.3f}", "—" if pd.isna(r.get("delta")) else f"{r.delta:+.3f}",
          "—" if pd.isna(r.get("delta_cluster_ci_lo")) else
          f"{r.delta_cluster_ci_lo:+.3f} to {r.delta_cluster_ci_hi:+.3f}",
          "—" if pd.isna(r.get("delong_p")) else f"{r.delong_p:.3f}",
          "—" if pd.isna(r.get("clustered_p")) else f"{r.clustered_p:.3f}"] for _, r in t10.iterrows()],
        widths=[3.8, 0.7, 0.6, 0.6, 0.6, 0.6, 0.7, 1.4, 0.7, 0.8], align=["l"] + ["c"] * 9, font_size=7.5)

t11 = rd("within_gene_discrimination.csv")
fm = lambda v: "—" if pd.isna(v) else f"{v:.2f}"
d.table("**Table S11. Within-gene discrimination of raw scores.** AUC_{W}: concordance over the "
        f"{V['pairs_w']} within-gene pairs; CI from a stratified bootstrap; p from stratified "
        "permutation, two-sided on |AUC_{W} − 0.5| (20,000 Monte Carlo replicates). CLR: conditional logistic regression stratified by "
        "gene, odds ratio per cohort SD. [within_gene_discrimination.csv]",
        [["Score", "AUC", "AUC_{B}", "AUC_{W}", "AUC_{W} 95% CI", "p", "PGK1", "GHRHR", "KAT8",
          "CLR OR (95% CI)", "CLR p"]],
        [[nm(r.score), f"{r.auc:.3f}", f"{r.auc_between:.3f}", f"{r.auc_within:.3f}",
          r.within_ci.replace("-", "–"), "< 0.0001" if r.within_p_two_sided == 0 else f"{r.within_p_two_sided:.4f}", fm(r.get("auc_PGK1")),
          fm(r.get("auc_GHRHR")), fm(r.get("auc_KAT8")),
          f"{r.clr_or_per_sd:.2f} ({r.clr_or_ci.replace('-', '–')})", f"{r.clr_p:.3f}"]
         for _, r in t11.iterrows()],
        widths=[2.6, 0.6, 0.6, 0.6, 1.1, 0.6, 0.6, 0.6, 0.6, 1.4, 0.6], align=["l"] + ["c"] * 10,
        font_size=7.5)

e = json.load(open(os.path.join(RES, "ensemble_audit.json")))
rows = [["Reported (in-sample blend)", f"{e['auc_reported']:.3f}"],
        ["Best of 3,006 candidate weightings, in sample", f"{e['auc_in_sample_best']:.3f}"],
        ["Equal weights", f"{e['auc_equal_weights']:.3f}"],
        ["Blend with weights nested in the leave-one-out loop", f"{e['auc_nested_blend']:.3f}"],
        ["Single best learner, nested choice", f"{e['auc_nested_best_single']:.3f}"]] + \
       [[f"Base learner: {k}", f"{v:.3f}"] for k, v in e["auc_base_models"].items()] + \
       [["Recovered weights", ", ".join(f"{k} {100 * v:.1f}%" for k, v in e["recovered_weights"].items())],
        ["Within-gene / between-gene AUC (reported)",
         f"{e['ensemble_decomposition']['auc_within']:.3f} / {e['ensemble_decomposition']['auc_between']:.3f}"],
        ["Gene-cluster 95% CI (reported)", f"{e['ensemble_cluster_ci'][0]:.3f}–{e['ensemble_cluster_ci'][1]:.3f}"],
        ["Versus gene-prevalence lookup: difference, DeLong p, clustered p",
         f"{e['ens_vs_prev']['delta']:+.3f}, {e['ens_vs_prev']['delong_p']:.2f}, {e['ens_vs_prev']['clustered_p']:.2f}"],
        ["Spearman correlation with gene prevalence (in sample)", f"{e['spearman_ens_prevalence_in_sample']:.2f}"]]
d.table("**Table S12. Audit of the exploratory ensemble.** [ensemble_audit.json]",
        [["Quantity", "Value"]], rows, widths=[5.0, 4.0], align=["l", "l"], font_size=8.5)

sim = rd("simulation.csv")
sel = sim[(sim.beta == 0.75) & (sim.gamma == 1.0)]
d.table("**Table S13. Simulation summary (within-gene effect β = 0.75, feature offset γ = 1).** "
        "Means over 150 replicates. kNN: seven-nearest-neighbour learner on two descriptors. "
        "Bias columns: estimator minus estimand. [simulation.csv; full grid including β = 0 and "
        "γ = 0 in the file]",
        [["Design", "Gene SD", "kNN variant-LOO", "kNN seen-gene", "kNN new-gene",
          "Variant-LOO − new", "Gene-LOO − new", "Gene-LOO (W) − within", "LR variant-LOO − new",
          "Lookup variant-LOO", "Lookup gene-LOO"]],
        [[r.design.replace(" x ", " × "), f"{r.gene_sd:.1f}", f"{r.knn_vloo_mean:.3f}",
          f"{r.knn_true_seen_mean:.3f}", f"{r.knn_true_new_mean:.3f}", f"{r.knn_infl_vloo_new_mean:+.3f}",
          f"{r.knn_bias_gloo_new_mean:+.3f}", f"{r.knn_bias_gloo_w_mean:+.3f}",
          f"{r.feat_infl_vloo_new_mean:+.3f}", f"{r.prev_vloo_mean:.3f}", f"{r.prev_gloo_mean:.3f}"]
         for _, r in sel.sort_values(["design", "gene_sd"]).iterrows()],
        widths=[1.9, 0.6, 0.8, 0.8, 0.8, 0.9, 0.9, 1.0, 0.9, 0.8, 0.8], align=["l"] + ["c"] * 10,
        font_size=7)

req = rd("design_requirements.csv")
d.table("**Table S14. Genes required for a stable estimate (five variants per gene).** "
        "[design_requirements.csv, design_simulation.csv]",
        [["Gene SD", "Within-gene AUC", "Estimator", "Half-width", "Genes", "Variants"]],
        [[f"{r.gene_sd:.0f}", f"{r.target_within_auc:.2f}", r.estimator, f"±{r.half_width:.2f}",
          f"{r.genes_required:,}", f"{r.variants_required:,}"] for _, r in req.iterrows()],
        widths=[0.8, 1.0, 1.0, 0.9, 0.9, 0.9], font_size=8)

cb = rd("calibration_by_design.csv")
d.table("**Table S15. Calibration of out-of-fold predictions under each design.** Intercept and "
        "slope of the logistic recalibration of the labels on the predicted log-odds; a "
        "well-calibrated set of predictions has intercept 0 and slope 1, and a negative slope "
        "indicates a reversed ordering. [calibration_by_design.csv]",
        [["Model", "Resampling", "AUC", "AUC_{W}", "AUC_{B}", "Cal. intercept", "Cal. slope",
          "Brier", "Mean prediction"]],
        [[r.model, "variant" if r.split == "variant" else "gene", f"{r.auc:.3f}",
          f"{r.auc_within:.3f}", f"{r.auc_between:.3f}", f"{r.cal_intercept:+.2f}",
          f"{r.cal_slope:+.2f}", f"{r.brier:.3f}", f"{r.mean_pred:.3f}"]
         for _, r in cb.iterrows()],
        widths=[2.0, 1.0, 0.8, 0.8, 0.8, 1.1, 1.0, 0.8, 1.2], align=["l", "l"] + ["c"] * 7,
        font_size=8,
        notes=f"Across the {V['n_genes']} gene-held-out folds the mean prediction of the "
              f"position-level model correlates with the training prevalence at "
              f"r = {V['fold_r_train']} and with the held-out gene's own prevalence at "
              f"r = {V['fold_r_gene']}; training and held-out prevalence correlate at "
              f"r = {V['fold_r_train_gene']} [fold_calibration.csv].")

fx = rd("fold_corrections.csv")
lpg = rd("leave_two_genes_out.csv").set_index("model")
d.table("**Table S16. Fold-wise corrections and leave-two-genes-out estimation.** Each correction "
        "is applied separately within every fold and is strictly increasing, so none of them can "
        "change AUC_{W}. Prior shift subtracts the fold's training prior from the predicted "
        "log-odds; fold rank replaces each prediction by its percentile among the training "
        "predictions of the same fold. [fold_corrections.csv, leave_two_genes_out.csv]",
        [["Model", "Correction", "AUC", "AUC_{W}", "AUC_{B}", "Cluster 95% CI"]],
        [[r.model, r.correction, f"{r.auc:.3f}", f"{r.auc_within:.3f}", f"{r.auc_between:.3f}",
          r.cluster_ci.replace("-", "–")] for _, r in fx.iterrows()],
        widths=[2.2, 1.2, 0.9, 0.9, 0.9, 1.4], align=["l", "l"] + ["c"] * 4, font_size=8,
        notes="Leave-two-genes-out between-gene AUC, in which every between-gene pair is ranked "
              "by one model fitted without either of its genes (" + V["n_gene_pairs"] +
              " gene pairs, " + V["lpgo_n_pairs"] + " variant pairs): " +
              "; ".join(f"{m}, {lpg.loc[m, 'auc_between_lpgo']:.3f} against "
                        f"{lpg.loc[m, 'auc_between_logo']:.3f} pooled leave-one-gene-out"
                        for m in lpg.index) + ".")

gdi = rd("gene_drop_influence.csv")
d.table("**Table S17. Influence of each gene.** The cohort is rebuilt without the named gene and "
        "every estimate recomputed. [gene_drop_influence.csv]",
        [["Gene removed", "Variants", "P/LP", "Within-gene pairs left", "Position-level AUC",
          "AUC_{W}", "AUC_{B}", "Within-gene AUC: AM", "DCS", "Conservation"]],
        [[r.dropped, int(r.n_removed), int(r.n_pathogenic_removed), int(r.n_within_pairs),
          f"{r.pos_auc:.3f}", f"{r.pos_within:.3f}", f"{r.pos_between:.3f}",
          f"{r['within_AlphaMissense (variant)']:.3f}", f"{r.within_DCS:.3f}",
          f"{r.within_Conservation:.3f}"] for _, r in gdi.iterrows()],
        widths=[1.3, 0.8, 0.6, 1.3, 1.2, 0.8, 0.8, 1.3, 0.8, 1.1],
        align=["l"] + ["c"] * 9, font_size=7.5)

pwr = rd("within_gene_power.csv")
pcm = rd("within_gene_paired.csv")
d.table("**Table S18. What the cohort can detect within genes.** Left, power of the two-sided "
        "stratified permutation test at level 0.05 against a score of known within-gene AUC in "
        "this cohort's gene and class structure (4,000 replicates per point). Right, observed "
        "differences in within-gene AUC between scores, with stratified bootstrap 95% intervals. "
        "[within_gene_power.csv, within_gene_paired.csv]",
        [["Target AUC_{W}", "Power", "", "Comparison", "AUC_{W} 1", "AUC_{W} 2", "Difference",
          "95% interval"]],
        [([f"{pwr.target_within_auc.iloc[i]:.2f}", f"{pwr.power.iloc[i]:.3f}", ""]
          if i < len(pwr) else ["", "", ""]) +
         ([pcm.comparison.iloc[i], f"{pcm.within_1.iloc[i]:.3f}", f"{pcm.within_2.iloc[i]:.3f}",
           f"{pcm.difference.iloc[i]:+.3f}",
           f"{pcm.ci_lo.iloc[i]:+.3f} to {pcm.ci_hi.iloc[i]:+.3f}"]
          if i < len(pcm) else ["", "", "", "", ""])
         for i in range(max(len(pwr), len(pcm)))],
        widths=[1.2, 0.8, 0.3, 2.6, 0.9, 0.9, 0.9, 1.5],
        align=["c", "c", "c", "l", "c", "c", "c", "c"], font_size=8,
        notes=f"The permutation null is distribution-free given the class counts: two "
              f"independent 20,000-replicate nulls agree to {V['null_df_check']} at every "
              f"quantile examined. Power reaches 80% at AUC_{{W}} = {V['mdes_within']}.")

dome = [
    ["Data", "Provenance", "ClinVar classifications for missense variants in proteins with ATLAS trajectories; "
     f"{V['n']} labelled variants, {V['n_genes']} genes; identifiers in Additional file 2"],
    ["Data", "Dataset splits", "Leave-one-variant-out and leave-one-gene-out; nested inner folds grouped by gene for tuning"],
    ["Data", "Redundancy between splits", f"{V['n_mixed']} genes contain both classes; {V['b_pct']}% of pathogenic–benign pairs between genes; reported throughout"],
    ["Data", "Availability", "Additional files 2 and 3; repository listed in the main text"],
    ["Optimisation", "Algorithms", "Logistic regression (L2), random forest, gradient boosting; exploratory ensemble of five learners"],
    ["Optimisation", "Data encoding", "Published per-protein scores; reconstructed components; within-protein encodings; 16 substitution descriptors"],
    ["Optimisation", "Parameters", "Nine selected features of 19; C tuned from 11 values; fixed tree hyperparameters"],
    ["Optimisation", "Features selected in fold", "Yes in controlled configurations; leakage quantified by factorial ablation"],
    ["Optimisation", "Fitting and regularisation", "Class-balanced weights; L2 penalty; nested tuning"],
    ["Model", "Interpretability", "Linear coefficients; AUC decomposition into within-gene and between-gene components"],
    ["Model", "Output", "Probabilities; out-of-fold predictions provided"],
    ["Evaluation", "Metrics", "AUC, average precision, within/between components, balanced accuracy for gene identity"],
    ["Evaluation", "Comparison", "Gene-prevalence baseline; AlphaMissense; permutation nulls"],
    ["Evaluation", "Confidence", "Gene-cluster bootstrap; Obuchowski clustered tests; Monte Carlo SE for permutation p"],
    ["Evaluation", "Availability", "All code, predictions and audit in Additional file 3"],
]
d.table("**Table S19. DOME checklist.**", [["Section", "Item", "Response"]], dome,
        widths=[1.2, 2.0, 7.0], align=["l", "l", "l"], font_size=8)

# ══════════════════════════════════════════════════════ FIGURES
d.section(landscape=False)
d.h(2, "Supplementary Figures")
d.figure(f"{FIG}/FigureS1.png",
         "**Figure S1. Mobility is normalised by each protein's most mobile residue.** Normalised "
         "mobility of every ClinVar-observed residue (blue) and of labelled variants (orange) in "
         "the 18 cohort proteins, ordered by the largest value observed. Numbers above give the "
         "maximum mean RMSF of the protein, which fixes the top of its scale.")
d.figure(f"{FIG}/FigureS2.png",
         "**Figure S2. Ablation for three model families with the AUC decomposition.** Pooled AUC "
         "(filled), within-gene component (open circle) and between-gene component (cross) for "
         "the eight configurations of Table S5. Labels: configuration, split (V, variant; G, gene) "
         "and scaling and selection on the full labelled set (L) or in fold (F).")
d.figure(f"{FIG}/FigureS3.png",
         "**Figure S3. Further simulation results.** Left, with no individual signal (β = 0), "
         "variant-level validation of the nearest-neighbour learner (solid) rises above 0.5 and "
         "gene-level validation (dashed) falls below it as between-gene heterogeneity grows. "
         "Centre, a logistic model on a single shifted feature cannot recognise genes and shows "
         "little bias. Right, the gene-prevalence lookup without individual signal.")
d.figure(f"{FIG}/FigureS4.png",
         "**Figure S4. Width of gene-cluster relative to naive bootstrap intervals** for the scores "
         "and models of Table S1 († gene-held-out models).", width_in=4.6)

d.save(os.path.join(HERE, "supplement.json"), os.path.join(HERE, "supplement.md"))
print("supplement written")
