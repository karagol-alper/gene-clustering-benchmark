#!/usr/bin/env python3
"""Generate the manuscript (JSON for build_docx.js, plus a Markdown copy used by
the audit) with every quantitative statement filled from results/."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mstools import Citer, Doc          # noqa: E402
from refs import REFS                  # noqa: E402
from values import load                # noqa: E402

V = load()
FIG = os.path.join(os.path.dirname(HERE), "results", "figures")
d = Doc(Citer(REFS), title="Gene-level clustering confounds cross-validated benchmarks of "
        "structure-informed missense variant classifiers",
        running_head="Gene clustering in variant classifier benchmarks", line_numbers=True,
        creator="T. Karagöl and A. Karagöl")

# ───────────────────────────────────────────────────────────── front matter
d.p("Gene-level clustering confounds cross-validated benchmarks of structure-informed "
    "missense variant classifiers", "title")
d.p("Taner Karagöl^{1,¶,*} and Alper Karagöl^{2,3,¶,*}", "authors")
d.p("^{1} McGovern Institute for Brain Research, Massachusetts Institute of Technology, "
    "Cambridge, MA 02139, USA", "affil")
d.p("^{2} Department of Mechanical Engineering, Massachusetts Institute of Technology, "
    "Cambridge, MA 02139, USA", "affil")
d.p("^{3} Broad Institute of MIT and Harvard, Cambridge, MA 02142, USA", "affil")
d.p("^{¶} These authors contributed equally. ^{*} Correspondence: karagolt@mit.edu (T.K.); "
    "akaragol@mit.edu (A.K.)", "affil")
d.p("**Running head:** Gene clustering in variant classifier benchmarks", "affil")

# ───────────────────────────────────────────────────────────── abstract
d.h(2, "Abstract")
d.p(f"**Background:** Structure- and dynamics-informed features are increasingly proposed as "
    f"complements to sequence-based missense pathogenicity predictors. Because clinically "
    f"labelled variants are distributed unevenly across genes, the design of a cross-validated "
    f"benchmark can determine what it measures. Dynamics-aware conservation scores serve as "
    f"a worked case.", "abstract")
d.p(f"**Methods:** We assembled {V['n']} ClinVar missense variants ({V['n_p']} pathogenic, "
    f"{V['n_b']} benign) in {V['n_genes']} genes with molecular dynamics trajectories, "
    f"reconstructed the components of two position-level scores, and evaluated classifiers in a "
    f"factorial design crossing the resampling unit with the placement of preprocessing, feature "
    f"selection and tuning. Every area under the ROC curve (AUC) was decomposed exactly into "
    f"within-gene and between-gene components and calibrated against gene-prevalence baselines, "
    f"structure-preserving permutation nulls, clustered tests against AlphaMissense, fold-wise "
    f"calibration diagnostics and simulations with known estimands.", "abstract")
d.p(f"**Results:** Fifteen of {V['n_genes']} genes contained one class, so {V['b_pct']}% of "
    f"pathogenic–benign pairs spanned genes. A lookup of each gene's "
    f"pathogenic fraction reached AUC {V['prev_variant_loo']} under leave-one-variant-out "
    f"validation, entirely between genes, and a tuned five-model ensemble reached {V['ens']}, "
    f"a difference not distinguishable from zero "
    f"(clustered p = {V['ens_vs_prev_p_clu']}). Against labels "
    f"permuted within genes, the leaked pipeline "
    f"was not significant (p = {V['perm_C1_within_gene_pu']}). The position-level scores "
    f"identified their protein almost as well before as after per-protein rescaling (balanced "
    f"accuracy {V['gid_labelled_unrescaled']} and {V['gid_labelled_published']}; chance 0.125). "
    f"Their within-gene AUC ({V['wg_dcs']}, p = {V['wg_dcs_p']}) matched that of their "
    f"conservation factor ({V['wg_cons']}); normalised mobility gave {V['wg_nd']}. Gene-held-out, "
    f"the scores gave {V['obs_POS']}, against {V['prev_gene_loo']} for a predictor knowing only "
    f"training prevalence; subtracting each fold's training prior raised this to "
    f"{V['fx_pos_prior']}, and ranking every between-gene pair with a model blind to both "
    f"genes raised the between-gene component from {V['lpgo_pos_logo']} to "
    f"{V['lpgo_pos']}. No "
    f"addition detectably improved gene-held-out AlphaMissense ({V['inc_am']}; all clustered "
    f"p ≥ {V['inc_min_pc']}). Simulations showed variant-level inflation and pooled gene-level "
    f"deflation growing with between-gene heterogeneity, while the within-gene component "
    f"tracked its estimand. The cohort reaches 80% power only against a within-gene AUC of "
    f"{V['mdes_within']}.", "abstract")
d.p(f"**Conclusions:** On gene-clustered cohorts the pooled AUC is predominantly a between-gene "
    f"statistic: variant-level validation rewards gene memorisation and pooled gene-level "
    f"validation penalises gene-level information, partly correctably. Within- and between-gene "
    f"components, with matched baselines, nulls and clustered inference, separate these effects. "
    f"Here the dynamics component added no within-gene discrimination beyond conservation.",
    "abstract")
d.p("**Keywords:** variant effect prediction; benchmarking; data leakage; cross-validation; "
    "clustered data; within-gene concordance; AlphaMissense; protein dynamics", "abstract")

# ───────────────────────────────────────────────────────────── background
d.h(2, "Background")
d.p("Computational prediction of missense variant effects is routine in clinical genomics, and "
    "published predictors far outnumber independent evaluations of them. The cohorts available "
    "for evaluation are small, unevenly distributed across genes and partly circular with "
    "respect to the predictors under test [@grimm2015,@livesey2023,@gunning2021]. Grimm and "
    "colleagues showed that variants from the same gene and overlap with training data both "
    "inflate apparent accuracy, and that correcting for them reorders published rankings "
    "[@grimm2015]; later assessments confirm that evaluation design accounts for much of the "
    "reported difference between tools [@livesey2023,@gunning2021].", "first")
d.p("Structural and dynamic descriptors are increasingly added to variant interpretation. "
    "Predictors draw on predicted structure [@jumper2021,@abramson2024], and features derived "
    "from molecular dynamics, conformational ensembles and coevolutionary coupling networks have "
    "been proposed [@xu2025,@sutto2015,@rodrigues2018]. Bahar and colleagues showed that "
    "elastic-network dynamics descriptors carry information about variant pathogenicity "
    "[@ponzoni2018] and built them into Rhapsody [@ponzoni2020]; its successor, Rhapsody-2, was "
    "trained on 117,525 variants in 12,094 proteins and reached an AUROC of 0.94 under "
    "cross-validation that excluded every variant of a test protein from training "
    "[@banerjee2025]. Holding out whole proteins is the design choice examined here.")
d.p("Many structure- and dynamics-derived descriptors are position-level quantities: they "
    "describe a residue rather than a substitution, take one value for every substitution at a "
    "site, and are often normalised within each protein so that residues can be ranked along a "
    "sequence. Two properties of clinically labelled cohorts interact with such features. First, "
    "genes differ greatly in the fraction of their labelled variants that are pathogenic, and "
    "many contribute variants of one class only. Second, the area under the ROC curve (AUC), the "
    "usual summary of discrimination, averages over all pathogenic–benign pairs, including pairs "
    "drawn from different genes. Together these mean that when most genes are class-pure, most "
    "pairs lie between genes and the AUC largely measures how well a score ranks genes rather "
    "than variants. The diagnostic-testing literature separates the two through the "
    "covariate-adjusted ROC curve, which compares cases and controls sharing a covariate value "
    "[@janes2008], but the distinction is rarely drawn in variant-effect benchmarking.")
d.p("Information leakage into cross-validated estimates is well documented. Fitting any part of a pipeline, whether scaling, "
    "feature selection or hyperparameter search, on data that include the evaluation fold "
    "produces optimistic estimates [@ambroise2002,@varma2006,@cawley2010,@wainer2021], and "
    "dependency in the data invalidates random splitting unless the resampling unit matches the "
    "unit of dependence [@roberts2017,@eid2021]. Kapoor and Narayanan catalogued leakage across "
    "294 papers in 17 fields [@kapoor2023], Whalen and colleagues described its genomic forms "
    "[@whalen2022], and the DOME recommendations ask authors to state the resampling unit "
    "explicitly [@walsh2021]. The converse problem has received less attention: pooled AUCs "
    "from leave-one-out designs are biased when held-out units shift the class balance of the "
    "training data [@parker2007,@airola2011,@forman2010], and leaving out a gene whose variants "
    "all belong to one class is an extreme instance of such a shift. Gene-level information is "
    "also not illegitimate in itself: the ACMG/AMP criteria PP2 and BP1 encode gene-level rates "
    "of pathogenic and benign missense variation [@richards2015], so a predictor that has "
    "learnt gene-level prevalence duplicates evidence that a clinical classification already "
    "uses.")
d.p("A concrete case motivates what follows. Having defined a dynamics-aware decomposition of "
    "evolutionary conservation into a Rigid Conserved Score (RCS) and a Dynamic Conserved Score "
    "(DCS) [@karagol2026], we tested whether the scores improve missense classification in a "
    f"supervised model. An exploratory five-model ensemble reached an AUC of {V['ens']} under "
    "leave-one-variant-out cross-validation, a value that did not survive gene-aware "
    "validation. That discrepancy is the object of study here. We decompose the AUC "
    "into within-gene and between-gene components; separate the pipeline stages that leak "
    "information in a factorial design that includes hyperparameter tuning and ensemble "
    "blending; reconstruct the components of the position-level scores to establish what they "
    "encode; calibrate gene-held-out estimates against null distributions that preserve the "
    "cohort's gene structure; and generalise the findings by simulation with known estimands. "
    "The dynamics-aware scores serve as a worked example, and the analysis applies to any "
    "feature evaluated on a gene-clustered cohort.")

# ───────────────────────────────────────────────────────────── methods
d.h(2, "Methods")
d.h(3, "Benchmark cohort")
d.p("Missense variants were retrieved from ClinVar [@landrum2014] for human proteins with "
    "replica molecular dynamics trajectories in the ATLAS database [@vandermeersche2024], as "
    "described previously [@karagol2026]. Classifications follow ACMG/AMP terminology "
    "[@richards2015]. Variants classified Pathogenic, Likely pathogenic or Pathogenic/Likely "
    "pathogenic were labelled positive (P/LP), and those classified Benign, Likely benign or "
    "Benign/Likely benign negative (B/LB). Variants of uncertain significance and those with "
    "conflicting classifications were excluded from the labelled cohort. Every ClinVar-observed "
    f"residue of these proteins ({V['n_all']} variants in {V['n_all_prot']} proteins, of which "
    f"{V['n_unlab']} carry no P/LP or B/LB classification) was retained without labels as the "
    "reference set used for reconstruction and for the label-free analyses. After "
    "deduplication on gene and protein change, the labelled "
    f"cohort comprised {V['n']} variants in {V['n_genes']} genes ({V['n_p']} P/LP, {V['n_b']} "
    f"B/LB; base rate {V['base_rate']}). ClinVar variation identifiers, UniProt accessions and "
    "all features are given in Additional file 2.", "first")
d.h(3, "Position-level scores and their components")
d.p("The dynamics-aware scores are defined per residue [@karagol2026]. Site conservation *c* "
    "was computed in the EVcouplings framework from the Shannon entropy of Jackhmmer alignment "
    "columns [@hopf2019]. Mobility was the backbone root-mean-square fluctuation (RMSF) averaged "
    "over three ATLAS replicas and min–max normalised within each protein to give *m* ∈ [0, 1]. "
    "The products *c*·*m* and *c*·(1 − *m*) were then min–max rescaled within each protein to "
    "0–10 to give DCS and RCS, and the normalised DCS is log(1 + DCS)/log 11, a monotone "
    "transform of DCS. Because each step is deterministic, every protein's rescaling can be "
    f"recovered. For {V['n_maps']} proteins, residues of uncertain significance with published "
    "ingredients define linear maps from the raw products to the published scores (maximum "
    f"residual below 10^{{−9}} score units); inverting the maps returns *c*, *m* and mean RMSF "
    f"for every scored residue. The reconstruction was exact on the {V['rec_n']} labelled "
    "residues whose ingredients are available directly (maximum absolute error below "
    "10^{−9}). The slope of each protein's map, the factor by which its raw products are "
    "multiplied to reach the 0–10 scale, is called its DCS multiplier below.", "first")
d.h(3, "Substitution descriptors and comparator")
d.p("Sixteen substitution descriptors were computed from the reference and alternate residues: "
    "signed and absolute changes in Kyte–Doolittle hydrophobicity [@kyte1982], formal charge and "
    "side-chain volume; signed changes in polarity, aromaticity and flexibility; reference and "
    "alternate hydrophobicity and flexibility; and indicators for reference or alternate glycine "
    "or proline and for loss of either. Variant-specific AlphaMissense scores and AlphaMissense "
    "scores averaged over all substitutions at a residue [@cheng2023], obtained through the "
    "HegeLab resource [@tordai2024], served as external comparators and were never inputs to the "
    "dynamics-aware models. AlphaMissense is among the strongest predictors in independent "
    "assessments [@livesey2023,@rastogi2025,@mahmood2017], and ClinGen has calibrated its score "
    "intervals for PP3/BP4 evidence [@bergquist2025].", "first")
d.h(3, "Within-protein encodings and transplant")
d.p("Position-level features were re-expressed relative to a label-free reference consisting "
    "of every distinct ClinVar-observed residue of the same protein, of any classification, as "
    "mid-rank percentiles, z-scores, centred values and min–max values. Cohort-internal versions "
    "used the labelled variants of the same gene as reference. To test whether per-protein "
    "rescaling damages a feature that is known to be informative, the same transformations were "
    "applied to AlphaMissense residue means, together with multiplication of each residue mean "
    "by its protein's DCS multiplier (scaled to unit mean).", "first")
d.h(3, "Exploratory ensemble")
d.p("The ensemble that motivated this study used twelve features: DCS, RCS and mean RMSF, plus "
    "nine features selected by univariate ANOVA F statistics from a pool of dynamics, "
    "conservation and substitution descriptors after median imputation and Yeo–Johnson "
    "transformation [@yeo2000], both fitted on the complete labelled set. Five learners "
    "(XGBoost, LightGBM, CatBoost, random forest and L2-penalised logistic regression) were each "
    "tuned by 75 trials of tree-structured Parzen estimator optimisation [@akiba2019] against "
    "repeated stratified five-fold AUC on the complete labelled set, evaluated by "
    "leave-one-variant-out cross-validation, and blended with weights chosen by a further 300 "
    "trials that maximised the leave-one-variant-out AUC of the blend. Its out-of-fold "
    "predictions are distributed in Additional file 3.", "first")
d.h(3, "Estimands and the AUC decomposition")
d.p("Three quantities are distinguished. *Seen-gene discrimination* is the AUC for new variants "
    "in genes that already contribute labelled variants; leave-one-variant-out "
    "cross-validation targets it, and a gene-level prior learnt from a variant's gene-mates is "
    "legitimately informative for it. *New-gene discrimination* is the AUC for variants in "
    "genes absent from training; leave-one-gene-out cross-validation with pooled out-of-fold "
    "predictions targets it. *Within-gene discrimination* is the probability that a pathogenic "
    "variant scores above a benign variant of the same gene; it is unaffected by any component "
    "of a score that is constant within genes, and is what a residue-level feature must "
    "deliver.", "first")
d.p("For any score, each of the pathogenic–benign pairs is either within a gene or between "
    "genes, and AUC = *w*·AUC_{W} + (1 − *w*)·AUC_{B}, where *w* is the fraction of "
    "within-gene pairs, AUC_{W} the concordance among within-gene pairs with each gene weighted "
    "by its number of such pairs (a pair-weighted stratified concordance, a within-covariate "
    "summary in the spirit of the covariate-adjusted ROC [@janes2008], which instead weights "
    "strata by their number of cases) and AUC_{B} the concordance "
    "among between-gene pairs (derivation in Additional file 1). A score that is constant within "
    "genes has AUC_{W} = 0.5, so no component of a score that varies only between genes can "
    "move it; a score that is not constant within genes, such as a leave-one-out gene "
    "prevalence, can take any value. Under leave-one-gene-out validation, within-gene "
    "pairs are scored by one fitted model, whereas every between-gene pair compares predictions "
    "from two different models, so AUC_{B} also depends on the calibration of the fold-specific "
    "models [@forman2010,@airola2011].")
d.h(3, "Factorial ablation")
d.p("Four binary factors were crossed, each a pipeline stage that either sees or does not see "
    "the held-out data: the resampling unit (leave-one-variant-out or leave-one-gene-out); "
    "preprocessing (median imputation and Yeo–Johnson transformation fitted on the complete "
    "labelled set or refitted in each training fold); feature selection (ANOVA F selection of "
    "nine of 19 features, ranked on the complete labelled set — after the complete-set "
    "transformation — or within each training fold); and tuning of the L2 penalty of logistic "
    "regression (C chosen from eleven log-spaced values between 10^{−3} and 10^{2}, either by "
    "repeated stratified five-fold AUC on the complete set, as in the exploratory pipeline, or "
    "by an inner five-fold cross-validation inside each training fold, grouped by gene whenever "
    "the outer split was by gene). All 2^{4} = 16 configurations were evaluated for "
    "class-weighted logistic regression. Logistic regression with C = 1, random forests (400 "
    "trees, minimum leaf size 3) [@breiman2001] and gradient boosting (60 trees of depth 2) "
    "[@friedman2001] were evaluated under the 2^{3} design without tuning. Main effects and "
    "interactions are standard factorial contrasts.", "first")
d.h(3, "Cluster diagnostics")
d.p("Gene–class association was tested by χ² with a Monte Carlo exact p-value "
    f"({V['chi2_n_mc']:,} label permutations), because {V['frac_exp_lt5']}% of expected cell "
    "counts were below five, and summarised by the between-gene standard deviation σ_{g} of the "
    "log-odds of pathogenicity in a random-intercept logistic model fitted by maximum likelihood, "
    "integrating the random effect of each gene numerically on a dense uniform grid (verified "
    "against a grid twice as fine and twice as wide), with a profile-likelihood interval. The "
    "gene-prevalence baseline scored each variant by the pathogenic fraction of its gene in "
    "sample; with the variant itself excluded (leave-one-variant-out; variants that are the "
    "only member of their gene receive the cohort base rate); or, for leave-one-gene-out, as "
    "the pathogenic fraction of the remaining genes. The last carries no information about the "
    "held-out gene and isolates the effect of the class-balance shift that removing a gene "
    "imposes on the training data.", "first")
d.p("Gene-identity recovery asked whether features identify the gene. A random forest (150 "
    "trees) predicted the gene of each of the 76 variants in the eight genes with at least "
    "three variants under five-fold cross-validation stratified by gene and grouped by residue, "
    "so that substitutions at one residue never straddled folds, repeated five times. The "
    "balanced accuracy reported is the mean over these five repeats; it was compared with a null "
    "of 50 permutations of the gene labels, each scored by a single such cross-validation. The "
    "analysis was repeated "
    f"without labels on all {V['gid_n_res_all']} ClinVar-observed residues of the "
    f"{V['gid_n_prot_all']} proteins with at least ten.")
d.h(3, "Permutation nulls")
d.p("Labels were permuted globally, which destroys all signal, or within genes, which destroys "
    "individual-variant signal while preserving every gene's class counts exactly. The complete "
    "pipeline, including feature selection, was refitted in every replicate "
    f"({V['perm_C1_global_n']} replicates per scheme) for the leaked variant-level "
    "configuration (C1) and the controlled gene-level configuration (C8). For gene-held-out "
    "evaluation of the position-level scores alone, a third null replaced each gene's scores "
    "with residues sampled from a randomly assigned donor protein other than its own, which "
    "preserves protein-specific feature distributions while removing any relation to the labels "
    f"({V['perm_POS_donor_n']} replicates); global and within-gene nulls were also computed for "
    "the position-level scores under both resampling units. Label-independent preprocessing was "
    "cached per "
    "split and logistic models were solved by Newton's method for the same penalised objective; "
    "observed values were recomputed in the same way and agree with the scikit-learn "
    "implementation to three decimals. Permutation p-values are (*b* + 1)/(*B* + 1) and are "
    "reported with Monte Carlo standard errors.", "first")
d.h(3, "Within-gene inference")
d.p("Within-gene discrimination of each raw score was tested by stratified permutation "
    "(labels shuffled within genes; 20,000 Monte Carlo replicates), with intervals from a stratified "
    "bootstrap that resamples variants within gene and class (5,000 replicates), and by "
    "conditional logistic regression stratified by gene [@breslow1980], to which only the three "
    "genes containing both classes contribute.", "first")
d.h(3, "Incremental value over AlphaMissense")
d.p("Logistic models were fitted under leave-one-gene-out validation with AlphaMissense alone "
    "and with each candidate feature set added, selecting nine features whenever more than nine "
    "candidates were present. Paired differences were tested with the DeLong test "
    "[@delong1988], which treats variants as independent, and with Obuchowski's clustered "
    "extension, which treats genes as the independent units [@obuchowski1997]; because the "
    "cohort has few genes, the clustered statistic was referred to a *t* distribution on one "
    "fewer degree of freedom than the number of genes rather than to the normal. 95% intervals "
    "for the difference come from a paired gene-cluster bootstrap.", "first")
d.h(3, "Calibration and cross-fold comparison")
d.p("Out-of-fold predictions were summarised by the intercept and slope of the logistic "
    "recalibration of the labels on the predicted log-odds and by the Brier score. Two "
    "corrections were then applied separately within each gene-held-out fold: a prior shift, "
    "which subtracts the fold's training log-odds of pathogenicity and adds the cohort's, and a "
    "fold rank, which replaces each prediction by its percentile among the predictions the same "
    "model makes on its own training set, ties being broken by the prediction itself. Both are "
    "strictly increasing within a fold and therefore leave the within-gene component unchanged "
    "by construction. A predictor scoring every variant by the pathogenic fraction of the "
    "training genes serves as a positive control, since the prior shift must return it to 0.5. "
    "To estimate between-gene discrimination without comparing predictions across folds, a model "
    "was fitted on all but two genes and used to score both, over all pairs of genes, so that "
    "every between-gene pair was ranked by a single model that had seen neither of its genes.",
    "first")
d.h(3, "Influence, power and paired within-gene comparisons")
d.p("Each analysis was repeated with every gene removed from the cohort in turn. For a "
    "continuous score the within-gene permutation null of AUC_{W} depends only on each gene's "
    "class counts, so it was computed once (20,000 replicates, verified against an independent "
    "run) and used as the reference distribution for power calculations: scores with a specified "
    "population within-gene AUC τ were simulated in the cohort's exact gene and class structure "
    "by adding √2·Φ^{−1}(τ) to the standard normal scores of pathogenic variants, and the "
    "proportion of 4,000 replicates falling outside the two-sided 5% critical values gives the "
    "power at τ. Differences in within-gene AUC between two scores were summarised by a "
    "stratified bootstrap that resamples variants within gene and class (5,000 replicates).",
    "first")
d.h(3, "Audit of the exploratory ensemble")
d.p("The blending weights were recovered from the distributed predictions by least squares. "
    "The blend was then re-estimated with weight selection nested inside the leave-one-out "
    "loop: for each held-out variant, the weight vector that maximised the AUC of the other 87 "
    "variants was chosen from 3,006 candidates (the five vertices, equal weights and 3,000 "
    "Dirichlet(1) draws) and applied to the held-out variant. A nested choice of the single "
    "best base learner was computed in the same way.", "first")
d.h(3, "Simulation")
d.p("Cohorts were generated from *y*_{ij} ~ Bernoulli(expit(*a* + *u*_{j} + β*x*_{ij})), with "
    "gene effect *u*_{j} ~ N(0, σ_{g}^{2}), variant property *x*_{ij} ~ N(0, 1), and *a* set so "
    "that the marginal pathogenic rate was 0.41. The observed feature *x̃*_{ij} = *x*_{ij} + "
    "γ*f*_{j} carried a gene-specific offset *f*_{j} ~ N(0, 1), unrelated to *u*_{j}, that "
    "mimics a protein-specific reference frame, and a second descriptor *h̃*_{ij} = *f*′_{j} + "
    "0.3*e*_{ij} was purely gene-level. Designs comprised 12 genes with 2, 5 or 10 variants, 48 "
    "genes with 5 variants, and the observed gene-size distribution (18 genes); σ_{g} ranged "
    "from 0 to 4, β was 0 or 0.75 and γ was 0 or 1. Three predictors were evaluated: the "
    "gene-prevalence lookup, logistic regression on *x̃* and a distance-weighted "
    "seven-nearest-neighbour vote on (*x̃*, *h̃*), a local learner that can identify genes from "
    "two position-level descriptors as a random forest can. Cross-validated AUCs were compared "
    "with estimands obtained by refitting each predictor on the full simulated cohort and "
    "applying it to 40 new variants in each existing gene (seen-gene AUC) and to 400 new genes "
    "of 10 variants each (new-gene AUC and its within-gene component); 150 replicates per "
    "cell.", "first")
d.h(3, "Design requirements")
d.p("Precision was evaluated by simulation from the same random-intercept model for a fixed, "
    "gene-agnostic score with a population within-gene AUC of 0.65 or 0.75, σ_{g} of 1, 2 or 4, "
    "five variants per gene and 20 to 1,280 genes (300 replicates per design). The number of "
    "genes needed for a 95% half-width of ±0.05 or ±0.03 was obtained from a log–log "
    "least-squares fit of the sampling standard deviation of the pooled and within-gene AUCs "
    "against the number of genes.", "first")
d.h(3, "Statistics and reporting")
d.p("Discrimination was summarised by AUC and by average precision, which is more sensitive to "
    "performance on the minority class [@saito2015]. Confidence intervals for AUCs come from a "
    "gene-cluster bootstrap with 5,000 resamples [@field2007]; naive variant-level intervals "
    "are reported for comparison, and intervals for cross-validated predictions hold the "
    "out-of-fold predictions fixed. Paired AUC comparisons are two-sided. For the pipeline "
    "permutation nulls both tails are reported (Additional file 1: Table S4): the upper tail "
    "tests the pre-specified direction that leakage inflates the AUC, and the two-sided value, "
    "twice the smaller tail, is quoted wherever the direction of departure was not specified in "
    "advance, as for the position-level scores. Within-gene permutation tests of individual "
    "scores are two-sided on the absolute deviation of AUC_{W} from 0.5. Reporting "
    "follows the DOME recommendations [@walsh2021] (Additional file 1). Analyses used Python "
    "3.13 with scikit-learn 1.9.1 [@pedregosa2011], SciPy 1.18.1, statsmodels 0.15.0 "
    "[@seabold2010], NumPy 2.5.3 and pandas 3.0.5, with fixed random seeds. Every number in "
    "the text is generated by the deposited scripts and checked against the manuscript by an "
    "automated audit (Additional file 3).", "first")

# The Results, Discussion and back matter are appended by results_text.py
from results_text import write_rest   # noqa: E402
write_rest(d, V, FIG)

out_dir = HERE
d.save(os.path.join(out_dir, "manuscript.json"), os.path.join(out_dir, "manuscript.md"))
print("manuscript JSON and Markdown written;", len(d.c.order), "references cited")
