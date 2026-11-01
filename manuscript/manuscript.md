Gene-level clustering confounds cross-validated benchmarks of structure-informed missense variant classifiers

Taner Karagöl^{1,¶,*} and Alper Karagöl^{2,3,¶,*}

^{1} McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA 02139, USA

^{2} Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139, USA

^{3} Broad Institute of MIT and Harvard, Cambridge, MA 02142, USA

^{¶} These authors contributed equally. ^{*} Correspondence: karagolt@mit.edu (T.K.); akaragol@mit.edu (A.K.)

**Running head:** Gene clustering in variant classifier benchmarks

## Abstract

**Background:** Structure- and dynamics-informed features are increasingly proposed as complements to sequence-based missense pathogenicity predictors. Because clinically labelled variants are distributed unevenly across genes, the design of a cross-validated benchmark can determine what it measures. Dynamics-aware conservation scores serve as a worked case.

**Methods:** We assembled 88 ClinVar missense variants (36 pathogenic, 52 benign) in 18 genes with molecular dynamics trajectories, reconstructed the components of two position-level scores, and evaluated classifiers in a factorial design crossing the resampling unit with the placement of preprocessing, feature selection and tuning. Every area under the ROC curve (AUC) was decomposed exactly into within-gene and between-gene components and calibrated against gene-prevalence baselines, structure-preserving permutation nulls, clustered tests against AlphaMissense, fold-wise calibration diagnostics and simulations with known estimands.

**Results:** Fifteen of 18 genes contained one class, so 84% of pathogenic–benign pairs spanned genes. A lookup of each gene's pathogenic fraction reached AUC 0.739 under leave-one-variant-out validation, entirely between genes, and a tuned five-model ensemble reached 0.756, a difference not distinguishable from zero (clustered p = 0.93). Against labels permuted within genes, the leaked pipeline was not significant (p = 0.170). The position-level scores identified their protein almost as well before as after per-protein rescaling (balanced accuracy 0.25 and 0.26; chance 0.125). Their within-gene AUC (0.681, p = 0.065) matched that of their conservation factor (0.668); normalised mobility gave 0.483. Gene-held-out, the scores gave 0.152, against 0.189 for a predictor knowing only training prevalence; subtracting each fold's training prior raised this to 0.297, and ranking every between-gene pair with a model blind to both genes raised the between-gene component from 0.126 to 0.216. No addition detectably improved gene-held-out AlphaMissense (0.874; all clustered p ≥ 0.09). Simulations showed variant-level inflation and pooled gene-level deflation growing with between-gene heterogeneity, while the within-gene component tracked its estimand. The cohort reaches 80% power only against a within-gene AUC of 0.76.

**Conclusions:** On gene-clustered cohorts the pooled AUC is predominantly a between-gene statistic: variant-level validation rewards gene memorisation and pooled gene-level validation penalises gene-level information, partly correctably. Within- and between-gene components, with matched baselines, nulls and clustered inference, separate these effects. Here the dynamics component added no within-gene discrimination beyond conservation.

**Keywords:** variant effect prediction; benchmarking; data leakage; cross-validation; clustered data; within-gene concordance; AlphaMissense; protein dynamics

## Background

Computational prediction of missense variant effects is routine in clinical genomics, and published predictors far outnumber independent evaluations of them. The cohorts available for evaluation are small, unevenly distributed across genes and partly circular with respect to the predictors under test [1–3]. Grimm and colleagues showed that variants from the same gene and overlap with training data both inflate apparent accuracy, and that correcting for them reorders published rankings [1]; later assessments confirm that evaluation design accounts for much of the reported difference between tools [2,3].

Structural and dynamic descriptors are increasingly added to variant interpretation. Predictors draw on predicted structure [4,5], and features derived from molecular dynamics, conformational ensembles and coevolutionary coupling networks have been proposed [6–8]. Bahar and colleagues showed that elastic-network dynamics descriptors carry information about variant pathogenicity [9] and built them into Rhapsody [10]; its successor, Rhapsody-2, was trained on 117,525 variants in 12,094 proteins and reached an AUROC of 0.94 under cross-validation that excluded every variant of a test protein from training [11]. Holding out whole proteins is the design choice examined here.

Many structure- and dynamics-derived descriptors are position-level quantities: they describe a residue rather than a substitution, take one value for every substitution at a site, and are often normalised within each protein so that residues can be ranked along a sequence. Two properties of clinically labelled cohorts interact with such features. First, genes differ greatly in the fraction of their labelled variants that are pathogenic, and many contribute variants of one class only. Second, the area under the ROC curve (AUC), the usual summary of discrimination, averages over all pathogenic–benign pairs, including pairs drawn from different genes. Together these mean that when most genes are class-pure, most pairs lie between genes and the AUC largely measures how well a score ranks genes rather than variants. The diagnostic-testing literature separates the two through the covariate-adjusted ROC curve, which compares cases and controls sharing a covariate value [12], but the distinction is rarely drawn in variant-effect benchmarking.

Information leakage into cross-validated estimates is well documented. Fitting any part of a pipeline, whether scaling, feature selection or hyperparameter search, on data that include the evaluation fold produces optimistic estimates [13–16], and dependency in the data invalidates random splitting unless the resampling unit matches the unit of dependence [17,18]. Kapoor and Narayanan catalogued leakage across 294 papers in 17 fields [19], Whalen and colleagues described its genomic forms [20], and the DOME recommendations ask authors to state the resampling unit explicitly [21]. The converse problem has received less attention: pooled AUCs from leave-one-out designs are biased when held-out units shift the class balance of the training data [22–24], and leaving out a gene whose variants all belong to one class is an extreme instance of such a shift. Gene-level information is also not illegitimate in itself: the ACMG/AMP criteria PP2 and BP1 encode gene-level rates of pathogenic and benign missense variation [25], so a predictor that has learnt gene-level prevalence duplicates evidence that a clinical classification already uses.

A concrete case motivates what follows. Having defined a dynamics-aware decomposition of evolutionary conservation into a Rigid Conserved Score (RCS) and a Dynamic Conserved Score (DCS) [26], we tested whether the scores improve missense classification in a supervised model. An exploratory five-model ensemble reached an AUC of 0.756 under leave-one-variant-out cross-validation, a value that did not survive gene-aware validation. That discrepancy is the object of study here. We decompose the AUC into within-gene and between-gene components; separate the pipeline stages that leak information in a factorial design that includes hyperparameter tuning and ensemble blending; reconstruct the components of the position-level scores to establish what they encode; calibrate gene-held-out estimates against null distributions that preserve the cohort's gene structure; and generalise the findings by simulation with known estimands. The dynamics-aware scores serve as a worked example, and the analysis applies to any feature evaluated on a gene-clustered cohort.

## Methods

### Benchmark cohort

Missense variants were retrieved from ClinVar [27] for human proteins with replica molecular dynamics trajectories in the ATLAS database [28], as described previously [26]. Classifications follow ACMG/AMP terminology [25]. Variants classified Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic were labelled positive (P/LP), and those classified Benign, Likely benign or Benign/Likely benign negative (B/LB). Variants of uncertain significance and those with conflicting classifications were excluded from the labelled cohort. Every ClinVar-observed residue of these proteins (768 variants in 28 proteins, of which 680 carry no P/LP or B/LB classification) was retained without labels as the reference set used for reconstruction and for the label-free analyses. After deduplication on gene and protein change, the labelled cohort comprised 88 variants in 18 genes (36 P/LP, 52 B/LB; base rate 0.409). ClinVar variation identifiers, UniProt accessions and all features are given in Additional file 2.

### Position-level scores and their components

The dynamics-aware scores are defined per residue [26]. Site conservation *c* was computed in the EVcouplings framework from the Shannon entropy of Jackhmmer alignment columns [29]. Mobility was the backbone root-mean-square fluctuation (RMSF) averaged over three ATLAS replicas and min–max normalised within each protein to give *m* ∈ [0, 1]. The products *c*·*m* and *c*·(1 − *m*) were then min–max rescaled within each protein to 0–10 to give DCS and RCS, and the normalised DCS is log(1 + DCS)/log 11, a monotone transform of DCS. Because each step is deterministic, every protein's rescaling can be recovered. For 26 proteins, residues of uncertain significance with published ingredients define linear maps from the raw products to the published scores (maximum residual below 10^{−9} score units); inverting the maps returns *c*, *m* and mean RMSF for every scored residue. The reconstruction was exact on the 40 labelled residues whose ingredients are available directly (maximum absolute error below 10^{−9}). The slope of each protein's map, the factor by which its raw products are multiplied to reach the 0–10 scale, is called its DCS multiplier below.

### Substitution descriptors and comparator

Sixteen substitution descriptors were computed from the reference and alternate residues: signed and absolute changes in Kyte–Doolittle hydrophobicity [30], formal charge and side-chain volume; signed changes in polarity, aromaticity and flexibility; reference and alternate hydrophobicity and flexibility; and indicators for reference or alternate glycine or proline and for loss of either. Variant-specific AlphaMissense scores and AlphaMissense scores averaged over all substitutions at a residue [31], obtained through the HegeLab resource [32], served as external comparators and were never inputs to the dynamics-aware models. AlphaMissense is among the strongest predictors in independent assessments [2,33,34], and ClinGen has calibrated its score intervals for PP3/BP4 evidence [35].

### Within-protein encodings and transplant

Position-level features were re-expressed relative to a label-free reference consisting of every distinct ClinVar-observed residue of the same protein, of any classification, as mid-rank percentiles, z-scores, centred values and min–max values. Cohort-internal versions used the labelled variants of the same gene as reference. To test whether per-protein rescaling damages a feature that is known to be informative, the same transformations were applied to AlphaMissense residue means, together with multiplication of each residue mean by its protein's DCS multiplier (scaled to unit mean).

### Exploratory ensemble

The ensemble that motivated this study used twelve features: DCS, RCS and mean RMSF, plus nine features selected by univariate ANOVA F statistics from a pool of dynamics, conservation and substitution descriptors after median imputation and Yeo–Johnson transformation [36], both fitted on the complete labelled set. Five learners (XGBoost, LightGBM, CatBoost, random forest and L2-penalised logistic regression) were each tuned by 75 trials of tree-structured Parzen estimator optimisation [37] against repeated stratified five-fold AUC on the complete labelled set, evaluated by leave-one-variant-out cross-validation, and blended with weights chosen by a further 300 trials that maximised the leave-one-variant-out AUC of the blend. Its out-of-fold predictions are distributed in Additional file 3.

### Estimands and the AUC decomposition

Three quantities are distinguished. *Seen-gene discrimination* is the AUC for new variants in genes that already contribute labelled variants; leave-one-variant-out cross-validation targets it, and a gene-level prior learnt from a variant's gene-mates is legitimately informative for it. *New-gene discrimination* is the AUC for variants in genes absent from training; leave-one-gene-out cross-validation with pooled out-of-fold predictions targets it. *Within-gene discrimination* is the probability that a pathogenic variant scores above a benign variant of the same gene; it is unaffected by any component of a score that is constant within genes, and is what a residue-level feature must deliver.

For any score, each of the pathogenic–benign pairs is either within a gene or between genes, and AUC = *w*·AUC_{W} + (1 − *w*)·AUC_{B}, where *w* is the fraction of within-gene pairs, AUC_{W} the concordance among within-gene pairs with each gene weighted by its number of such pairs (a pair-weighted stratified concordance, a within-covariate summary in the spirit of the covariate-adjusted ROC [12], which instead weights strata by their number of cases) and AUC_{B} the concordance among between-gene pairs (derivation in Additional file 1). A score that is constant within genes has AUC_{W} = 0.5, so no component of a score that varies only between genes can move it; a score that is not constant within genes, such as a leave-one-out gene prevalence, can take any value. Under leave-one-gene-out validation, within-gene pairs are scored by one fitted model, whereas every between-gene pair compares predictions from two different models, so AUC_{B} also depends on the calibration of the fold-specific models [23,24].

### Factorial ablation

Four binary factors were crossed, each a pipeline stage that either sees or does not see the held-out data: the resampling unit (leave-one-variant-out or leave-one-gene-out); preprocessing (median imputation and Yeo–Johnson transformation fitted on the complete labelled set or refitted in each training fold); feature selection (ANOVA F selection of nine of 19 features, ranked on the complete labelled set — after the complete-set transformation — or within each training fold); and tuning of the L2 penalty of logistic regression (C chosen from eleven log-spaced values between 10^{−3} and 10^{2}, either by repeated stratified five-fold AUC on the complete set, as in the exploratory pipeline, or by an inner five-fold cross-validation inside each training fold, grouped by gene whenever the outer split was by gene). All 2^{4} = 16 configurations were evaluated for class-weighted logistic regression. Logistic regression with C = 1, random forests (400 trees, minimum leaf size 3) [38] and gradient boosting (60 trees of depth 2) [39] were evaluated under the 2^{3} design without tuning. Main effects and interactions are standard factorial contrasts.

### Cluster diagnostics

Gene–class association was tested by χ² with a Monte Carlo exact p-value (100,000 label permutations), because 86% of expected cell counts were below five, and summarised by the between-gene standard deviation σ_{g} of the log-odds of pathogenicity in a random-intercept logistic model fitted by maximum likelihood, integrating the random effect of each gene numerically on a dense uniform grid (verified against a grid twice as fine and twice as wide), with a profile-likelihood interval. The gene-prevalence baseline scored each variant by the pathogenic fraction of its gene in sample; with the variant itself excluded (leave-one-variant-out; variants that are the only member of their gene receive the cohort base rate); or, for leave-one-gene-out, as the pathogenic fraction of the remaining genes. The last carries no information about the held-out gene and isolates the effect of the class-balance shift that removing a gene imposes on the training data.

Gene-identity recovery asked whether features identify the gene. A random forest (150 trees) predicted the gene of each of the 76 variants in the eight genes with at least three variants under five-fold cross-validation stratified by gene and grouped by residue, so that substitutions at one residue never straddled folds, repeated five times. The balanced accuracy reported is the mean over these five repeats; it was compared with a null of 50 permutations of the gene labels, each scored by a single such cross-validation. The analysis was repeated without labels on all 623 ClinVar-observed residues of the 15 proteins with at least ten.

### Permutation nulls

Labels were permuted globally, which destroys all signal, or within genes, which destroys individual-variant signal while preserving every gene's class counts exactly. The complete pipeline, including feature selection, was refitted in every replicate (2,000 replicates per scheme) for the leaked variant-level configuration (C1) and the controlled gene-level configuration (C8). For gene-held-out evaluation of the position-level scores alone, a third null replaced each gene's scores with residues sampled from a randomly assigned donor protein other than its own, which preserves protein-specific feature distributions while removing any relation to the labels (1,000 replicates); global and within-gene nulls were also computed for the position-level scores under both resampling units. Label-independent preprocessing was cached per split and logistic models were solved by Newton's method for the same penalised objective; observed values were recomputed in the same way and agree with the scikit-learn implementation to three decimals. Permutation p-values are (*b* + 1)/(*B* + 1) and are reported with Monte Carlo standard errors.

### Within-gene inference

Within-gene discrimination of each raw score was tested by stratified permutation (labels shuffled within genes; 20,000 Monte Carlo replicates), with intervals from a stratified bootstrap that resamples variants within gene and class (5,000 replicates), and by conditional logistic regression stratified by gene [40], to which only the three genes containing both classes contribute.

### Incremental value over AlphaMissense

Logistic models were fitted under leave-one-gene-out validation with AlphaMissense alone and with each candidate feature set added, selecting nine features whenever more than nine candidates were present. Paired differences were tested with the DeLong test [41], which treats variants as independent, and with Obuchowski's clustered extension, which treats genes as the independent units [42]; because the cohort has few genes, the clustered statistic was referred to a *t* distribution on one fewer degree of freedom than the number of genes rather than to the normal. 95% intervals for the difference come from a paired gene-cluster bootstrap.

### Calibration and cross-fold comparison

Out-of-fold predictions were summarised by the intercept and slope of the logistic recalibration of the labels on the predicted log-odds and by the Brier score. Two corrections were then applied separately within each gene-held-out fold: a prior shift, which subtracts the fold's training log-odds of pathogenicity and adds the cohort's, and a fold rank, which replaces each prediction by its percentile among the predictions the same model makes on its own training set, ties being broken by the prediction itself. Both are strictly increasing within a fold and therefore leave the within-gene component unchanged by construction. A predictor scoring every variant by the pathogenic fraction of the training genes serves as a positive control, since the prior shift must return it to 0.5. To estimate between-gene discrimination without comparing predictions across folds, a model was fitted on all but two genes and used to score both, over all pairs of genes, so that every between-gene pair was ranked by a single model that had seen neither of its genes.

### Influence, power and paired within-gene comparisons

Each analysis was repeated with every gene removed from the cohort in turn. For a continuous score the within-gene permutation null of AUC_{W} depends only on each gene's class counts, so it was computed once (20,000 replicates, verified against an independent run) and used as the reference distribution for power calculations: scores with a specified population within-gene AUC τ were simulated in the cohort's exact gene and class structure by adding √2·Φ^{−1}(τ) to the standard normal scores of pathogenic variants, and the proportion of 4,000 replicates falling outside the two-sided 5% critical values gives the power at τ. Differences in within-gene AUC between two scores were summarised by a stratified bootstrap that resamples variants within gene and class (5,000 replicates).

### Audit of the exploratory ensemble

The blending weights were recovered from the distributed predictions by least squares. The blend was then re-estimated with weight selection nested inside the leave-one-out loop: for each held-out variant, the weight vector that maximised the AUC of the other 87 variants was chosen from 3,006 candidates (the five vertices, equal weights and 3,000 Dirichlet(1) draws) and applied to the held-out variant. A nested choice of the single best base learner was computed in the same way.

### Simulation

Cohorts were generated from *y*_{ij} ~ Bernoulli(expit(*a* + *u*_{j} + β*x*_{ij})), with gene effect *u*_{j} ~ N(0, σ_{g}^{2}), variant property *x*_{ij} ~ N(0, 1), and *a* set so that the marginal pathogenic rate was 0.41. The observed feature *x̃*_{ij} = *x*_{ij} + γ*f*_{j} carried a gene-specific offset *f*_{j} ~ N(0, 1), unrelated to *u*_{j}, that mimics a protein-specific reference frame, and a second descriptor *h̃*_{ij} = *f*′_{j} + 0.3*e*_{ij} was purely gene-level. Designs comprised 12 genes with 2, 5 or 10 variants, 48 genes with 5 variants, and the observed gene-size distribution (18 genes); σ_{g} ranged from 0 to 4, β was 0 or 0.75 and γ was 0 or 1. Three predictors were evaluated: the gene-prevalence lookup, logistic regression on *x̃* and a distance-weighted seven-nearest-neighbour vote on (*x̃*, *h̃*), a local learner that can identify genes from two position-level descriptors as a random forest can. Cross-validated AUCs were compared with estimands obtained by refitting each predictor on the full simulated cohort and applying it to 40 new variants in each existing gene (seen-gene AUC) and to 400 new genes of 10 variants each (new-gene AUC and its within-gene component); 150 replicates per cell.

### Design requirements

Precision was evaluated by simulation from the same random-intercept model for a fixed, gene-agnostic score with a population within-gene AUC of 0.65 or 0.75, σ_{g} of 1, 2 or 4, five variants per gene and 20 to 1,280 genes (300 replicates per design). The number of genes needed for a 95% half-width of ±0.05 or ±0.03 was obtained from a log–log least-squares fit of the sampling standard deviation of the pooled and within-gene AUCs against the number of genes.

### Statistics and reporting

Discrimination was summarised by AUC and by average precision, which is more sensitive to performance on the minority class [43]. Confidence intervals for AUCs come from a gene-cluster bootstrap with 5,000 resamples [44]; naive variant-level intervals are reported for comparison, and intervals for cross-validated predictions hold the out-of-fold predictions fixed. Paired AUC comparisons are two-sided. For the pipeline permutation nulls both tails are reported (Additional file 1: Table S4): the upper tail tests the pre-specified direction that leakage inflates the AUC, and the two-sided value, twice the smaller tail, is quoted wherever the direction of departure was not specified in advance, as for the position-level scores. Within-gene permutation tests of individual scores are two-sided on the absolute deviation of AUC_{W} from 0.5. Reporting follows the DOME recommendations [21] (Additional file 1). Analyses used Python 3.13 with scikit-learn 1.9.1 [45], SciPy 1.18.1, statsmodels 0.15.0 [46], NumPy 2.5.3 and pandas 3.0.5, with fixed random seeds. Every number in the text is generated by the deposited scripts and checked against the manuscript by an automated audit (Additional file 3).

## Results

### Most pathogenic–benign comparisons in the cohort are between genes

The 88 labelled variants are distributed very unevenly (Fig. 1a). PGK1 contributes 34 (39%), GBP1 contributes 13 benign variants and no pathogenic ones, and HMGCS2 contributes nine pathogenic variants and no benign ones. Of the 18 genes, 15 (holding 44 variants) contain a single class; only PGK1, GHRHR and KAT8 contain both. Gene and class are strongly associated (χ² = 44.7, df = 17, Monte Carlo p ≤ 1 × 10^{−5}; Cramér's V = 0.71), and a random-intercept model places the between-gene standard deviation of the log-odds of pathogenicity at 4.3 (profile-likelihood 95% interval 1.6–16.2; latent intraclass correlation 0.85, 0.44–0.99).

Of the 1,872 pathogenic–benign pairs over which the AUC averages, 298 (16%) lie within a gene (288 in PGK1, 6 in GHRHR and 4 in KAT8) and 1,574 (84%) lie between genes (Fig. 1b). An AUC computed on this cohort is therefore predominantly a statement about how a score orders genes. The two components separate scores that pooled AUCs conflate (Fig. 1c). AlphaMissense discriminates within genes as well as between them (within-gene AUC 0.896 against 0.924 between genes); the position-level scores discriminate only within genes (DCS 0.681 within, 0.511 between); and the gene-prevalence lookup discriminates only between genes (0.966), which is enough for an in-sample AUC of 0.892.

This concentration is not peculiar to the present cohort. It follows from requiring both a ClinVar classification and an existing molecular dynamics trajectory, and similar concentration is a general property of ClinVar-derived benchmarks, in which a few intensively studied genes dominate and classification quality varies with how heavily a gene has been studied [1,47,48].

### A gene-prevalence lookup matches a tuned ensemble under variant-level validation

A classifier that ignores variants and scores each by the pathogenic fraction of its gene reaches AUC 0.892 in sample and 0.739 when the scored variant is excluded from its own gene's fraction (Fig. 2a; Table 1). Excluding the variant reverses the ordering inside every mixed gene, so the within-gene component of this baseline is exactly 0.000; all of its performance is between-gene (0.879).

The exploratory ensemble reached 0.756 under the same design. Its weights, recovered exactly from the distributed predictions, were 57.0% logistic regression, 30.2% LightGBM, 5.8% XGBoost, 5.7% random forest and 1.2% CatBoost; the five base learners alone gave 0.701–0.731 and an equal-weight blend 0.723. When the weights were chosen without the held-out variant, the blend gave 0.732, and a nested choice of the single best learner gave 0.685. Choosing the weights without the held-out variant costs 0.021 against the best in-sample choice from the same candidate set (0.753); the reported 0.756 exceeds that best candidate, so the distributed weights are not reproducible from this set. The ensemble exceeded the lookup by only 0.017, a difference the data cannot distinguish from zero (clustered standard error 0.203; clustered p = 0.93; DeLong p = 0.80). Unlike the lookup, whose within-gene component is 0.000 because excluding a variant from its own gene's fraction inverts the ordering inside each gene, the ensemble retains a within-gene component of 0.718. That advantage applies to the 16% of pairs that lie within a gene; over the remaining 84% the lookup is the better of the two (0.879 against 0.763), and the pooled AUC is the pair-weighted average of the two.

### The resampling unit acts on the between-gene component

Four pipeline stages were crossed with the resampling unit. At a fixed penalty (C = 1), across three model families (Additional file 1: Table S5), the fully leaked configuration — variant-level split with preprocessing and feature selection fitted on the complete labelled set — gave 0.680 for logistic regression, and the fully controlled configuration, with a gene-level split and nothing fitted outside the training fold, gave 0.594, a difference of 0.086. Random forests (0.648 to 0.607) and gradient boosting (0.615 to 0.586) moved in the same direction, and for every model family every variant-level configuration exceeded every gene-level one.

Adding the placement of penalty tuning as a fourth factor gives the 16 logistic-regression configurations of Table 2 and Fig. 2b. The size of the contrast depends on how the penalty is tuned. In the only configuration in which no stage sees the held-out gene (L16, with C tuned by a gene-grouped inner cross-validation inside each fold), the AUC was 0.664 (median C = 0.1), whereas tuning on the complete set by variant-level folds, as in the exploratory pipeline, chose C = 0.316 and gave 0.597. Averaged over the other factors, variant-level splitting raised the AUC by 0.043, leaked feature selection by 0.018, and leaked tuning changed it by −0.012; the largest interaction was between the split and leaked selection (−0.028), because selection on the complete set matters most when the split is otherwise correct (Additional file 1: Table S6).

The decomposition shows where these differences arise. Across the eight gene-level logistic configurations the within-gene component stayed within 0.564–0.607, whereas the between-gene component spanned 0.597–0.677 and carried almost all of the variation in pooled AUC (0.592–0.666). Changing the resampling unit acts chiefly on the between-gene component — where gene memorisation operates under variant-level splitting and cross-fold calibration under gene-level splitting — but it does not leave the within-gene component untouched: for logistic regression the variant-level configurations raise it (0.631–0.701 against 0.564–0.607 across the gene-level configurations). The resampling unit is the cause, not the leaked preprocessing: averaged over the other factors the split moves AUC_{W} by +0.078, whereas scaling, selection and tuning on the complete set move it by −0.003, +0.005 and −0.027. Under leave-one-variant-out the model is shown the held-out variant's gene-mates and their labels, which is itself gene-level information acting on the within-gene ordering. For random forests the within-gene component is instead slightly higher under gene-level splitting (0.678–0.698 against 0.627–0.644).

### Permutation nulls that preserve gene structure

Permuting labels globally destroys both individual and gene-level signal; permuting them within genes destroys individual signal while holding every gene's class counts fixed (Fig. 2c; Table 1; Additional file 1: Table S4). Over 2,000 replicates of the leaked variant-level pipeline, the global null had mean 0.554 (95% range 0.410–0.686), against which the observed 0.680 gives p = 0.030 (Monte Carlo SE 0.004). The within-gene null had mean 0.619 (0.487–0.736), against which the observed value is not significant (p = 0.170; SE 0.008). A pipeline given labels that carry no individual-variant information at all thus returns 0.619 under variant-level resampling, and the observed 0.680 does not depart significantly from it.

The same nulls behave differently under gene-level resampling. For the controlled configuration C8, the within-gene null had mean 0.466, below chance, and the global null 0.492; the observed 0.594 gives p = 0.086 and p = 0.112, respectively. Neither design provides significant evidence of individual-variant discrimination by the 19-feature pipeline, and under gene-level resampling the reference point for 'no signal' is not 0.5.

### What the position-level scores encode

Reconstruction makes the per-protein rescaling explicit (Fig. 3a). The DCS multiplier ranged from 10.0 to 31.7 across the 18 cohort proteins, so identical products of conservation and mobility map to scores that differ up to 3.2-fold between proteins. The multiplier is set by each protein's most extreme residue: mobility is normalised by the protein's range of RMSF, whose maximum lay between 3.1 and 18.7 Å, and in 11 of 18 proteins no ClinVar-observed residue reached half of that range (median maximum normalised mobility 0.42).

Whether rescaling is what makes the scores protein-specific was tested by asking a classifier to recover the gene from each encoding (Fig. 3b; Additional file 1: Table S8). The published scores identified the gene with balanced accuracy 0.26 (chance 0.125; permutation null 0.13, p = 0.020, the smallest attainable with 50 permutations), and the same products before rescaling did so almost as well (0.25). Conservation and normalised mobility gave 0.20 and 0.21 individually, and absolute RMSF 0.20. Only the within-protein percentile, which gives every protein the same marginal distribution, removed identity altogether (0.07). The pattern was the same without labels across 623 ClinVar-observed residues of the 15 proteins with at least ten (published 0.14, before rescaling 0.13, percentile 0.02; chance 0.067). Rescaling therefore adds a protein signature, as multiplying AlphaMissense residue means by the same multipliers shows (balanced accuracy 0.19 raw, 0.26 multiplied), but most of the identity is present before rescaling: it is carried by the joint distribution of conservation and mobility, which differs between proteins whatever the scale, and min–max rescaling, which equalises only two order statistics, cannot remove it.

Within genes, DCS and RCS ranked pathogenic above benign variants in 0.681 and 0.668 of the 298 within-gene pairs (stratified permutation p = 0.065 and 0.088; Fig. 3c; Additional file 1: Table S11). This signal is that of the conservation factor (0.668, p = 0.088): normalised mobility alone gave 0.483 (p = 0.87) and absolute RMSF 0.483. It came almost entirely from PGK1 (DCS 0.70, conservation 0.69), was reversed in GHRHR and at chance in KAT8, and was not supported by conditional logistic regression (DCS odds ratio per standard deviation 0.84, 95% CI 0.37–1.93; conservation 1.13, 0.60–2.13), although that statistic is sensitive to the skew of the DCS scale: the normalised DCS, a monotone transform with the same within-gene AUC, gives 1.04, and a rank-based model 1.12 per decile. AlphaMissense, by contrast, discriminated within genes unambiguously (variant score 0.896, p < 0.0001; odds ratio 6.47, 2.16–19.39), as did its residue mean (0.772, p = 0.004).

### Below-chance gene-held-out AUCs reflect cohort structure

Evaluated gene-held-out, a model using only DCS, RCS and normalised DCS gave a pooled AUC of 0.152, a value that would conventionally be read as strong anti-prediction. Three reference points calibrate it (Fig. 4b). A predictor that knows nothing except the pathogenic rate of the training genes scores 0.189 under the same design, because removing a class-pure gene shifts the class balance of the training set against that gene. Refitting the pipeline after permuting labels within genes, which keeps each gene's class counts and the gene-level relationship between features and labels, gives a null centred on 0.291 (95% range 0.171–0.480). Replacing each gene's scores with residues of a randomly chosen donor protein, which keeps the gene–class structure but breaks the link between the real feature values and the labels, gives 0.471, close to the global null 0.477. Random protein fingerprints therefore do not depress the pooled estimate; the depression below chance comes from the actual position-level features' own relationship to the labels across the cohort's genes.

The observed value nevertheless lies below the within-gene null (two-sided p = 0.002), and the decomposition shows that both components contribute. The within-gene component (0.289, against a null expectation of 0.498) shows that the direction learnt from the other 17 genes reverses the ordering inside PGK1, the only gene with appreciable within-gene information (within-PGK1 AUC of DCS 0.70). The between-gene component (0.126, against 0.252) shows the converse: the direction learnt with PGK1 in training runs against the ordering of the remaining genes by prevalence. Both are expressions of one reversal between the within-gene association in PGK1 and the across-gene association in the rest of the cohort, the variant-level analogue of an ecological fallacy [49], superimposed on a structural deficit that the cohort imposes on any gene-informative predictor.

The same decomposition clarifies the alternative encodings (Fig. 4a; Additional file 1: Table S9). Before rescaling, the products gave 0.209; conservation alone gave 0.230 and normalised mobility 0.360. Within-protein percentiles gave 0.254 with a label-free reference and 0.336 when ranked among the labelled variants of each gene. For AlphaMissense residue means, by contrast, per-protein min–max rescaling (0.847) and multiplication by the published DCS multipliers (0.858) left transfer intact relative to the raw scale (0.845), whereas within-protein ranking lowered it to 0.792, because it discards genuine between-protein differences in constraint (between-gene component 0.858 raw against 0.795 ranked). Rescaling does not damage a feature that carries between-gene information; it cannot create information that is not there. The deficit is also not specific to gene-level splitting: under leave-one-variant-out validation the position-level scores gave 0.278, inside the range of their within-gene null (mean 0.387; two-sided p = 0.324), because in a cohort of this structure leave-one-out predictions from weak features are anti-correlated with the held-out labels [22,23].

### Cross-fold calibration accounts for part of the gene-held-out deficit

Because every between-gene pair compares predictions from two different fitted models, AUC_{B} depends on how those models are calibrated relative to one another, and the calibration of the out-of-fold predictions shows that they are not (Fig. 6a; Additional file 1: Table S15). AlphaMissense is well calibrated under either design (recalibration slope 0.91 gene-held-out, 1.05 variant-level, against 1 for perfect calibration). The 19-feature model is overconfident (0.28 against 0.50), and the position-level model has a negative slope (−5.53), the calibration signature of a predictor whose ordering is reversed.

The source of the miscalibration is the class balance of each training set. Across the 18 gene-held-out folds, the mean prediction tracks the prevalence of the training genes (*r* = 0.74) and therefore runs against the prevalence of the gene being predicted (*r* = −0.45), because removing a gene shifts the training balance away from that gene's own composition (*r* = −0.63). HMGCS2 is the extreme case: all nine of its variants are pathogenic, so holding it out leaves the lowest training prevalence of any fold (0.34) and produces the lowest mean prediction (0.30) for the gene that should score highest.

Two corrections test how much of the deficit this explains (Fig. 6b; Additional file 1: Table S16). Both are applied separately within each fold and are strictly increasing, so neither can change AUC_{W}. Subtracting each fold's training prior from its predicted log-odds raises the pooled gene-held-out AUC of the position-level scores from 0.152 to 0.297, entirely through the between-gene component (0.126 to 0.299); applied to a predictor that knows only the training prevalence, the same correction returns exactly 0.500 from 0.189, as it must. Replacing each prediction by its percentile among the training predictions of the same fold leaves the estimate unchanged (0.150), so what separates the folds is the location shift imposed by the prior and not a difference of scale.

A design that never compares predictions across folds points the same way. Fitting a model on all but two genes and using it to score both, repeated over all 153 pairs of genes so that each of the 1,574 between-gene pairs is ranked by a single model that saw neither of its genes, raises the between-gene component of the position-level scores from 0.126 to 0.216 (+0.090) and that of AlphaMissense from 0.870 to 0.924; for the 19-feature model it moves the other way (0.602 to 0.554). The prior-shift correction recovers two fifths of the distance between the pooled gene-held-out estimate and 0.5, and the leave-two-genes-out design a quarter of the corresponding distance for the between-gene component. What the correction leaves is a gap of 0.203 to 0.5, of which the still-depressed between-gene component contributes 0.169 and the within-gene reversal 0.034; the deficit is therefore mostly, but not only, a between-gene phenomenon. Leave-pair-out designs remove the artefact at the cost of one fit per pair of held-out units [23], which is affordable at this scale but not for a proteome-wide benchmark, where the prior-shift correction is the practical option.

### No feature set detectably improves on AlphaMissense

AlphaMissense alone reached 0.874 in gene-held-out logistic regression (Fig. 4c; Additional file 1: Table S10). Adding the position-level scores as published lowered this to 0.854 (−0.020). The decrease is borderline: it is significant when variants are treated as independent (DeLong p = 0.047) and its paired gene-cluster bootstrap interval excludes zero (−0.100 to −0.002), but Obuchowski's clustered test, which treats the 18 genes as the units, does not reject (p = 0.129). Conservation alone behaved the same way (−0.018; clustered p = 0.093). No addition detectably improved on AlphaMissense: within-protein percentiles gave −0.019 (clustered p = 0.464), within-gene percentiles among the labelled variants −0.022 (0.371), substitution descriptors +0.012 (0.405) and the full feature set −0.001 (0.944).

The percentile result depends on a detail that illustrates the general problem. Percentiles here are mid-ranks, (rank − 1/2)/*n*, whose mean is 1/2 in every gene. The common convention rank/*n* has mean (*n* + 1)/2*n*, which is 1 for a gene with one variant and approaches 1/2 only for large genes, so it encodes gene size; with that convention AlphaMissense plus within-gene percentiles reached 0.892, an apparent gain of +0.018 that disappears when the encoding no longer carries gene-level information.

### Simulation separates inflation, deflation and the within-gene estimand

In simulated cohorts with known estimands (Fig. 5; Additional file 1: Table S13), a nearest-neighbour learner that can recognise genes from two position-level descriptors was evaluated as the cohort was. Its variant-level AUC exceeded the new-gene AUC by an amount that grew with between-gene heterogeneity and with the number of variants per gene: with twelve genes of ten variants the excess was −0.005 without heterogeneity, +0.103 at a between-gene SD of 2 and +0.174 at 4, against +0.065 at SD 2 with five variants per gene and +0.107 on the observed gene-size distribution. Relative to the seen-gene AUC, which it is designed to estimate, the variant-level estimate was close to unbiased (−0.016 at SD 2 with ten variants per gene): variant-level validation answers a legitimate question, but not the question asked of a new gene. On uneven gene-size distributions the variant-level and seen-gene AUCs diverge even for this learner, because leave-one-variant-out weights genes by size whereas the seen-gene estimand weights them equally. A logistic model on a single shifted feature, which cannot recognise genes, showed little inflation (+0.024 at SD 4).

Pooled gene-level estimates were biased in the opposite direction (−0.042 at SD 2 and −0.059 at SD 4 with twelve genes of ten variants; −0.095 on the observed gene sizes at SD 4), and the bias shrank as genes were added (−0.016 with 48 genes at SD 4). For both the nearest-neighbour and the logistic learner, the within-gene component of the gene-held-out AUC tracked its estimand (largest absolute mean bias 0.033 across the designs with five or ten variants in every gene, 0.042 on the observed gene-size distribution and 0.075 with two variants per gene). Without any individual signal, the gene-prevalence lookup on the observed gene sizes scored 0.693–0.779 under variant-level and 0.144–0.204 under gene-level validation at between-gene SDs of 2–3, bracketing the values observed in the real cohort (0.739 and 0.189) and consistent with the interval 1.6–16.2 estimated from the random-intercept model.

### Design requirements

Because only genes that contain both classes contribute within-gene pairs, the number of such genes, not of variants, limits precision (Fig. 5c; Additional file 1: Table S14). With five variants per gene, about 62% of genes contain both classes at a between-gene SD of 2 and 38% at an SD of 4. Estimating a within-gene AUC of 0.65 to ±0.05 then requires about 220 genes at SD 2 and 378 at SD 4 (607 and 1,012 for ±0.03), compared with 172 at SD 1; a within-gene AUC of 0.75 at SD 4 requires 270. The pooled AUC of the same score is estimated more precisely (103 and 110 genes for ±0.05) but estimates a different quantity: between-gene heterogeneity attenuates it towards 0.5, to 0.558 for a score whose within-gene AUC is 0.65 and 0.607 for one whose within-gene AUC is 0.75 at SD 4. The present cohort, with 18 genes of which 3 contain both classes, falls far short of these requirements.

### What this cohort can and cannot resolve

The limits can be stated exactly for the cohort as it stands. Because the within-gene permutation null of AUC_{W} depends only on each gene's class counts and not on the score, it can be computed once and used to obtain the power of the test against scores of known within-gene discrimination in this cohort's exact structure (Fig. 6c; Additional file 1: Table S18). The test holds its size (0.054 at AUC_{W} = 0.5) and its null occupies 0.305–0.691. Power reaches 80% only at a within-gene AUC of 0.76; a score with a true within-gene AUC of 0.65 is detected 0.35 of the time and one with 0.70, 0.55. The within-gene AUC observed for DCS (0.681) lies below that threshold, so its non-significance is a statement about the cohort and not about the score.

Direct comparisons between scores behave accordingly (Additional file 1: Table S18). AlphaMissense exceeds DCS within genes by +0.215 (stratified bootstrap 95% interval +0.007 to +0.419) and conservation by +0.228 (+0.081 to +0.406), both intervals excluding zero. The comparisons that matter for the dynamics question do not: DCS exceeds conservation by +0.013 (−0.164 to +0.201), and normalised mobility falls below it by 0.185 (−0.480 to +0.141). The cohort therefore resolves the gap between AlphaMissense and the dynamics-aware scores but not the gap between those scores and their own conservation factor.

Every within-gene estimate here is, moreover, an estimate about one gene. Removing PGK1 leaves 10 of the 298 within-gene pairs, and the within-gene AUCs collapse (DCS 0.681 to 0.20, conservation 0.668 to 0.00, AlphaMissense 0.896 to 0.80). The pooled gene-held-out AUC is similarly concentrated: over the 18 single-gene deletions it ranges from 0.141 to 0.548, and the upper end is reached by deleting HMGCS2, the one gene whose 9 variants are all pathogenic (between-gene component 0.558; Additional file 1: Table S17). A single class-pure gene thus accounts for most of the departure from chance.

## Discussion

On a gene-clustered cohort the pooled AUC is predominantly a between-gene statistic. Here 84% of the pathogenic–benign pairs compare variants from different genes and 15 of 18 genes contribute a single class, and the other findings follow from that fact. Variant-level validation lets a model reach a high AUC by recognising each variant's gene and applying its prevalence, so that a lookup using no variant information matches a tuned five-model ensemble. Gene-level validation removes that route but brings its own artefact: between-gene pairs then compare predictions from different fitted models, and removing a class-pure gene shifts the training class balance against it, so the pooled estimate is pulled below its target, to 0.189 for a predictor that knows only the training prevalence. The within-gene component isolates residue-level discrimination from both of these between-gene artefacts, because a score that is constant within genes has AUC_{W} = 0.5 whatever the design. It is not, however, invariant to the resampling unit once a model is fitted: under leave-one-variant-out the model is shown the held-out variant's gene-mates and their labels, which raised AUC_{W} by 0.078 here, so the within-gene component should be read from gene-held-out predictions.

For the dynamics-aware scores, the within-gene evidence is modest (AUC 0.681, p = 0.065, almost entirely from PGK1) and is the evidence carried by the conservation factor; normalised mobility shows none. This agrees with the original description of the scores, in which RCS tracked conservation against AlphaMissense classes while DCS did not [26]. It does not show that dynamics is irrelevant to variant interpretation: Rhapsody-2's protein-held-out performance and the enrichment of disease mutations at hinge-adjacent residues are evidence to the contrary [9,11,50], and three genes with both classes cannot test a mechanism. What it shows is that this cohort cannot detect a dynamics contribution, and that a validation design resting on variant-level resampling measures gene identity. The result also underlines how strong the sequence-based baseline has become: AlphaMissense discriminated within genes as well as between them, and no addition detectably improved its gene-held-out performance (the largest change, −0.020 from the position-level scores, was significant by the variant-level DeLong test and its bootstrap interval but not by the clustered test).

Per-protein normalisation is a modelling decision, but it was not the decision that caused the problems here. It multiplies each protein's raw values by a factor set by that protein's most extreme residue (10.0 to 31.7 for DCS in this cohort) and so distorts cross-protein comparisons, yet the products before rescaling identify proteins almost as well, and applying the published multipliers to AlphaMissense added protein identity without degrading its transfer to held-out genes. Protein identity is carried by the joint distribution of conservation and mobility, which differs between proteins whatever the scaling. Only an encoding that equalises the whole marginal distribution, the within-protein rank, removes it, and removal has a cost: ranking AlphaMissense residue means within proteins lowered their gene-held-out AUC from 0.845 to 0.792 because between-protein differences in constraint are real signal. The reference frame of a feature should therefore follow the estimand: within-protein ranks when within-gene discrimination is the target, proteome-calibrated scales when between-gene information is wanted.

The problem differs from the circularity described by Grimm and colleagues [1], which concerns overlap between the data used to train a published predictor and the data used to evaluate it. Here the training set is the evaluation cohort itself; their type 1 circularity, variants of one gene on both sides of the split, is the closest analogue, with the addition that the gene relationship is legible from the features as well as from the labels. In the taxonomy of Kapoor and Narayanan [19] it lies between the absence of an independent test set and the use of illegitimate features. The complementary bias of pooled gene-held-out estimates is the cluster-level counterpart of the stratification bias described for leave-one-out designs [22–24]. Gene-level prevalence is not illegitimate information in itself, since PP2 and BP1 encode it [25], but a predictor that has learnt it should not be combined with those criteria without accounting for the overlap, and its benchmark performance says little about residue-level discrimination.

These observations translate into a short protocol (Table 3). The composition of pathogenic–benign pairs and the within/between decomposition should accompany any pooled AUC on a clustered cohort. A gene-prevalence baseline, evaluated under the same design, shows how much of a reported value gene identity alone can produce, and a within-gene permutation null tests whether anything remains. Gene-level resampling, with feature selection and tuning nested and grouped by gene, should be the default whenever features describe positions or proteins, with the within-gene component reported beside the pooled AUC, because the pooled gene-held-out value has its own downward bias. That bias is partly removable: subtracting each fold's training prior from its predicted log-odds costs nothing, cannot affect the within-gene component, and recovered two fifths of the deficit here, while leave-pair-out designs remove the cross-fold comparison entirely where the number of held-out units makes them affordable. A null result should be accompanied by the effect it could have detected, which for a within-gene comparison follows from the class counts alone. Intervals and paired tests should treat genes as the independent units [42,44]; here the clustered and variant-level tests disagreed on four of the thirteen paired comparisons examined, and the clustered test was the more significant of the two in three of them, so it is not uniformly conservative (Additional file 1: Table S2). For position-level or protein-derived features, a gene-identity recovery test and an explicit statement of the reference frame belong in the methods. Benchmarks should be planned by counting genes that contain both classes, the quantity that limits within-gene precision. Each of these steps is inexpensive, requires no additional data, and together they separate residue-level signal from gene-level structure before a pooled AUC is interpreted.

The limitations are those of the cohort. Only three genes contain both classes, so within-gene estimates are dominated by PGK1 and every interval is wide. ClinVar classifications draw on computational evidence that includes conservation [1,51], so conservation-derived features are evaluated optimistically. The reconstruction depends on residues of uncertain significance in the same proteins, which were available for every labelled protein here but need not be in general. The audit of the exploratory ensemble could re-estimate only its blending step, because some of its input features are not reconstructible from the published scores. The simulations establish generality only within the generative model specified. Deep mutational scanning resources [52–54] evaluate predictors within proteins and so avoid between-gene confounding by construction, although assay scores and clinical classifications correspond imperfectly [55]; a replication there, and in ClinVar-scale cohorts with gene-held-out designs, is the natural next step. The present analysis is a precondition for such work rather than a substitute for it: performance reported on clustered cohorts is only as informative as the baseline it is compared with [56,57].

## Conclusions

On a gene-clustered ClinVar cohort, the AUC of a structure-informed classifier was dominated by comparisons between genes. Variant-level cross-validation rewarded gene memorisation, so that a prevalence lookup (0.739) matched a tuned ensemble (0.756), and pooled gene-held-out validation penalised gene-level information, producing below-chance values (0.152) that are largely explained by the cohort's gene structure. Decomposing the AUC into within-gene and between-gene components, together with a gene-prevalence baseline, a within-gene permutation null and clustered inference, separates these effects. Applied here, they show that the dynamics-aware scores carry no within-gene information beyond that of their conservation factor and no detectable improvement on AlphaMissense, and that resolving an incremental benefit in either direction requires cohorts with many genes that contain both classes.

## Abbreviations

ACMG/AMP: American College of Medical Genetics and Genomics/Association for Molecular Pathology; AP: average precision; AUC: area under the receiver operating characteristic curve; B/LB: benign or likely benign; DCS: Dynamic Conserved Score; LOO: leave-one-out; MD: molecular dynamics; P/LP: pathogenic or likely pathogenic; RCS: Rigid Conserved Score; RMSF: root-mean-square fluctuation; ROC: receiver operating characteristic.

## Declarations

**Ethics approval and consent to participate.** Not applicable. This computational study used publicly available data and did not involve human participants, animal subjects or identifiable data.

**Consent for publication.** Not applicable.

**Availability of data and materials.** ClinVar classifications are available at https://www.ncbi.nlm.nih.gov/clinvar/. Molecular dynamics trajectories are available from ATLAS (https://www.dsimb.inserm.fr/ATLAS). AlphaMissense scores are available at https://alphamissense.hegelab.org. The dynamics-aware scores were computed as described previously [26] and can be regenerated with ADEPT (https://www.karagolresearch.com/adept). The labelled cohort with ClinVar variation identifiers, reconstructed components and all features is provided as Additional file 2. All analysis code, intermediate results, out-of-fold predictions, permutation distributions, simulation outputs and the consistency audit are provided as Additional file 3 and are deposited at https://github.com/karagol-taner/Dynamics-aware-Evolutionary-Profiling.

**Competing interests.** The authors declare no competing interests.

**Funding.** The authors received no specific funding for this work.

**Authors' contributions.** T.K. and A.K. contributed equally. Conceptualisation, T.K., A.K.; methodology, T.K., A.K.; software, T.K., A.K.; formal analysis, T.K., A.K.; investigation, T.K., A.K.; visualisation, T.K., A.K.; writing, original draft, T.K., A.K.; writing, review and editing, T.K., A.K. Both authors read and approved the final manuscript.

**Acknowledgements.** Not applicable.

## Tables

**Table 1. Baselines, null distributions and comparators.** Values are AUCs with their within-gene (AUC_{W}) and between-gene (AUC_{B}) components; null rows give the mean over replicates.

Baseline or control | Resampling | AUC | AUC_{W} | AUC_{B} | Note

Gene prevalence, in sample | — | 0.892 | 0.500 | 0.966 | Uses no variant-level information

Gene prevalence, variant excluded | Variant | 0.739 | 0.000 | 0.879 | Cluster 95% CI 0.340–0.980

Gene prevalence, gene excluded | Gene | 0.189 | 0.500 | 0.130 | Training prevalence only; prior-shift reference

Exploratory ensemble, as reported | Variant | 0.756 | 0.718 | 0.763 | Cluster 95% CI 0.547–0.865

Exploratory ensemble, nested blend | Variant | 0.732 | 0.685 | 0.741 | Weights chosen without the held-out variant

C1 pipeline, labels permuted globally | Variant | 0.554 | 0.550 | 0.555 | Observed 0.680: p = 0.030

C1 pipeline, labels permuted within gene | Variant | 0.619 | 0.535 | 0.635 | Observed 0.680: p = 0.170

C8 pipeline, labels permuted within gene | Gene | 0.466 | 0.502 | 0.459 | Observed 0.594: p = 0.086

Position-level scores | Gene | 0.152 | 0.289 | 0.126 | DCS, RCS, normalised DCS

  same, labels permuted within gene | Gene | 0.291 | 0.498 | 0.252 | Observed below null: two-sided p = 0.002

  same, donor-protein scores | Gene | 0.471 | 0.499 | 0.466 | Protein fingerprint without signal

Substitution descriptors | Gene | 0.635 | 0.607 | 0.640 | 16 descriptors, nine selected in fold

AlphaMissense, raw score | — | 0.919 | 0.896 | 0.924 | Cluster 95% CI 0.858–0.989; no model fitted

AlphaMissense, logistic model | Gene | 0.874 | 0.896 | 0.870 | Matched comparator for incremental value

AlphaMissense + position-level scores | Gene | 0.854 | 0.889 | 0.848 | DeLong p = 0.047; clustered p = 0.129

**Table 2. Logistic regression under every combination of four validation choices.** All configurations use the same 19 features (three position-level scores and 16 substitution descriptors) on the same 88 variants in 18 genes. *Split*: leave-one-variant-out (Variant) or leave-one-gene-out (Gene). *Scaling*, *selection* and *tuning*: median imputation with Yeo–Johnson transformation, ANOVA F selection of nine features, and choice of the L2 penalty C, each fitted on the full labelled set or inside every training fold. AUC_{W} and AUC_{B}: within-gene and between-gene components of the AUC.

Config. | Split | Scaling | Selection | Tuning | AUC | AP | AUC_{W} | AUC_{B} | Median C

L1 | Variant | full set | full set | full set | 0.685 | 0.599 | 0.631 | 0.695 | 0.316

L2 | Variant | full set | full set | in fold | 0.673 | 0.578 | 0.685 | 0.671 | 0.0316

L3 | Variant | full set | in fold | full set | 0.691 | 0.586 | 0.681 | 0.693 | 0.316

L4 | Variant | full set | in fold | in fold | 0.696 | 0.562 | 0.701 | 0.694 | 0.316

L5 | Variant | in fold | full set | full set | 0.684 | 0.598 | 0.631 | 0.694 | 0.316

L6 | Variant | in fold | full set | in fold | 0.677 | 0.595 | 0.701 | 0.673 | 0.0316

L7 | Variant | in fold | in fold | full set | 0.676 | 0.574 | 0.634 | 0.684 | 0.316

L8 | Variant | in fold | in fold | in fold | 0.696 | 0.566 | 0.698 | 0.695 | 0.316

L9 | Gene | full set | full set | full set | 0.666 | 0.587 | 0.607 | 0.677 | 0.316

L10 | Gene | full set | full set | in fold | 0.661 | 0.525 | 0.604 | 0.672 | 0.1

L11 | Gene | full set | in fold | full set | 0.592 | 0.512 | 0.564 | 0.597 | 0.316

L12 | Gene | full set | in fold | in fold | 0.620 | 0.513 | 0.567 | 0.630 | 0.316

L13 | Gene | in fold | full set | full set | 0.663 | 0.587 | 0.607 | 0.673 | 0.316

L14 | Gene | in fold | full set | in fold | 0.665 | 0.578 | 0.607 | 0.676 | 0.1

L15 | Gene | in fold | in fold | full set | 0.597 | 0.517 | 0.587 | 0.599 | 0.316

L16 | Gene | in fold | in fold | in fold | 0.664 | 0.518 | 0.597 | 0.677 | 0.1

L16 is the only configuration in which no stage sees the held-out gene. The 2^{3} design without tuning for logistic regression (C = 1), random forest and gradient boosting is given in Additional file 1: Table S5; factorial effects in Table S6. Main effects (variant minus gene, full set minus in fold): split +0.043, scaling −0.005, selection +0.018, tuning −0.012.

**Table 3. Diagnostic protocol for benchmarks on gene-clustered cohorts.**

Diagnostic | Question | Computation | This cohort

Pair composition | What fraction of the AUC compares variants of different genes? | Count pathogenic–benign pairs within and between genes | 84% between genes; 3 of 18 genes mixed

AUC decomposition | Does a score discriminate within genes, between them, or both? | AUC = w·AUC_{W} + (1 − w)·AUC_{B} | Ensemble 0.718 within, 0.763 between

Gene-prevalence baseline | How much can gene identity alone produce? | Score each variant by its gene's pathogenic fraction under the evaluated design | 0.739 (variant-level); 0.189 (gene-level)

Within-gene permutation null | Is there signal beyond gene-level structure? | Permute labels within genes; refit the full pipeline | p = 0.170 (variant-level pipeline)

Nested, gene-grouped selection and tuning | Do preprocessing, selection, tuning or blending see the held-out data? | Refit every data-dependent step inside each training fold | Blending: 0.756 → 0.732

Clustered inference | Are intervals and paired tests anticonservative? | Gene-cluster bootstrap; Obuchowski test on gene-level df | AM vs AM + scores: DeLong p = 0.047, clustered p = 0.129

Fold calibration | Are the fold-specific models on a common scale? | Calibration slope per design; subtract each fold's training prior | Slope -5.53; AUC_{B} 0.126 → 0.299

Detectable effect | Is a null result informative? | Power of the within-gene test against known AUC_{W}, at the observed class counts | 80% power only at AUC_{W} = 0.76

Gene-identity recovery | Do the features encode the protein? | Predict gene from features, residue-grouped CV, balanced accuracy | 0.26 (chance 0.125); ranks 0.07

Mixed-gene count | Can the cohort estimate within-gene discrimination? | Count genes with both classes; plan by simulation | 3 of 18 contain both classes; hundreds of mixed genes needed

## References

1. Grimm DG, Azencott CA, Aicheler F, Gieraths U, MacArthur DG, Samocha KE, et al. The evaluation of tools used to predict the impact of missense variants is hindered by two types of circularity. *Hum Mutat*. 2015;36(5):513–23.

2. Livesey BJ, Marsh JA. Updated benchmarking of variant effect predictors using deep mutational scanning. *Mol Syst Biol*. 2023;19(8):e11474.

3. Gunning AC, Fryer V, Fasham J, Crosby AH, Ellard S, Baple EL, et al. Assessing performance of pathogenicity predictors using clinically relevant variant datasets. *J Med Genet*. 2021;58(8):547–55.

4. Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596(7873):583–9.

5. Abramson J, Adler J, Dunger J, Evans R, Green T, Pritzel A, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*. 2024;630(8016):493–500.

6. Xu M, Dantu SC, Garnett JA, Bonomo RA, Pandini A, Haider S. Functionally important residues from graph analysis of coevolved dynamic couplings. *eLife*. 2025;14:RP105005.

7. Sutto L, Marsili S, Valencia A, Gervasio FL. From residue coevolution to protein conformational ensembles and functional dynamics. *Proc Natl Acad Sci U S A*. 2015;112(44):13567–72.

8. Rodrigues CH, Pires DE, Ascher DB. DynaMut: predicting the impact of mutations on protein conformation, flexibility and stability. *Nucleic Acids Res*. 2018;46(W1):W350–5.

9. Ponzoni L, Bahar I. Structural dynamics is a determinant of the functional significance of missense variants. *Proc Natl Acad Sci U S A*. 2018;115(16):4164–9.

10. Ponzoni L, Peñaherrera DA, Oltvai ZN, Bahar I. Rhapsody: predicting the pathogenicity of human missense variants. *Bioinformatics*. 2020;36(10):3084–92.

11. Banerjee A, Bogetti AT, Bahar I. Accurate identification and mechanistic evaluation of pathogenic missense variants with Rhapsody-2. *Proc Natl Acad Sci U S A*. 2025;122(18):e2418100122.

12. Janes H, Pepe MS. Adjusting for covariates in studies of diagnostic, screening, or prognostic markers: an old concept in a new setting. *Am J Epidemiol*. 2008;168(1):89–97.

13. Ambroise C, McLachlan GJ. Selection bias in gene extraction on the basis of microarray gene-expression data. *Proc Natl Acad Sci U S A*. 2002;99(10):6562–6.

14. Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*. 2006;7:91.

15. Cawley GC, Talbot NLC. On over-fitting in model selection and subsequent selection bias in performance evaluation. *J Mach Learn Res*. 2010;11:2079–107.

16. Wainer J, Cawley G. Nested cross-validation when selecting classifiers is overzealous for most practical applications. *Expert Syst Appl*. 2021;182:115222.

17. Roberts DR, Bahn V, Ciuti S, Boyce MS, Elith J, Guillera-Arroita G, et al. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*. 2017;40(8):913–29.

18. Eid FE, Elmarakeby HA, Chan YA, Fornelos N, ElHefnawi M, Van Allen EM, et al. Systematic auditing is essential to debiasing machine learning in biology. *Commun Biol*. 2021;4:183.

19. Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. *Patterns (N Y)*. 2023;4(9):100804.

20. Whalen S, Schreiber J, Noble WS, Pollard KS. Navigating the pitfalls of applying machine learning in genomics. *Nat Rev Genet*. 2022;23(3):169–81.

21. Walsh I, Fishman D, Garcia-Gasulla D, Titma T, Pollastri G, ELIXIR Machine Learning Focus Group, et al. DOME: recommendations for supervised machine learning validation in biology. *Nat Methods*. 2021;18(10):1122–7.

22. Parker BJ, Günter S, Bedo J. Stratification bias in low signal microarray studies. *BMC Bioinformatics*. 2007;8:326.

23. Airola A, Pahikkala T, Waegeman W, De Baets B, Salakoski T. An experimental comparison of cross-validation techniques for estimating the area under the ROC curve. *Comput Stat Data Anal*. 2011;55(4):1828–44.

24. Forman G, Scholz M. Apples-to-apples in cross-validation studies: pitfalls in classifier performance measurement. *ACM SIGKDD Explor Newsl*. 2010;12(1):49–57.

25. Richards S, Aziz N, Bale S, Bick D, Das S, Gastier-Foster J, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med*. 2015;17(5):405–24.

26. Karagöl T, Karagöl A. Dynamics-aware evolutionary profiling uncouples structural rigidity from functional motion to enable enhanced variant interpretation. *Funct Integr Genomics*. 2026;26:228.

27. Landrum MJ, Lee JM, Riley GR, Jang W, Rubinstein WS, Church DM, et al. ClinVar: public archive of relationships among sequence variation and human phenotype. *Nucleic Acids Res*. 2014;42(D1):D980–5.

28. Vander Meersche Y, Cretin G, Gheeraert A, Gelly JC, Galochkina T. ATLAS: protein flexibility description from atomistic molecular dynamics simulations. *Nucleic Acids Res*. 2024;52(D1):D384–92.

29. Hopf TA, Green AG, Schubert B, Mersmann S, Schärfe CPI, Ingraham JB, et al. The EVcouplings Python framework for coevolutionary sequence analysis. *Bioinformatics*. 2019;35(9):1582–4.

30. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. *J Mol Biol*. 1982;157(1):105–32.

31. Cheng J, Novati G, Pan J, Bycroft C, Žemgulytė A, Applebaum T, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381(6664):eadg7492.

32. Tordai H, Torres O, Csepi M, Padányi R, Lukács GL, Hegedűs T. Analysis of AlphaMissense data in different protein groups and structural context. *Sci Data*. 2024;11(1):495.

33. Rastogi R, Chung R, Li S, Li C, Lee K, Woo J, et al. Critical assessment of missense variant effect predictors on disease-relevant variant data. *Hum Genet*. 2025;144(2–3):281–93.

34. Mahmood K, Jung CH, Philip G, Georgeson P, Chung J, Pope BJ, et al. Variant effect prediction tools assessed using independent, functional assay-based datasets: implications for discovery and diagnostics. *Hum Genomics*. 2017;11:10.

35. Bergquist T, Stenton SL, Nadeau EAW, Byrne AB, Greenblatt MS, Harrison SM, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med*. 2025;27(6):101402.

36. Yeo IK, Johnson RA. A new family of power transformations to improve normality or symmetry. *Biometrika*. 2000;87(4):954–9.

37. Akiba T, Sano S, Yanase T, Ohta T, Koyama M. Optuna: a next-generation hyperparameter optimization framework. In: Proceedings of the 25th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining. New York: ACM; 2019. p. 2623–31.

38. Breiman L. Random forests. *Mach Learn*. 2001;45(1):5–32.

39. Friedman JH. Greedy function approximation: a gradient boosting machine. *Ann Stat*. 2001;29(5):1189–232.

40. Breslow NE, Day NE. Statistical methods in cancer research. Volume I: the analysis of case-control studies. *IARC Sci Publ*. 1980;(32):5–338.

41. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. *Biometrics*. 1988;44(3):837–45.

42. Obuchowski NA. Nonparametric analysis of clustered ROC curve data. *Biometrics*. 1997;53(2):567–78.

43. Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS One*. 2015;10(3):e0118432.

44. Field CA, Welsh AH. Bootstrapping clustered data. *J R Stat Soc Series B Stat Methodol*. 2007;69(3):369–90.

45. Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, et al. Scikit-learn: machine learning in Python. *J Mach Learn Res*. 2011;12:2825–30.

46. Seabold S, Perktold J. Statsmodels: econometric and statistical modeling with Python. In: Proceedings of the 9th Python in Science Conference. Austin: SciPy; 2010. p. 92–6.

47. Shah N, Hou YCC, Yu HC, Sainger R, Caskey CT, Venter JC, et al. Identification of misclassified ClinVar variants via disease population prevalence. *Am J Hum Genet*. 2018;102(4):609–19.

48. Karczewski KJ, Francioli LC, Tiao G, Cummings BB, Alföldi J, Wang Q, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*. 2020;581(7809):434–43.

49. Robinson WS. Ecological correlations and the behavior of individuals. *Am Sociol Rev*. 1950;15(3):351–7.

50. Sayılgan JF, Haliloğlu T, Gönen M. Protein dynamics analysis reveals that missense mutations in cancer-related genes appear frequently on hinge-neighboring residues. *Proteins*. 2019;87(6):512–9.

51. Pejaver V, Byrne AB, Feng BJ, Pagel KA, Mooney SD, Karchin R, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109(12):2163–77.

52. Notin P, Kollasch AW, Ritter D, van Niekerk L, Paul S, Spinner H, et al. ProteinGym: large-scale benchmarks for protein fitness prediction and design. *Adv Neural Inf Process Syst*. 2023;36:64331–79.

53. Esposito D, Weile J, Shendure J, Starita LM, Papenfuss AT, Roth FP, et al. MaveDB: an open-source platform to distribute and interpret data from multiplexed assays of variant effect. *Genome Biol*. 2019;20(1):223.

54. Tsuboyama K, Dauparas J, Chen J, Laine E, Mohseni Behbahani Y, Weinstein JJ, et al. Mega-scale experimental analysis of protein folding stability in biology and design. *Nature*. 2023;620(7973):434–44.

55. Livesey BJ, Marsh JA. Variant effect predictor correlation with functional assays is reflective of clinical classification performance. *Genome Biol*. 2025;26:104.

56. Poldrack RA, Huckins G, Varoquaux G. Establishment of best practices for evidence for prediction: a review. *JAMA Psychiatry*. 2020;77(5):534–40.

57. Vihinen M. How to evaluate performance of prediction methods? Measures and their interpretation in variation effect analysis. *BMC Genomics*. 2012;13(Suppl 4):S2.

## Figures

**Figure 1. The cohort is clustered by gene, and most of the AUC compares genes.** **a** Labelled variants per gene by class; only three genes contain both classes. **b** Composition of the 1,872 pathogenic–benign pairs over which the AUC averages: 1,574 (84%) are between genes. **c** Within-gene against between-gene component of the AUC for raw scores (circles; vertical bars, stratified bootstrap 95% intervals of the within-gene component), the gene-prevalence lookup (diamonds) and the exploratory ensemble (square). Scores on the horizontal line at 0.5 carry no within-gene information.

**Figure 2. Gene memorisation under variant-level validation.** **a** ROC curves of the gene-prevalence lookup, in sample and with the scored variant excluded, and of the exploratory ensemble as reported and with its blending weights chosen without the held-out variant. **b** Factorial ablation of logistic regression: pooled AUC (filled) and its within-gene component (open) for the 16 combinations of split (V, variant; G, gene) and of scaling, selection and tuning on the full labelled set (L) or in fold (F). **c** Null distributions of the leaked variant-level AUC over 2,000 replicates with labels permuted globally or within genes; the vertical line marks the observed value (0.680).

**Figure 3. What the position-level scores encode.** **a** Published DCS against the product of conservation and normalised mobility for all ClinVar-observed residues of the cohort proteins; each protein lies on its own line, with multipliers from 10.0 to 31.7. **b** Balanced accuracy of recovering the gene from each encoding (filled, labelled cohort, eight genes; open, all ClinVar-observed residues); grey bars, 95% range of a label-permutation null; dotted line, chance for eight genes. **c** Within-gene AUC with stratified bootstrap 95% intervals and per-gene AUCs; p-values from stratified permutation.

**Figure 4. Gene-held-out evaluation and its calibration.** **a** Pooled gene-held-out AUC (filled) with its within-gene (up-triangle) and between-gene (down-triangle) components for alternative encodings of the position-level scores and of AlphaMissense residue means; † ranked or standardised among the labelled variants of each gene. **b** Null distributions of the gene-held-out AUC of the position-level scores with labels permuted globally, with each gene's scores replaced by those of a random donor protein, and with labels permuted within genes; horizontal lines mark the observed value and the AUC of a predictor that knows only the training prevalence. **c** Change in gene-held-out AUC when each feature set is added to AlphaMissense, with paired gene-cluster bootstrap 95% intervals and clustered (Obuchowski) p-values.

**Figure 5. Simulation and design.** **a** Mean difference between the variant-level cross-validated AUC and the new-gene AUC for a nearest-neighbour learner that can recognise genes, as a function of the between-gene SD of the log-odds, for four designs (within-gene effect β = 0.75; feature offset γ = 1; 150 replicates per point). **b** Mean difference between the pooled gene-held-out AUC and the new-gene AUC (solid), and between its within-gene component and the within-gene AUC, averaged over designs (dotted). **c** Genes required (five variants each) to estimate the within-gene and pooled AUC of a score with population within-gene AUC 0.65 to a 95% half-width of ±0.05 or ±0.03; the horizontal line marks the present cohort. Designs with two variants per gene and β = 0 are in Additional file 1 (Table S13, Figure S3).

**Figure 6. Cross-fold calibration and the limits of the cohort.** **a** Mean out-of-fold prediction for each held-out gene against the pathogenic fraction of the genes it was trained on, for the position-level model; colour marks whether the held-out gene contains one class or both, and the dashed line is the least-squares fit. **b** Between-gene component of the AUC under pooled leave-one-gene-out, after subtracting each fold's training prior, and under leave-two-genes-out, in which every between-gene pair is ranked by one model that saw neither of its genes; the dotted line marks 0.5. **c** Power of the two-sided stratified permutation test of AUC_{W} at level 0.05, against scores of known within-gene AUC simulated in this cohort's gene and class structure (4,000 replicates per point); vertical lines mark the within-gene AUCs observed for AlphaMissense, DCS and conservation.

## Additional files

**Additional file 1 (PDF).** Supplementary Methods S1–S10 (derivation of the AUC decomposition; reconstruction of the position-level scores; clustered paired test; random-intercept model; simulation and design calculations; exploratory ensemble; calibration and cross-fold estimation; power of the within-gene permutation test; software; contents of Additional file 2), Supplementary Tables S1–S19 (including the DOME checklist) and Supplementary Figures S1–S4.

**Additional file 2 (CSV).** Labelled cohort (88 rows): gene, UniProt accession, ATLAS identifier, ClinVar variation identifier, protein change, classification and label; the three position-level scores with their reconstructed conservation, normalised mobility, mean RMSF and per-protein multiplier; within-protein percentiles; both AlphaMissense scores; and the 16 substitution descriptors.

**Additional file 3 (ZIP).** Analysis code (Python package and scripts), requirements, the complete results directory, figure source files and the automated consistency audit.