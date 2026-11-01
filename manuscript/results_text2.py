"""Second half: gene-held-out results, simulation, design, discussion, back matter."""
import os

import pandas as pd


def write_rest2(d, V, FIG, rd):
    sim = rd("simulation.csv")
    req = rd("design_requirements.csv")
    dsim = rd("design_simulation.csv")

    def S(design, sd, beta, gamma, col):
        r = sim[(sim.design == design) & (sim.gene_sd == sd) & (sim.beta == beta) & (sim.gamma == gamma)]
        return float(r[col].iloc[0])

    def G(sd, target, est, h):
        r = req[(req.gene_sd == sd) & (req.target_within_auc == target) & (req.estimator == est)
                & (req.half_width == h)]
        return int(r.genes_required.iloc[0])

    # R6 ─────────────────────────────────────────────────────────────
    d.h(3, "Below-chance gene-held-out AUCs reflect cohort structure")
    d.p(f"Evaluated gene-held-out, a model using only DCS, RCS and normalised DCS gave a pooled "
        f"AUC of {V['obs_POS']}, a value that would conventionally be read as strong "
        f"anti-prediction. Three reference points calibrate it (Fig. 4b). A predictor that knows "
        f"nothing except the pathogenic rate of the training genes scores {V['prev_gene_loo']} "
        f"under the same design, because removing a class-pure gene shifts the class balance of "
        f"the training set against that gene. Refitting the pipeline after permuting labels "
        f"within genes, which keeps each gene's class counts and the gene-level relationship "
        f"between features and labels, gives a null centred on {V['perm_POS_within_gene_mean']} "
        f"(95% range {V['perm_POS_within_gene_lo']}–{V['perm_POS_within_gene_hi']}). Replacing "
        f"each gene's scores with residues of a randomly chosen donor protein, which keeps the "
        f"gene–class structure but breaks the link between the real feature values and the "
        f"labels, gives {V['perm_POS_donor_mean']}, close to the global null "
        f"{V['perm_POS_global_mean']}. Random protein fingerprints therefore do not depress the "
        f"pooled estimate; the depression below chance comes from the actual position-level "
        f"features' own relationship to the labels across the cohort's genes.", "first")
    d.p(f"The observed value nevertheless lies below the within-gene null (two-sided p = "
        f"{V['perm_POS_within_gene_p2']}), and the decomposition shows that both components "
        f"contribute. The within-gene component ({V['obs_POS_w']}, against a null expectation of "
        f"{V['perm_POS_within_gene_w']}) shows that the direction learnt from the other 17 genes "
        f"reverses the ordering inside PGK1, the only gene with appreciable within-gene "
        f"information (within-PGK1 AUC of DCS {V['wg_dcs_PGK1']}). The between-gene component "
        f"({V['obs_POS_b']}, against {V['perm_POS_within_gene_b']}) shows the converse: the "
        f"direction learnt with PGK1 in training runs against the ordering of the remaining "
        f"genes by prevalence. Both are expressions of one reversal between the within-gene "
        f"association in PGK1 and the across-gene association in the rest of the cohort, the "
        f"variant-level analogue of an ecological fallacy [@robinson1950], superimposed on a "
        f"structural deficit that the cohort imposes on any gene-informative predictor.")
    d.p(f"The same decomposition clarifies the alternative encodings (Fig. 4a; Additional file "
        f"1: Table S9). Before rescaling, the products gave {V['tr_unrescaled']}; conservation "
        f"alone gave {V['tr_conservation']} and normalised mobility {V['tr_mobility']}. "
        f"Within-protein percentiles gave {V['tr_pct_protein']} with a label-free reference and "
        f"{V['tr_pct_cohort']} when ranked among the labelled variants of each gene. For "
        f"AlphaMissense residue means, by contrast, per-protein min–max rescaling "
        f"({V['tr_am_res_minmax']}) and multiplication by the published DCS multipliers "
        f"({V['tr_am_res_xdcs']}) left transfer intact relative to the raw scale "
        f"({V['tr_am_res']}), whereas within-protein ranking lowered it to "
        f"{V['tr_am_res_pct']}, because it discards genuine between-protein differences in "
        f"constraint (between-gene component {V['tr_am_res_b']} raw against "
        f"{V['tr_am_res_pct_b']} ranked). Rescaling does not damage a feature that carries "
        f"between-gene information; it cannot create information that is not there. The "
        f"deficit is also not specific to gene-level splitting: under leave-one-variant-out "
        f"validation the position-level scores gave {V['obs_POSV']}, inside the range of their "
        f"within-gene null (mean {V['perm_POSV_within_gene_mean']}; two-sided p = "
        f"{V['perm_POSV_within_gene_p2']}), because in a cohort of this structure leave-one-out "
        f"predictions from weak features are anti-correlated with the held-out labels "
        f"[@parker2007,@airola2011].")

    # R6b ────────────────────────────────────────────────────────────
    d.h(3, "Cross-fold calibration accounts for part of the gene-held-out deficit")
    d.p(f"Because every between-gene pair compares predictions from two different fitted models, "
        f"AUC_{{B}} depends on how those models are calibrated relative to one another, and the "
        f"calibration of the out-of-fold predictions shows that they are not (Fig. 6a; "
        f"Additional file 1: Table S15). AlphaMissense is well calibrated under either design (recalibration slope "
        f"{V['cal_am_gene_slope']} gene-held-out, {V['cal_am_variant_slope']} variant-level, "
        f"against 1 for perfect calibration). The 19-feature model is overconfident "
        f"({V['cal_full_gene_slope']} against {V['cal_full_variant_slope']}), and the "
        f"position-level model has a negative slope ({V['cal_pos_gene_slope']}), the calibration "
        f"signature of a predictor whose ordering is reversed.", "first")
    d.p(f"The source of the miscalibration is the class balance of each training set. Across the "
        f"{V['n_genes']} gene-held-out folds, the mean prediction tracks the prevalence of the "
        f"training genes (*r* = {V['fold_r_train']}) and therefore runs against the prevalence "
        f"of the gene being predicted (*r* = {V['fold_r_gene']}), because removing a gene shifts "
        f"the training balance away from that gene's own composition "
        f"(*r* = {V['fold_r_train_gene']}). HMGCS2 is the extreme case: all nine of its variants "
        f"are pathogenic, so holding it out leaves the lowest training prevalence of any fold "
        f"({V['fold_HMGCS2_train']}) and produces the lowest mean prediction "
        f"({V['fold_HMGCS2_pred']}) for the gene that should score highest.")
    d.p(f"Two corrections test how much of the deficit this explains (Fig. 6b; Additional file 1: "
        f"Table S16). Both are applied separately "
        f"within each fold and are strictly increasing, so neither can change AUC_{{W}}. "
        f"Subtracting each fold's training prior from its predicted log-odds raises the pooled "
        f"gene-held-out AUC of the position-level scores from {V['fx_pos_none']} to "
        f"{V['fx_pos_prior']}, entirely through the between-gene component "
        f"({V['fx_pos_none_b']} to {V['fx_pos_prior_b']}); applied to a predictor that knows only "
        f"the training prevalence, the same correction returns exactly "
        f"{V['fx_prev_prior']} from {V['fx_prev_none']}, as it must. Replacing each prediction by "
        f"its percentile among the training predictions of the same fold leaves the estimate "
        f"unchanged ({V['fx_pos_rank']}), so what separates the folds is the location shift "
        f"imposed by the prior and not a difference of scale.")
    d.p(f"A design that never compares predictions across folds points the same way. Fitting a "
        f"model on all but two genes and using it to score both, repeated over all "
        f"{V['n_gene_pairs']} pairs of genes so that each of the {V['lpgo_n_pairs']} between-gene "
        f"pairs is ranked by a single model that saw neither of its genes, raises the "
        f"between-gene component of the position-level scores from {V['lpgo_pos_logo']} to "
        f"{V['lpgo_pos']} ({V['lpgo_pos_d']}) and that of AlphaMissense from "
        f"{V['lpgo_am_logo']} to {V['lpgo_am']}; for the 19-feature model it moves the other way "
        f"({V['lpgo_full_logo']} to {V['lpgo_full']}). The prior-shift correction recovers two "
        f"fifths of the distance between the pooled gene-held-out estimate and 0.5, and the "
        f"leave-two-genes-out design a quarter of the corresponding distance for the between-gene "
        f"component. What the correction leaves is a gap of {V['fx_pos_gap']} to 0.5, of which "
        f"the still-depressed between-gene component contributes {V['fx_pos_gap_b']} and the "
        f"within-gene reversal {V['fx_pos_gap_w']}; the deficit is therefore mostly, but not "
        f"only, a between-gene phenomenon. Leave-pair-out "
        f"designs remove the artefact at the cost of one fit per pair of held-out units "
        f"[@airola2011], which is affordable at this scale but not for a proteome-wide benchmark, "
        f"where the prior-shift correction is the practical option.")

    # R7 ─────────────────────────────────────────────────────────────
    d.h(3, "No feature set detectably improves on AlphaMissense")
    d.p(f"AlphaMissense alone reached {V['inc_am']} in gene-held-out logistic regression (Fig. 4c; "
        f"Additional file 1: Table S10). Adding the position-level scores as published lowered "
        f"this to {V['inc_am_pub']} ({V['inc_am_pub_d']}). The decrease is borderline: it is "
        f"significant when variants are treated as independent (DeLong p = {V['inc_am_pub_pd']}) "
        f"and its paired gene-cluster bootstrap interval excludes zero ({V['inc_am_pub_ci']}), "
        f"but Obuchowski's clustered test, which treats the {V['n_genes']} genes as the units, "
        f"does not reject (p = {V['inc_am_pub_pc']}). Conservation alone behaved the same way "
        f"({V['inc_am_cons_d']}; clustered p = {V['inc_am_cons_pc']}). No addition detectably "
        f"improved on AlphaMissense: within-protein percentiles gave {V['inc_am_pctp_d']} (clustered "
        f"p = {V['inc_am_pctp_pc']}), within-gene percentiles among the labelled variants "
        f"{V['inc_am_pctc_d']} ({V['inc_am_pctc_pc']}), substitution descriptors "
        f"{V['inc_am_sub_d']} ({V['inc_am_sub_pc']}) and the full feature set "
        f"{V['inc_am_full_d']} ({V['inc_am_full_pc']}).", "first")
    d.p(f"The percentile result depends on a detail that illustrates the general problem. "
        f"Percentiles here are mid-ranks, (rank − 1/2)/*n*, whose mean is 1/2 in every gene. "
        f"The common convention rank/*n* has mean (*n* + 1)/2*n*, which is 1 for a gene with one "
        f"variant and approaches 1/2 only for large genes, so it encodes gene size; with that "
        f"convention AlphaMissense plus within-gene percentiles reached {V['inc_am_pctc_ranknative']}, "
        f"an apparent gain of {V['inc_am_pctc_ranknative_d']} that disappears when the encoding "
        f"no longer carries gene-level information.")

    # R8 ─────────────────────────────────────────────────────────────
    d.h(3, "Simulation separates inflation, deflation and the within-gene estimand")
    s10, s5, s48, s2 = "12 genes x 10", "12 genes x 5", "48 genes x 5", "12 genes x 2"
    cm = "cohort-matched (18 genes)"
    infl = lambda dsg, sd: S(dsg, sd, 0.75, 1.0, "knn_infl_vloo_new_mean")
    gbias = lambda dsg, sd: S(dsg, sd, 0.75, 1.0, "knn_bias_gloo_new_mean")
    seen2 = S(s10, 2.0, 0.75, 1.0, "knn_bias_vloo_seen_mean")
    wcols = ["knn_bias_gloo_w_mean", "feat_bias_gloo_w_mean"]
    bw = lambda dsgs: float(sim[(sim.beta == 0.75) & (sim.design.isin(dsgs))][wcols]
                            .abs().to_numpy().max())
    allw5 = bw([s5, s10, s48])      # equal gene sizes, five or more per gene
    allwcm = bw([cm])               # observed gene-size distribution
    allw2 = bw([s2])                # two variants per gene
    lr_infl4 = S(s10, 4.0, 0.75, 1.0, "feat_infl_vloo_new_mean")
    pv = {sd: S(cm, sd, 0.0, 1.0, "prev_vloo_mean") for sd in (2.0, 3.0)}
    pg = {sd: S(cm, sd, 0.0, 1.0, "prev_gloo_mean") for sd in (2.0, 3.0)}
    d.p(f"In simulated cohorts with known estimands (Fig. 5; Additional file 1: Table S13), a "
        f"nearest-neighbour learner that can recognise genes from two position-level "
        f"descriptors was evaluated as the cohort was. Its variant-level AUC exceeded the "
        f"new-gene AUC by an amount that grew with between-gene heterogeneity and with the "
        f"number of variants per gene: with twelve genes of ten variants the excess was "
        f"{infl(s10, 0.0):+.3f} without heterogeneity, {infl(s10, 2.0):+.3f} at a between-gene SD "
        f"of 2 and {infl(s10, 4.0):+.3f} at 4, against {infl(s5, 2.0):+.3f} at SD 2 with five "
        f"variants per gene and {infl(cm, 2.0):+.3f} on the observed gene-size distribution. "
        f"Relative to the seen-gene AUC, which it is designed to estimate, the variant-level "
        f"estimate was close to unbiased ({seen2:+.3f} at SD 2 with ten variants per gene): "
        f"variant-level validation answers a legitimate question, but not the question asked of "
        f"a new gene. On uneven gene-size distributions the variant-level and seen-gene AUCs "
        f"diverge even for this learner, because leave-one-variant-out weights genes by size "
        f"whereas the seen-gene estimand weights them equally. A logistic model on a single "
        f"shifted feature, which cannot recognise genes, showed little inflation "
        f"({lr_infl4:+.3f} at SD 4).", "first")
    d.p(f"Pooled gene-level estimates were biased in the opposite direction "
        f"({gbias(s10, 2.0):+.3f} at SD 2 and {gbias(s10, 4.0):+.3f} at SD 4 with twelve genes of "
        f"ten variants; {gbias(cm, 4.0):+.3f} on the observed gene sizes at SD 4), and the bias "
        f"shrank as genes were added ({gbias(s48, 4.0):+.3f} with 48 genes at SD 4). For both "
        f"the nearest-neighbour and the logistic learner, the within-gene component of the "
        f"gene-held-out AUC tracked its estimand (largest absolute mean bias {allw5:.3f} across "
        f"the designs with five or ten variants in every gene, {allwcm:.3f} on the observed "
        f"gene-size distribution and {allw2:.3f} with two variants per gene). Without any "
        f"individual "
        f"signal, the gene-prevalence lookup on "
        f"the observed gene sizes scored {pv[2.0]:.3f}–{pv[3.0]:.3f} under variant-level and "
        f"{min(pg.values()):.3f}–{max(pg.values()):.3f} under gene-level validation at between-gene SDs of 2–3, "
        f"bracketing the values observed in the real cohort ({V['prev_variant_loo']} and "
        f"{V['prev_gene_loo']}) and consistent with the interval {V['sigma_ci']} estimated "
        f"from the random-intercept model.")

    # R9 ─────────────────────────────────────────────────────────────
    d.h(3, "Design requirements")
    g1w, g2w, g4w = G(1.0, 0.65, "within", 0.05), G(2.0, 0.65, "within", 0.05), G(4.0, 0.65, "within", 0.05)
    g2w3, g4w3 = G(2.0, 0.65, "within", 0.03), G(4.0, 0.65, "within", 0.03)
    g2p, g4p = G(2.0, 0.65, "pooled", 0.05), G(4.0, 0.65, "pooled", 0.05)
    g4w75 = G(4.0, 0.75, "within", 0.05)
    big_ = dsim[dsim.n_genes == dsim.n_genes.max()]
    mix = lambda sd: float((big_[(big_.gene_sd == sd) & (big_.target_within_auc == 0.65)].mixed_genes_mean
                            / big_.n_genes.max()).iloc[0])
    pooled_mean = lambda sd, t: float(big_[(big_.gene_sd == sd) & (big_.target_within_auc == t)].pooled_mean.iloc[0])
    d.p(f"Because only genes that contain both classes contribute within-gene pairs, the number "
        f"of such genes, not of variants, limits precision (Fig. 5c; Additional file 1: Table "
        f"S14). With five variants per gene, about {100 * mix(2.0):.0f}% of genes contain both "
        f"classes at a between-gene SD of 2 and {100 * mix(4.0):.0f}% at an SD of 4. Estimating "
        f"a within-gene AUC of 0.65 to ±0.05 then requires about {g2w:,} genes at SD 2 and "
        f"{g4w:,} at SD 4 ({g2w3:,} and {g4w3:,} for ±0.03), compared with {g1w:,} at SD 1; a "
        f"within-gene AUC of 0.75 at SD 4 requires {g4w75:,}. The pooled AUC of the same score "
        f"is estimated more precisely ({g2p:,} and {g4p:,} genes for ±0.05) but estimates a "
        f"different quantity: between-gene heterogeneity attenuates it towards 0.5, to "
        f"{pooled_mean(4.0, 0.65):.3f} for a score whose within-gene AUC is 0.65 and "
        f"{pooled_mean(4.0, 0.75):.3f} for one whose within-gene AUC is 0.75 at SD 4. The present "
        f"cohort, with {V['n_genes']} genes of which {V['n_mixed']} contain both classes, falls "
        f"far short of these requirements.", "first")

    # R9b ────────────────────────────────────────────────────────────
    d.h(3, "What this cohort can and cannot resolve")
    d.p(f"The limits can be stated exactly for the cohort as it stands. Because the within-gene "
        f"permutation null of AUC_{{W}} depends only on each gene's class counts and not on the "
        f"score, it can be computed once and used to obtain the power of the test against scores "
        f"of known within-gene discrimination in this cohort's exact structure (Fig. 6c; "
        f"Additional file 1: Table S18). The test holds its size ({V['null_size']} at AUC_{{W}} = 0.5) and its "
        f"null occupies {V['within_null_lo']}–{V['within_null_hi']}. Power reaches 80% only at a "
        f"within-gene AUC of {V['mdes_within']}; a score with a true within-gene AUC of 0.65 is "
        f"detected {V['power_065']} of the time and one with 0.70, {V['power_070']}. The "
        f"within-gene AUC observed for DCS ({V['wg_dcs']}) lies below that threshold, so its "
        f"non-significance is a statement about the cohort and not about the score.", "first")
    d.p(f"Direct comparisons between scores behave accordingly (Additional file 1: Table S18). "
        f"AlphaMissense exceeds DCS within genes by {V['pd_am_dcs']} (stratified bootstrap 95% "
        f"interval {V['pd_am_dcs_ci']}) and conservation by {V['pd_am_cons']} "
        f"({V['pd_am_cons_ci']}), both intervals excluding zero. The comparisons that matter for "
        f"the dynamics question do not: DCS exceeds conservation by {V['pd_dcs_cons']} "
        f"({V['pd_dcs_cons_ci']}), and normalised mobility falls below it by "
        f"{V['pd_nd_cons'].lstrip('+-')} ({V['pd_nd_cons_ci']}). The cohort therefore resolves "
        f"the gap between AlphaMissense and the dynamics-aware scores but not the gap between "
        f"those scores and their own conservation factor.")
    d.p(f"Every within-gene estimate here is, moreover, an estimate about one gene. Removing PGK1 "
        f"leaves {V['drop_pgk1_pairs']} of the {V['pairs_w']} within-gene pairs, and the "
        f"within-gene AUCs collapse (DCS {V['wg_dcs']} to {V['drop_pgk1_dcs']}, conservation "
        f"{V['wg_cons']} to {V['drop_pgk1_cons']}, AlphaMissense {V['wg_amv']} to "
        f"{V['drop_pgk1_am']}). The pooled gene-held-out AUC is similarly concentrated: over the "
        f"{V['n_genes']} single-gene deletions it ranges from {V['pos_drop_lo']} to "
        f"{V['pos_drop_hi']}, and the upper end is reached by deleting {V['pos_drop_worst']}, "
        f"the one gene whose {V['drop_hmgcs2_n']} variants are all pathogenic "
        f"(between-gene component {V['drop_hmgcs2_b']}; Additional file 1: Table S17). A single "
        f"class-pure gene thus accounts for most of the departure from chance.")

    # ════════════════════════════════════════════════════════ DISCUSSION
    d.h(2, "Discussion")
    d.p(f"On a gene-clustered cohort the pooled AUC is predominantly a between-gene statistic. "
        f"Here {V['b_pct']}% of the pathogenic–benign pairs compare variants from different "
        f"genes and {V['n_pure']} of {V['n_genes']} genes contribute a single class, and the "
        f"other findings follow from that fact. Variant-level validation lets a model reach a "
        f"high AUC by recognising each variant's gene and applying its prevalence, so that a "
        f"lookup using no variant information matches a tuned five-model ensemble. Gene-level "
        f"validation removes that route but brings its own artefact: between-gene pairs then "
        f"compare predictions from different fitted models, and removing a class-pure gene shifts "
        f"the training class balance against it, so the pooled estimate is pulled below its "
        f"target, to {V['prev_gene_loo']} for a predictor that knows only the training "
        f"prevalence. The within-gene component isolates residue-level discrimination from both of "
        f"these between-gene artefacts, because a score that is constant within genes has "
        f"AUC_{{W}} = 0.5 whatever the design. It is not, however, invariant to the resampling "
        f"unit once a model is fitted: under leave-one-variant-out the model is shown the "
        f"held-out variant's gene-mates and their labels, which raised AUC_{{W}} by "
        f"{V['eff_split_w'].lstrip('+')} here, so the within-gene component should be read from "
        f"gene-held-out predictions.", "first")
    d.p(f"For the dynamics-aware scores, the within-gene evidence is modest (AUC {V['wg_dcs']}, "
        f"p = {V['wg_dcs_p']}, almost entirely from PGK1) and is the evidence carried by the "
        f"conservation factor; normalised mobility shows none. This agrees with the original "
        f"description of the scores, in which RCS tracked conservation against AlphaMissense "
        f"classes while DCS did not [@karagol2026]. It does not show that dynamics is irrelevant "
        f"to variant interpretation: Rhapsody-2's protein-held-out performance and the enrichment "
        f"of disease mutations at hinge-adjacent residues are evidence to the contrary "
        f"[@banerjee2025,@ponzoni2018,@sayilgan2019], and three genes with both classes cannot "
        f"test a mechanism. What it shows is that this cohort cannot detect a dynamics "
        f"contribution, and that a validation design resting on variant-level resampling measures "
        f"gene identity. The result also underlines how strong the sequence-based baseline has "
        f"become: AlphaMissense discriminated within genes as well as between them, and no "
        f"addition detectably improved its gene-held-out performance (the largest change, "
        f"{V['inc_am_pub_d']} from the position-level scores, was significant by the "
        f"variant-level DeLong test and its bootstrap interval but not by the clustered test).")
    d.p(f"Per-protein normalisation is a modelling decision, but it was not the decision that "
        f"caused the problems here. It multiplies each protein's raw values by a factor set by "
        f"that protein's most extreme residue ({V['mult_min']} to {V['mult_max']} for DCS in this "
        f"cohort) and so distorts cross-protein comparisons, yet the products before rescaling "
        f"identify proteins almost as well, and applying the published multipliers to "
        f"AlphaMissense added protein identity without degrading its transfer to held-out genes. "
        f"Protein identity is carried by the joint distribution "
        f"of conservation and mobility, which differs between proteins whatever the scaling. Only "
        f"an encoding that equalises the whole marginal distribution, the within-protein rank, "
        f"removes it, and removal has a cost: ranking AlphaMissense residue means within "
        f"proteins lowered their gene-held-out AUC from {V['tr_am_res']} to {V['tr_am_res_pct']} "
        f"because between-protein differences in constraint are real signal. The reference frame "
        f"of a feature should therefore follow the estimand: within-protein ranks when "
        f"within-gene discrimination is the target, proteome-calibrated scales when between-gene "
        f"information is wanted.")
    d.p("The problem differs from the circularity described by Grimm and colleagues "
        "[@grimm2015], which concerns overlap between the data used to train a published "
        "predictor and the data used to evaluate it. Here the training set is the evaluation "
        "cohort itself; their type 1 circularity, variants of one gene on both sides of the "
        "split, is the closest analogue, with the addition that the gene relationship is legible "
        "from the features as well as from the labels. In the taxonomy of Kapoor and Narayanan "
        "[@kapoor2023] it lies between the absence of an independent test set and the use of "
        "illegitimate features. The complementary bias of pooled gene-held-out estimates is the "
        "cluster-level counterpart of the stratification bias described for leave-one-out "
        "designs [@parker2007,@airola2011,@forman2010]. Gene-level prevalence is not "
        "illegitimate information in itself, since PP2 and BP1 encode it [@richards2015], but a "
        "predictor that has learnt it should not be combined with those criteria without "
        "accounting for the overlap, and its benchmark performance says little about "
        "residue-level discrimination.")
    d.p("These observations translate into a short protocol (Table 3). The composition of "
        "pathogenic–benign pairs and the within/between decomposition should accompany any "
        "pooled AUC on a clustered cohort. A gene-prevalence baseline, evaluated under the same "
        "design, shows how much of a reported value gene identity alone can produce, and a "
        "within-gene permutation null tests whether anything remains. Gene-level resampling, "
        "with feature selection and tuning nested and grouped by gene, should be the default "
        "whenever features describe positions or proteins, with the within-gene component "
        "reported beside the pooled AUC, because the pooled gene-held-out value has its own "
        "downward bias. That bias is partly removable: subtracting each fold's training prior "
        "from its predicted log-odds costs nothing, cannot affect the within-gene component, and "
        "recovered two fifths of the deficit here, while leave-pair-out designs remove the "
        "cross-fold comparison entirely where the number of held-out units makes them affordable. "
        "A null result should be accompanied by the effect it could have detected, which for a "
        "within-gene comparison follows from the class counts alone. Intervals and paired tests "
        "should treat genes as the independent units "
        f"[@field2007,@obuchowski1997]; here the clustered and variant-level tests disagreed on "
        f"{V['n_test_disagree']} of the {V['n_tests']} paired comparisons examined, and the "
        f"clustered test was the more significant of the two in "
        f"{V['n_clustered_smaller']} of them, so it is not uniformly conservative (Additional "
        f"file 1: Table S2). For position-level or "
        "protein-derived features, a "
        "gene-identity recovery test and an explicit statement of the reference frame belong "
        "in the methods. Benchmarks should be planned by counting genes that contain both "
        "classes, the quantity that limits within-gene precision. Each of these steps is "
        "inexpensive, requires no additional data, and together they separate residue-level "
        "signal from gene-level structure before a pooled AUC is interpreted.")
    d.p("The limitations are those of the cohort. Only three genes contain both classes, so "
        "within-gene estimates are dominated by PGK1 and every interval is wide. ClinVar "
        "classifications draw on computational evidence that includes conservation "
        "[@grimm2015,@pejaver2022], so conservation-derived features are evaluated "
        "optimistically. The reconstruction depends on residues of uncertain "
        "significance in the same proteins, which were available for every labelled protein "
        "here but need not be in general. The audit of the exploratory ensemble could "
        "re-estimate only its blending step, because some of its input features are not "
        "reconstructible from the published scores. The simulations establish generality only "
        "within the generative model specified. Deep mutational "
        "scanning resources [@notin2023,@esposito2019,@tsuboyama2023] evaluate predictors "
        "within proteins and so avoid between-gene confounding by construction, although assay "
        "scores and clinical classifications correspond imperfectly [@livesey2025]; a "
        "replication there, and in ClinVar-scale cohorts with gene-held-out designs, is the "
        "natural next step. The present analysis is a precondition for such work rather than a "
        "substitute for it: performance reported on clustered cohorts is only as informative as "
        "the baseline it is compared with [@poldrack2020,@vihinen2012].")

    # ════════════════════════════════════════════════════════ CONCLUSIONS
    d.h(2, "Conclusions")
    d.p(f"On a gene-clustered ClinVar cohort, the AUC of a structure-informed classifier was "
        f"dominated by comparisons between genes. Variant-level cross-validation rewarded gene "
        f"memorisation, so that a prevalence lookup ({V['prev_variant_loo']}) matched a tuned "
        f"ensemble ({V['ens']}), and pooled gene-held-out validation penalised gene-level "
        f"information, producing below-chance values ({V['obs_POS']}) that are largely explained "
        f"by the cohort's gene structure. Decomposing the AUC into within-gene and between-gene "
        f"components, together with a gene-prevalence baseline, a within-gene permutation null "
        f"and clustered inference, separates these effects. Applied here, they show that the "
        f"dynamics-aware scores carry no within-gene information beyond that of their conservation "
        f"factor and no detectable improvement on AlphaMissense, and that resolving an "
        f"incremental benefit in either direction requires cohorts with many genes that contain "
        f"both classes.", "first")

    # ════════════════════════════════════════════════════════ BACK MATTER
    d.h(2, "Abbreviations")
    d.p("ACMG/AMP: American College of Medical Genetics and Genomics/Association for Molecular "
        "Pathology; AP: average precision; AUC: area under the receiver operating characteristic "
        "curve; B/LB: benign or likely benign; DCS: Dynamic Conserved Score; LOO: leave-one-out; "
        "MD: molecular dynamics; P/LP: pathogenic or likely pathogenic; RCS: Rigid Conserved "
        "Score; RMSF: root-mean-square fluctuation; ROC: receiver operating characteristic.",
        "first")
    d.h(2, "Declarations")
    d.p("**Ethics approval and consent to participate.** Not applicable. This computational study "
        "used publicly available data and did not involve human participants, animal subjects or "
        "identifiable data.", "first")
    d.p("**Consent for publication.** Not applicable.", "first")
    d.p("**Availability of data and materials.** ClinVar classifications are available at "
        "https://www.ncbi.nlm.nih.gov/clinvar/. Molecular dynamics trajectories are available "
        "from ATLAS (https://www.dsimb.inserm.fr/ATLAS). AlphaMissense scores are available at "
        "https://alphamissense.hegelab.org. The dynamics-aware scores were computed as described "
        "previously [@karagol2026] and can be regenerated with ADEPT "
        "(https://www.karagolresearch.com/adept). The labelled cohort with ClinVar variation "
        "identifiers, reconstructed components and all features is provided as Additional file "
        "2. All analysis code, intermediate results, out-of-fold predictions, permutation "
        "distributions, simulation outputs and the consistency audit are provided as Additional "
        "file 3 and are deposited at "
        "https://github.com/karagol-taner/Dynamics-aware-Evolutionary-Profiling.", "first")
    d.p("**Competing interests.** The authors declare no competing interests.", "first")
    d.p("**Funding.** The authors received no specific funding for this work.", "first")
    d.p("**Authors' contributions.** T.K. and A.K. contributed equally. Conceptualisation, T.K., "
        "A.K.; methodology, T.K., A.K.; software, T.K., A.K.; formal analysis, T.K., A.K.; "
        "investigation, T.K., A.K.; visualisation, T.K., A.K.; writing, original draft, T.K., "
        "A.K.; writing, review and editing, T.K., A.K. Both authors read and approved the final "
        "manuscript.", "first")
    d.p("**Acknowledgements.** Not applicable.", "first")

    # tables before references so that any citation inside them is numbered in order
    tables(d, V, rd)
    d.h(2, "References")
    d.refs()
    legends(d, V, FIG)
    d.h(2, "Additional files")
    d.p("**Additional file 1 (PDF).** Supplementary Methods S1–S10 (derivation of the AUC "
        "decomposition; reconstruction of the position-level scores; clustered paired test; "
        "random-intercept model; simulation and design calculations; exploratory ensemble; "
        "calibration and cross-fold estimation; power of the within-gene permutation test; "
        "software; contents of Additional file 2), Supplementary Tables S1–S19 (including the "
        "DOME checklist) and Supplementary Figures S1–S4.", "first")
    d.p("**Additional file 2 (CSV).** Labelled cohort (88 rows): gene, UniProt accession, ATLAS "
        "identifier, ClinVar variation identifier, protein change, classification and label; the "
        "three position-level scores with their reconstructed conservation, normalised mobility, "
        "mean RMSF and per-protein multiplier; within-protein percentiles; both AlphaMissense "
        "scores; and the 16 substitution descriptors.", "first")
    d.p("**Additional file 3 (ZIP).** Analysis code (Python package and scripts), requirements, "
        "the complete results directory, figure source files and the automated consistency "
        "audit.", "first")


def tables(d, V, rd):
    d.pagebreak()
    d.h(2, "Tables")
    rows = [
        ["Gene prevalence, in sample", "—", V["prev_in_sample"], V["prev_in_sample_w"], V["prev_in_sample_b"], "Uses no variant-level information"],
        ["Gene prevalence, variant excluded", "Variant", V["prev_variant_loo"], V["prev_variant_loo_w"], V["prev_variant_loo_b"], f"Cluster 95% CI {V['prev_loo_ci']}"],
        ["Gene prevalence, gene excluded", "Gene", V["prev_gene_loo"], V["prev_gene_loo_w"], V["prev_gene_loo_b"], "Training prevalence only; prior-shift reference"],
        ["Exploratory ensemble, as reported", "Variant", V["ens"], V["ens_within"], V["ens_between"], f"Cluster 95% CI {V['ens_ci']}"],
        ["Exploratory ensemble, nested blend", "Variant", V["ens_nested"], V["ens_nested_w"], V["ens_nested_b"], "Weights chosen without the held-out variant"],
        ["C1 pipeline, labels permuted globally", "Variant", V["perm_C1_global_mean"], V["perm_C1_global_w"], V["perm_C1_global_b"], f"Observed {V['obs_C1']}: p = {V['perm_C1_global_pu']}"],
        ["C1 pipeline, labels permuted within gene", "Variant", V["perm_C1_within_gene_mean"], V["perm_C1_within_gene_w"], V["perm_C1_within_gene_b"], f"Observed {V['obs_C1']}: p = {V['perm_C1_within_gene_pu']}"],
        ["C8 pipeline, labels permuted within gene", "Gene", V["perm_C8_within_gene_mean"], V["perm_C8_within_gene_w"], V["perm_C8_within_gene_b"], f"Observed {V['obs_C8']}: p = {V['perm_C8_within_gene_pu']}"],
        ["Position-level scores", "Gene", V["obs_POS"], V["obs_POS_w"], V["obs_POS_b"], "DCS, RCS, normalised DCS"],
        ["  same, labels permuted within gene", "Gene", V["perm_POS_within_gene_mean"], V["perm_POS_within_gene_w"], V["perm_POS_within_gene_b"], f"Observed below null: two-sided p = {V['perm_POS_within_gene_p2']}"],
        ["  same, donor-protein scores", "Gene", V["perm_POS_donor_mean"], V["perm_POS_donor_w"], V["perm_POS_donor_b"], "Protein fingerprint without signal"],
        ["Substitution descriptors", "Gene", V["tr_substitution"], V["tr_substitution_w"], V["tr_substitution_b"], "16 descriptors, nine selected in fold"],
        ["AlphaMissense, raw score", "—", V["wg_amv_auc"], V["wg_amv"], V["wg_amv_b"], f"Cluster 95% CI {V['wg_amv_aucci']}; no model fitted"],
        ["AlphaMissense, logistic model", "Gene", V["inc_am"], V["inc_am_w"], V["inc_am_b"], "Matched comparator for incremental value"],
        ["AlphaMissense + position-level scores", "Gene", V["inc_am_pub"], V["inc_am_pub_w"], V["inc_am_pub_b"], f"DeLong p = {V['inc_am_pub_pd']}; clustered p = {V['inc_am_pub_pc']}"],
    ]
    d.table("**Table 1. Baselines, null distributions and comparators.** Values are AUCs with "
            "their within-gene (AUC_{W}) and between-gene (AUC_{B}) components; null rows give "
            "the mean over replicates.",
            [["Baseline or control", "Resampling", "AUC", "AUC_{W}", "AUC_{B}", "Note"]], rows,
            widths=[2.6, 0.9, 0.7, 0.7, 0.7, 2.6], align=["l", "c", "c", "c", "c", "l"],
            font_size=8.5)

    lr = rd("ablation_lr16.csv")
    rows = []
    for _, r in lr.iterrows():
        rows.append([r.config, r.split.capitalize(), r.scaling, r.selection, r.tuning,
                     f"{r.auc:.3f}", f"{r.ap:.3f}", f"{r.auc_within:.3f}", f"{r.auc_between:.3f}",
                     f"{r.median_C:.3g}"])
    d.table("**Table 2. Logistic regression under every combination of four validation "
            "choices.** All configurations use the same 19 features (three position-level scores "
            f"and 16 substitution descriptors) on the same {V['n']} variants in {V['n_genes']} "
            "genes. *Split*: leave-one-variant-out (Variant) or leave-one-gene-out (Gene). "
            "*Scaling*, *selection* and *tuning*: median imputation with Yeo–Johnson "
            "transformation, ANOVA F selection of nine features, and choice of the L2 penalty C, "
            "each fitted on the full labelled set or inside every training fold. AUC_{W} and "
            "AUC_{B}: within-gene and between-gene components of the AUC.",
            [["Config.", "Split", "Scaling", "Selection", "Tuning", "AUC", "AP", "AUC_{W}",
              "AUC_{B}", "Median C"]], rows,
            widths=[0.8, 0.9, 1.0, 1.0, 1.0, 0.75, 0.75, 0.75, 0.75, 0.9],
            notes=f"L16 is the only configuration in which no stage sees the held-out gene. The "
                  f"2^{{3}} design without tuning for logistic regression (C = 1), random forest "
                  f"and gradient boosting is given in Additional file 1: Table S5; factorial "
                  f"effects in Table S6. Main effects (variant minus gene, full set minus in fold): "
                  f"split {V['eff_split']}, scaling {V['eff_scaling']}, selection "
                  f"{V['eff_selection']}, tuning {V['eff_tuning']}.", font_size=8.5)

    rows = [
        ["Pair composition", "What fraction of the AUC compares variants of different genes?",
         "Count pathogenic–benign pairs within and between genes", f"{V['b_pct']}% between genes; {V['n_mixed']} of {V['n_genes']} genes mixed"],
        ["AUC decomposition", "Does a score discriminate within genes, between them, or both?",
         "AUC = w·AUC_{W} + (1 − w)·AUC_{B}", f"Ensemble {V['ens_within']} within, {V['ens_between']} between"],
        ["Gene-prevalence baseline", "How much can gene identity alone produce?",
         "Score each variant by its gene's pathogenic fraction under the evaluated design", f"{V['prev_variant_loo']} (variant-level); {V['prev_gene_loo']} (gene-level)"],
        ["Within-gene permutation null", "Is there signal beyond gene-level structure?",
         "Permute labels within genes; refit the full pipeline", f"p = {V['perm_C1_within_gene_pu']} (variant-level pipeline)"],
        ["Nested, gene-grouped selection and tuning", "Do preprocessing, selection, tuning or blending see the held-out data?",
         "Refit every data-dependent step inside each training fold", f"Blending: {V['ens']} → {V['ens_nested']}"],
        ["Clustered inference", "Are intervals and paired tests anticonservative?",
         "Gene-cluster bootstrap; Obuchowski test on gene-level df", f"AM vs AM + scores: DeLong p = {V['inc_am_pub_pd']}, clustered p = {V['inc_am_pub_pc']}"],
        ["Fold calibration", "Are the fold-specific models on a common scale?",
         "Calibration slope per design; subtract each fold's training prior", f"Slope {V['cal_pos_gene_slope']}; AUC_{{B}} {V['fx_pos_none_b']} → {V['fx_pos_prior_b']}"],
        ["Detectable effect", "Is a null result informative?",
         "Power of the within-gene test against known AUC_{W}, at the observed class counts", f"80% power only at AUC_{{W}} = {V['mdes_within']}"],
        ["Gene-identity recovery", "Do the features encode the protein?",
         "Predict gene from features, residue-grouped CV, balanced accuracy", f"{V['gid_labelled_published']} (chance 0.125); ranks {V['gid_labelled_pct_protein']}"],
        ["Mixed-gene count", "Can the cohort estimate within-gene discrimination?",
         "Count genes with both classes; plan by simulation", f"{V['n_mixed']} of {V['n_genes']} contain both classes; hundreds of mixed genes needed"],
    ]
    d.table("**Table 3. Diagnostic protocol for benchmarks on gene-clustered cohorts.**",
            [["Diagnostic", "Question", "Computation", "This cohort"]], rows,
            widths=[1.5, 2.2, 2.4, 2.0], align=["l", "l", "l", "l"], font_size=8.5)


def legends(d, V, FIG):
    d.pagebreak()
    d.h(2, "Figures")
    d.figure(f"{FIG}/Figure1.png",
             f"**Figure 1. The cohort is clustered by gene, and most of the AUC compares genes.** "
             f"**a** Labelled variants per gene by class; only three genes contain both classes. "
             f"**b** Composition of the {V['pairs']} pathogenic–benign pairs over which the AUC "
             f"averages: {V['pairs_b']} ({V['b_pct']}%) are between genes. **c** Within-gene "
             f"against between-gene component of the AUC for raw scores (circles; vertical bars, "
             f"stratified bootstrap 95% intervals of the within-gene component), the "
             f"gene-prevalence lookup (diamonds) and the exploratory ensemble (square). Scores on "
             f"the horizontal line at 0.5 carry no within-gene information.", alt="Figure 1")
    d.figure(f"{FIG}/Figure2.png",
             f"**Figure 2. Gene memorisation under variant-level validation.** **a** ROC curves of "
             f"the gene-prevalence lookup, in sample and with the scored variant excluded, and of "
             f"the exploratory ensemble as reported and with its blending weights chosen without "
             f"the held-out variant. **b** Factorial ablation of logistic regression: pooled AUC "
             f"(filled) and its within-gene component (open) for the 16 combinations of split "
             f"(V, variant; G, gene) and of scaling, selection and tuning on the full labelled set "
             f"(L) or in fold (F). **c** Null distributions of the leaked variant-level AUC over "
             f"{V['perm_C1_global_n']} replicates with labels permuted globally or within genes; "
             f"the vertical line marks the observed value ({V['obs_C1']}).", alt="Figure 2")
    d.figure(f"{FIG}/Figure3.png",
             f"**Figure 3. What the position-level scores encode.** **a** Published DCS against "
             f"the product of conservation and normalised mobility for all ClinVar-observed "
             f"residues of the cohort proteins; each protein lies on its own line, with "
             f"multipliers from {V['mult_min']} to {V['mult_max']}. **b** Balanced accuracy of "
             f"recovering the gene from each encoding (filled, labelled cohort, eight genes; open, "
             f"all ClinVar-observed residues); grey bars, 95% range of a label-permutation null; "
             f"dotted line, chance for eight genes. **c** Within-gene AUC with stratified "
             f"bootstrap 95% intervals and per-gene AUCs; p-values from stratified "
             f"permutation.", alt="Figure 3")
    d.figure(f"{FIG}/Figure4.png",
             f"**Figure 4. Gene-held-out evaluation and its calibration.** **a** Pooled "
             f"gene-held-out AUC (filled) with its within-gene (up-triangle) and between-gene "
             f"(down-triangle) components for alternative encodings of the position-level scores "
             f"and of AlphaMissense residue means; † ranked or standardised among the labelled "
             f"variants of each gene. **b** Null distributions of the gene-held-out AUC of the "
             f"position-level scores with labels permuted globally, with each gene's scores "
             f"replaced by those of a random donor protein, and with labels permuted within genes; "
             f"horizontal lines mark the observed value and the AUC of a predictor that knows only "
             f"the training prevalence. **c** Change in gene-held-out AUC when each feature set is "
             f"added to AlphaMissense, with paired gene-cluster bootstrap 95% intervals and "
             f"clustered (Obuchowski) p-values.", alt="Figure 4")
    d.figure(f"{FIG}/Figure5.png",
             "**Figure 5. Simulation and design.** **a** Mean difference between the variant-level "
             "cross-validated AUC and the new-gene AUC for a nearest-neighbour learner that can "
             "recognise genes, as a function of the between-gene SD of the log-odds, for four "
             "designs (within-gene effect β = 0.75; feature offset γ = 1; 150 replicates per point). "
             "**b** Mean difference between the pooled gene-held-out AUC and the new-gene AUC "
             "(solid), and between its within-gene component and the within-gene AUC, averaged over "
             "designs (dotted). **c** Genes required (five variants each) to estimate the within-gene "
             "and pooled AUC of a score with population within-gene AUC 0.65 to a 95% half-width of "
             "±0.05 or ±0.03; the horizontal line marks the present cohort. Designs with two variants "
             "per gene and β = 0 are in Additional file 1 (Table S13, Figure S3).", alt="Figure 5")
    d.figure(f"{FIG}/Figure6.png",
             f"**Figure 6. Cross-fold calibration and the limits of the cohort.** **a** Mean "
             f"out-of-fold prediction for each held-out gene against the pathogenic fraction of "
             f"the genes it was trained on, for the position-level model; colour marks whether "
             f"the held-out gene contains one class or both, and the dashed line is the "
             f"least-squares fit. **b** Between-gene component of the AUC under pooled "
             f"leave-one-gene-out, after subtracting each fold's training prior, and under "
             f"leave-two-genes-out, in which every between-gene pair is ranked by one model that "
             f"saw neither of its genes; the dotted line marks 0.5. **c** Power of the two-sided "
             f"stratified permutation test of AUC_{{W}} at level 0.05, against scores of known "
             f"within-gene AUC simulated in this cohort's gene and class structure (4,000 "
             f"replicates per point); vertical lines mark the within-gene AUCs observed for "
             f"AlphaMissense, DCS and conservation.", alt="Figure 6")
