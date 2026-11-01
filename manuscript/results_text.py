"""Results, Discussion, Conclusions, Declarations, tables and figure legends."""
import os

import pandas as pd


def write_rest(d, V, FIG):
    R = os.path.join(os.path.dirname(FIG))
    rd = lambda n: pd.read_csv(os.path.join(R, n))

    # ════════════════════════════════════════════════════════ RESULTS
    d.h(2, "Results")
    d.h(3, "Most pathogenic–benign comparisons in the cohort are between genes")
    d.p(f"The {V['n']} labelled variants are distributed very unevenly (Fig. 1a). PGK1 contributes "
        f"{V['n_pgk1']} ({V['pct_pgk1']}%), GBP1 contributes 13 benign variants and no pathogenic "
        f"ones, and HMGCS2 contributes nine pathogenic variants and no benign ones. Of the "
        f"{V['n_genes']} genes, {V['n_pure']} (holding {V['n_var_pure']} variants) contain a "
        f"single class; only PGK1, GHRHR and KAT8 contain both. Gene and class are strongly "
        f"associated (χ² = {V['chi2']}, df = {V['chi2_df']}, Monte Carlo p {V['chi2_p_mc']}; "
        f"Cramér's V = {V['cramer']}), and a random-intercept model places the between-gene "
        f"standard deviation of the log-odds of pathogenicity at {V['sigma_g']} (profile-likelihood "
        f"95% interval {V['sigma_ci']}; latent intraclass correlation {V['icc']}, {V['icc_ci']}).",
        "first")
    d.p(f"Of the {V['pairs']} pathogenic–benign pairs over "
        f"which the AUC averages, {V['pairs_w']} ({V['w_pct']}%) lie within a gene "
        f"({V['pairs_pgk1']} in PGK1, {V['pairs_ghrhr']} in GHRHR and {V['pairs_kat8']} in KAT8) "
        f"and {V['pairs_b']} ({V['b_pct']}%) lie between genes (Fig. 1b). An AUC computed on this "
        f"cohort is therefore predominantly a statement about how a score orders genes. "
        f"The two components separate scores that pooled AUCs conflate (Fig. 1c). AlphaMissense "
        f"discriminates within genes as well as between them (within-gene AUC {V['wg_amv']} "
        f"against {V['wg_amv_b']} between genes); the position-level scores discriminate only "
        f"within genes (DCS {V['wg_dcs']} within, {V['wg_dcs_b']} between); and the "
        f"gene-prevalence lookup discriminates only between genes ({V['prev_in_sample_b']}), "
        f"which is enough for an in-sample AUC of {V['prev_in_sample']}.")
    d.p("This concentration is not peculiar to the present cohort. It follows from requiring both a "
        "ClinVar classification and an existing molecular dynamics trajectory, and similar "
        "concentration is a general property of ClinVar-derived benchmarks, in which a few "
        "intensively studied genes dominate and classification quality varies with how heavily a "
        "gene has been studied [@grimm2015,@shah2018,@karczewski2020].")

    d.h(3, "A gene-prevalence lookup matches a tuned ensemble under variant-level validation")
    d.p(f"A classifier that ignores variants and scores each by the pathogenic fraction of its "
        f"gene reaches AUC {V['prev_in_sample']} in sample and {V['prev_variant_loo']} when the "
        f"scored variant is excluded from its own gene's fraction (Fig. 2a; Table 1). Excluding the "
        f"variant reverses the ordering inside every mixed gene, so the within-gene component of "
        f"this baseline is exactly {V['prev_variant_loo_w']}; all of its performance is "
        f"between-gene ({V['prev_variant_loo_b']}).", "first")
    w = V["ens_w"]
    d.p(f"The exploratory ensemble reached {V['ens']} under the same design. Its weights, "
        f"recovered exactly from the distributed predictions, were {100 * w['Logistic']:.1f}% "
        f"logistic regression, {100 * w['LightGBM']:.1f}% LightGBM, {100 * w['XGBoost']:.1f}% "
        f"XGBoost, {100 * w['RandomForest']:.1f}% random forest and {100 * w['CatBoost']:.1f}% "
        f"CatBoost; the five base learners alone gave {V['ens_base_rng']} and an equal-weight "
        f"blend {V['ens_eq']}. When the weights were chosen without the held-out variant, the "
        f"blend gave {V['ens_nested']}, and a nested choice of the single best learner gave "
        f"{V['ens_nested_single']}. Choosing the weights without the held-out variant costs "
        f"{V['ens_nested_gap']} against the best in-sample choice from the same candidate set "
        f"({V['ens_insample_best']}); the reported {V['ens']} exceeds that best candidate, so "
        f"the distributed weights are not reproducible from this set. The ensemble "
        f"exceeded the lookup by only {V['ens_vs_prev_delta']}, a difference the data cannot "
        f"distinguish from zero (clustered standard error {V['ens_vs_prev_se']}; clustered "
        f"p = {V['ens_vs_prev_p_clu']}; DeLong p = {V['ens_vs_prev_p_delong']}). Unlike the "
        f"lookup, whose within-gene component is {V['prev_variant_loo_w']} because excluding a "
        f"variant from its own gene's fraction inverts the ordering inside each gene, the "
        f"ensemble retains a within-gene component of {V['ens_within']}. That advantage applies "
        f"to the {V['w_pct']}% of pairs that lie within a gene; over the remaining "
        f"{V['b_pct']}% the lookup is the better of the two ({V['prev_variant_loo_b']} against "
        f"{V['ens_between']}), and the pooled AUC is the pair-weighted average of the two.")

    d.h(3, "The resampling unit acts on the between-gene component")
    d.p(f"Four pipeline stages were crossed with the resampling unit. At a fixed penalty (C = 1), "
        f"across three model families (Additional file 1: Table S5), the fully leaked "
        f"configuration — variant-level split with preprocessing and feature selection fitted on "
        f"the complete labelled set — gave {V['C1_lr']} for logistic regression, and the fully "
        f"controlled configuration, with a gene-level split and nothing fitted outside the "
        f"training fold, gave {V['C8_lr']}, a difference of {V['C1_C8_lr']}. Random forests "
        f"({V['C1_rf']} to {V['C8_rf']}) and gradient boosting ({V['C1_gb']} to {V['C8_gb']}) "
        f"moved in the same direction, and for every model family every variant-level "
        f"configuration exceeded every gene-level one.", "first")
    d.p(f"Adding the placement of penalty tuning as a fourth factor gives the 16 "
        f"logistic-regression configurations of Table 2 and Fig. 2b. The size of the contrast "
        f"depends on how the penalty is tuned. In the only configuration in which no stage sees "
        f"the held-out gene (L16, with C tuned by a gene-grouped inner cross-validation inside "
        f"each fold), the AUC was {V['L16']} (median C = {V['L16_C']}), whereas tuning on the "
        f"complete set by variant-level folds, as in the exploratory pipeline, chose "
        f"C = {V['L15_C']} and gave {V['L15']}. Averaged over the other factors, variant-level "
        f"splitting raised the AUC by {V['eff_split'].lstrip('+')}, leaked feature selection by "
        f"{V['eff_selection'].lstrip('+')}, and leaked tuning changed it by {V['eff_tuning']}; the "
        f"largest interaction was between the split and leaked selection "
        f"({V['eff_split_x_selection']}), because selection on the complete set matters most "
        f"when the split is otherwise correct (Additional file 1: Table S6).")
    d.p(f"The decomposition shows where these differences arise. Across the eight gene-level "
        f"logistic configurations the within-gene component stayed within {V['lr16_gene_w_rng']}, "
        f"whereas the between-gene component spanned {V['lr16_gene_b_rng']} and carried almost all "
        f"of the variation in pooled AUC ({V['lr16_gene_rng']}). Changing the resampling unit acts "
        f"chiefly on the between-gene component — where gene memorisation operates under "
        f"variant-level splitting and cross-fold calibration under gene-level splitting — but it "
        f"does not leave the within-gene component untouched: for logistic regression the "
        f"variant-level configurations raise it ({V['lr16_var_w_rng']} against "
        f"{V['lr16_gene_w_rng']} across the gene-level configurations). The resampling unit is "
        f"the cause, not the leaked preprocessing: averaged over the other factors the split "
        f"moves AUC_{{W}} by {V['eff_split_w']}, whereas scaling, selection and tuning on the "
        f"complete set move it by {V['eff_scaling_w']}, {V['eff_selection_w']} and "
        f"{V['eff_tuning_w']}. Under leave-one-variant-out the model is shown the held-out "
        f"variant's gene-mates and their labels, which is itself gene-level information acting "
        f"on the within-gene ordering. For random forests the within-gene component is instead "
        f"slightly higher under gene-level splitting ({V['rf_w_gene_rng']} against "
        f"{V['rf_w_var_rng']}).")

    d.h(3, "Permutation nulls that preserve gene structure")
    d.p(f"Permuting labels globally destroys both individual and gene-level signal; permuting "
        f"them within genes destroys individual signal while holding every gene's class counts "
        f"fixed (Fig. 2c; Table 1; Additional file 1: Table S4). Over {V['perm_C1_global_n']} replicates "
        f"of the leaked variant-level pipeline, the global null had mean "
        f"{V['perm_C1_global_mean']} (95% range {V['perm_C1_global_lo']}–"
        f"{V['perm_C1_global_hi']}), against which the observed {V['obs_C1']} gives "
        f"p = {V['perm_C1_global_pu']} (Monte Carlo SE {V['perm_C1_global_mcse']}). The "
        f"within-gene null had mean {V['perm_C1_within_gene_mean']} "
        f"({V['perm_C1_within_gene_lo']}–{V['perm_C1_within_gene_hi']}), against which the "
        f"observed value is not significant (p = {V['perm_C1_within_gene_pu']}; SE "
        f"{V['perm_C1_within_gene_mcse']}). A pipeline given labels that carry no "
        f"individual-variant information at all thus returns {V['perm_C1_within_gene_mean']} "
        f"under variant-level resampling, and the observed {V['obs_C1']} does not depart "
        f"significantly from it.", "first")
    d.p(f"The same nulls behave differently under gene-level resampling. For the controlled "
        f"configuration C8, the within-gene null had mean {V['perm_C8_within_gene_mean']}, below "
        f"chance, and the global null {V['perm_C8_global_mean']}; the observed {V['obs_C8']} "
        f"gives p = {V['perm_C8_within_gene_pu']} and p = {V['perm_C8_global_pu']}, "
        f"respectively. Neither design provides significant evidence of individual-variant "
        f"discrimination by the 19-feature pipeline, and under gene-level resampling the "
        f"reference point for 'no signal' is not 0.5.")

    # R5
    d.h(3, "What the position-level scores encode")
    d.p(f"Reconstruction makes the per-protein rescaling explicit (Fig. 3a). The DCS multiplier "
        f"ranged from {V['mult_min']} to {V['mult_max']} across the {V['n_cohort_prot']} cohort "
        f"proteins, so identical products of conservation and mobility map to scores that differ "
        f"up to {float(V['mult_max']) / float(V['mult_min']):.1f}-fold between proteins. The "
        f"multiplier is set by each protein's most extreme residue: mobility is normalised by "
        f"the protein's range of RMSF, whose maximum lay between {V['rmsfmax_min']} and "
        f"{V['rmsfmax_max']} Å, and in {V['n_nd_lt_half']} of {V['n_cohort_prot']} proteins no "
        f"ClinVar-observed residue reached half of that range (median maximum normalised "
        f"mobility {V['median_max_nd']}).", "first")
    d.p(f"Whether rescaling is what makes the scores protein-specific was tested by asking a "
        f"classifier to recover the gene from each encoding (Fig. 3b; Additional file 1: Table "
        f"S8). The published scores identified the gene with balanced accuracy "
        f"{V['gid_labelled_published']} (chance 0.125; permutation null "
        f"{V['gid_labelled_published_null']}, p = {V['gid_labelled_published_p']}, the smallest "
        f"attainable with 50 permutations), and the same products before rescaling did so almost "
        f"as well ({V['gid_labelled_unrescaled']}). Conservation and normalised mobility gave "
        f"{V['gid_labelled_conservation']} and {V['gid_labelled_mobility']} individually, and "
        f"absolute RMSF {V['gid_labelled_rmsf']}. Only the within-protein percentile, which gives "
        f"every protein the same marginal distribution, removed identity altogether "
        f"({V['gid_labelled_pct_protein']}). The pattern was the same without labels across "
        f"{V['gid_n_res_all']} ClinVar-observed residues of the {V['gid_n_prot_all']} proteins "
        f"with at least ten (published {V['gid_all_clinvar_published']}, before rescaling "
        f"{V['gid_all_clinvar_unrescaled']}, percentile {V['gid_all_clinvar_pct_protein']}; "
        f"chance {1 / int(V['gid_n_prot_all']):.3f}). Rescaling therefore adds a protein "
        f"signature, as multiplying AlphaMissense residue means by the same multipliers shows "
        f"(balanced accuracy {V['gid_labelled_am_res']} raw, {V['gid_labelled_am_res_xdcs']} "
        f"multiplied), but most of the identity is present before rescaling: it is carried by "
        f"the joint distribution of conservation and mobility, which differs between proteins "
        f"whatever the scale, and min–max rescaling, which equalises only two order statistics, "
        f"cannot remove it.")
    d.p(f"Within genes, DCS and RCS ranked pathogenic above benign variants in {V['wg_dcs']} and "
        f"{V['wg_rcs']} of the {V['pairs_w']} within-gene pairs (stratified permutation "
        f"p = {V['wg_dcs_p']} and {V['wg_rcs_p']}; Fig. 3c; Additional file 1: Table S11). This "
        f"signal is that of the conservation factor ({V['wg_cons']}, p = {V['wg_cons_p']}): "
        f"normalised mobility alone gave {V['wg_nd']} (p = {V['wg_nd_p']}) and absolute RMSF "
        f"{V['wg_rmsf']}. It came almost entirely from PGK1 (DCS {V['wg_dcs_PGK1']}, "
        f"conservation {V['wg_cons_PGK1']}), was reversed in GHRHR and at chance in KAT8, and was not "
        f"supported by conditional logistic regression (DCS odds ratio per standard deviation "
        f"{V['wg_dcs_or']}, 95% CI {V['wg_dcs_orci']}; conservation {V['wg_cons_or']}, "
        f"{V['wg_cons_orci']}), although that statistic is sensitive to the skew of the DCS "
        f"scale: the normalised DCS, a monotone transform with the same within-gene AUC, gives "
        f"{V['wg_ndcs_or']}, and a rank-based model {V['wg_dcs_clrrank_or']} per decile. "
        f"AlphaMissense, by contrast, discriminated within genes "
        f"unambiguously (variant score {V['wg_amv']}, p = {V['wg_amv_p']}; odds ratio "
        f"{V['wg_amv_or']}, {V['wg_amv_orci']}), as did its residue mean ({V['wg_amr']}, "
        f"p = {V['wg_amr_p']}).")
    # R6, R7, R8, R9 and the remainder come from the second module
    from results_text2 import write_rest2
    write_rest2(d, V, FIG, rd)
