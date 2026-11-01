"""
Collect every number quoted in the manuscript from results/ so that the text
is generated from the analysis outputs, never typed by hand.
"""
import json
import os

import numpy as np
import pandas as pd

R = os.environ.get("DEPBENCH_RESULTS", os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "results"))


def rd(name):
    return pd.read_csv(os.path.join(R, name))


def js(name):
    return json.load(open(os.path.join(R, name)))


def f3(x):
    return f"{float(x):.3f}"


def f2(x):
    return f"{float(x):.2f}"


def pfmt(p, mc_n=None):
    p = float(p)
    if mc_n and p <= 1.0 / (mc_n + 1) + 1e-12:
        return f"≤ {1.0 / mc_n:.0e}".replace("e-0", " × 10^{−").replace("e-", " × 10^{−") + "}"
    if p <= 0:
        return "< 0.0001"
    if p < 0.001:
        e = int(np.floor(np.log10(p)))
        m = p / 10 ** e
        return f"{m:.1f} × 10^{{−{-e}}}"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.2f}" if p >= 0.1 else f"{p:.3f}"


def load():
    V = {}
    coh = rd("cohort_full.csv")
    V["n"] = len(coh)
    V["n_p"] = int(coh.label.sum())
    V["n_b"] = int((coh.label == 0).sum())
    V["n_genes"] = coh.Gene_Name.nunique()
    V["base_rate"] = f3(coh.label.mean())
    allv = rd("clinvar_all_components.csv")
    V["n_all"] = len(allv)
    V["n_unlab"] = len(allv) - len(coh)
    V["n_all_prot"] = allv.Protein_ID.nunique()
    V["n_pgk1"] = int((coh.Gene_Name == "PGK1").sum())
    V["pct_pgk1"] = f"{100 * (coh.Gene_Name == 'PGK1').mean():.0f}"

    rec = js("reconstruction_validation.json")
    maps = rd("protein_scaling_maps.csv")
    V["n_maps"] = len(maps)
    V["map_resid"] = maps.max_abs_residual.max()
    V["rec_n"] = rec["n_validated"]
    V["rec_err"] = max(rec["max_err_conservation"], rec["max_err_norm_dynamics"], rec["max_err_rmsf"])

    d = js("diagnostics.json")
    a = d["association"]
    V["chi2"] = f"{a['chi2']:.1f}"
    V["chi2_df"] = a["df"]
    V["chi2_p_mc"] = pfmt(a["p_monte_carlo"], a["n_mc"])
    V["chi2_n_mc"] = a["n_mc"]
    V["frac_exp_lt5"] = f"{100 * a['frac_expected_lt5']:.0f}"
    V["cramer"] = f2(a["cramers_v"])
    V["sigma_g"] = f"{a['ri_sigma']:.1f}"
    V["sigma_lo"] = f"{a['ri_sigma_lo']:.1f}"
    V["sigma_hi"] = f"{a['ri_sigma_hi']:.1f}" if np.isfinite(a['ri_sigma_hi']) else "∞"
    V["sigma_ci"] = f"{a['ri_sigma_lo']:.1f}–" + (f"{a['ri_sigma_hi']:.1f}" if np.isfinite(a['ri_sigma_hi']) else "∞")
    V["icc"] = f2(a["ri_icc_latent"])
    V["icc_ci"] = f"{a['ri_icc_lo']:.2f}–{a['ri_icc_hi']:.2f}"
    pr = d["pairs"]
    V["pairs"] = f"{pr['n_pairs']:,}"
    V["pairs_w"] = pr["n_within"]
    V["pairs_b"] = f"{pr['n_between']:,}"
    V["w_pct"] = f"{100 * pr['w_within']:.0f}"
    V["b_pct"] = f"{100 * (1 - pr['w_within']):.0f}"
    V["pairs_pgk1"] = pr["per_gene"]["PGK1"]
    V["pairs_ghrhr"] = pr["per_gene"]["GHRHR"]
    V["pairs_kat8"] = pr["per_gene"]["KAT8"]
    V["n_pure"] = pr["n_single_class_genes"]
    V["n_mixed"] = pr["n_mixed_genes"]
    V["n_var_pure"] = pr["n_variants_in_single_class_genes"]
    for k, v in d["prevalence"].items():
        V[f"prev_{k}"] = f3(v["auc"])
        V[f"prev_{k}_w"] = f3(v["auc_within"])
        V[f"prev_{k}_b"] = f3(v["auc_between"])

    ps = rd("permutation_summary.csv")
    for _, r in ps.iterrows():
        k = f"perm_{r.design}_{r.null}"
        V[k + "_mean"] = f3(r["mean"])
        V[k + "_sd"] = f3(r["sd"])
        V[k + "_lo"] = f3(r.pct2_5)
        V[k + "_hi"] = f3(r.pct97_5)
        V[k + "_pu"] = f3(r.p_upper)
        V[k + "_pl"] = f3(r.p_lower)
        V[k + "_p2"] = f3(r.p_two_sided) if "p_two_sided" in r else f3(min(1.0, 2*min(r.p_upper, r.p_lower)))
        V[k + "_mcse"] = f3(r.mc_se_upper)
        V[k + "_n"] = f"{int(r.n):,}"
        V[k + "_obs"] = f3(r.observed)
        V[k + "_w"] = f3(r.mean_within)
        V[k + "_b"] = f3(r.mean_between)
    ob = d["observed"]
    for k in ("C1", "C8", "POS", "POSV"):
        V[f"obs_{k}"] = f3(ob[k]["auc"])
        V[f"obs_{k}_w"] = f3(ob[k]["auc_within"])
        V[f"obs_{k}_b"] = f3(ob[k]["auc_between"])

    ab = rd("ablation_full.csv").set_index("config")
    for c in ab.index:
        for m in ("lr", "rf", "gb"):
            V[f"{c}_{m}"] = f3(ab.loc[c, m])
            V[f"{c}_{m}_ap"] = f3(ab.loc[c, f"{m}_ap"])
            V[f"{c}_{m}_w"] = f3(ab.loc[c, f"{m}_within"])
            V[f"{c}_{m}_b"] = f3(ab.loc[c, f"{m}_between"])
    V["C1_C8_lr"] = f3(ab.loc["C1", "lr"] - ab.loc["C8", "lr"])
    V["C1_C8_rf"] = f3(ab.loc["C1", "rf"] - ab.loc["C8", "rf"])
    V["C1_C8_gb"] = f3(ab.loc["C1", "gb"] - ab.loc["C8", "gb"])
    var = ab[ab.split == "variant"]
    gen = ab[ab.split == "gene"]
    V["min_var_all"] = min(var[m].min() for m in ("lr", "rf", "gb"))
    V["max_gene_all"] = max(gen[m].max() for m in ("lr", "rf", "gb"))
    V["sep_ok"] = all(var[m].min() > gen[m].max() for m in ("lr", "rf", "gb"))
    for m in ("lr", "rf", "gb"):
        V[f"{m}_w_var_rng"] = f"{var[m + '_within'].min():.3f}–{var[m + '_within'].max():.3f}"
        V[f"{m}_w_gene_rng"] = f"{gen[m + '_within'].min():.3f}–{gen[m + '_within'].max():.3f}"
        V[f"{m}_b_var_rng"] = f"{var[m + '_between'].min():.3f}–{var[m + '_between'].max():.3f}"
        V[f"{m}_b_gene_rng"] = f"{gen[m + '_between'].min():.3f}–{gen[m + '_between'].max():.3f}"

    lr16 = rd("ablation_lr16.csv").set_index("config")
    for c in lr16.index:
        V[f"{c}"] = f3(lr16.loc[c, "auc"])
        V[f"{c}_w"] = f3(lr16.loc[c, "auc_within"])
        V[f"{c}_b"] = f3(lr16.loc[c, "auc_between"])
        V[f"{c}_C"] = f"{lr16.loc[c, 'median_C']:.3g}"
    eff = rd("ablation_effects.csv").set_index("effect")
    for k in eff.index:
        V["eff_" + k.replace(" x ", "_x_")] = f"{eff.loc[k, 'auc_effect']:+.3f}"
    # main effects on the within-gene component
    for f_, hi in (("split", "variant"), ("scaling", "full set"),
                   ("selection", "full set"), ("tuning", "full set")):
        e_ = (lr16[lr16[f_] == hi].auc_within.mean()
              - lr16[lr16[f_] != hi].auc_within.mean())
        V[f"eff_{f_}_w"] = f"{e_:+.3f}"
    gl = lr16[lr16.split == "gene"]
    V["lr16_gene_w_rng"] = f"{gl.auc_within.min():.3f}–{gl.auc_within.max():.3f}"
    V["lr16_gene_b_rng"] = f"{gl.auc_between.min():.3f}–{gl.auc_between.max():.3f}"
    V["lr16_gene_rng"] = f"{gl.auc.min():.3f}–{gl.auc.max():.3f}"
    vl = lr16[lr16.split == "variant"]
    V["lr16_var_rng"] = f"{vl.auc.min():.3f}–{vl.auc.max():.3f}"
    V["lr16_var_w_rng"] = f"{vl.auc_within.min():.3f}–{vl.auc_within.max():.3f}"

    e = js("ensemble_audit.json")
    V["ens"] = f3(e["auc_reported"])
    V["ens_eq"] = f3(e["auc_equal_weights"])
    V["ens_nested"] = f3(e["auc_nested_blend"])
    V["ens_insample_best"] = f3(e["auc_in_sample_best"])
    V["ens_nested_single"] = f3(e["auc_nested_best_single"])
    V["ens_nested_gap"] = f3(e["auc_in_sample_best"] - e["auc_nested_blend"])
    V["ens_w"] = e["recovered_weights"]
    V["ens_base_rng"] = f"{min(e['auc_base_models'].values()):.3f}–{max(e['auc_base_models'].values()):.3f}"
    V["ens_within"] = f3(e["ensemble_decomposition"]["auc_within"])
    V["ens_between"] = f3(e["ensemble_decomposition"]["auc_between"])
    V["ens_ci"] = f"{e['ensemble_cluster_ci'][0]:.3f}–{e['ensemble_cluster_ci'][1]:.3f}"
    V["ens_nested_ci"] = f"{e['nested_blend_cluster_ci'][0]:.3f}–{e['nested_blend_cluster_ci'][1]:.3f}"
    V["prev_loo_ci"] = f"{e['prevalence_loo_cluster_ci'][0]:.3f}–{e['prevalence_loo_cluster_ci'][1]:.3f}"
    V["ens_vs_prev_delta"] = f3(e["ens_vs_prev"]["delta"])
    V["ens_vs_prev_p_delong"] = f2(e["ens_vs_prev"]["delong_p"])
    V["ens_vs_prev_p_clu"] = f2(e["ens_vs_prev"]["clustered_p"])
    V["ens_vs_prev_se"] = f3(e["ens_vs_prev"].get("clustered_se", float("nan")))
    V["ens_rho"] = f2(e["spearman_ens_prevalence_in_sample"])
    V["ens_within_p"] = f3(e["ensemble_within_gene"]["p_two_sided"])
    V["ens_nested_w"] = f3(e["nested_blend_decomposition"]["auc_within"])
    V["ens_nested_b"] = f3(e["nested_blend_decomposition"]["auc_between"])

    nj = js("normalisation.json")["provenance"]
    V["mult_min"] = f"{nj['dcs_mult_min']:.1f}"
    V["mult_max"] = f"{nj['dcs_mult_max']:.1f}"
    V["rcs_mult_min"] = f"{nj['rcs_mult_min']:.1f}"
    V["rcs_mult_max"] = f"{nj['rcs_mult_max']:.1f}"
    V["rmsfmax_min"] = f"{nj['rmsf_max_min']:.1f}"
    V["rmsfmax_max"] = f"{nj['rmsf_max_max']:.1f}"
    V["n_nd_lt_half"] = nj["n_max_nd_below_0_5"]
    V["median_max_nd"] = f2(nj["median_max_nd"])
    V["n_cohort_prot"] = nj["n_cohort_proteins"]
    gj = js("normalisation.json")["gene_identity"]
    V["gid_n_prot_all"] = gj["n_proteins_all"]
    V["gid_n_res_all"] = gj["n_residues_all"]
    gid = rd("gene_identity_recovery.csv")
    for _, r in gid.iterrows():
        k = f"gid_{r.set}_{r.encoding}"
        V[k] = f2(r.bacc)
        V[k + "_acc"] = f2(r.acc)
        V[k + "_p"] = pfmt(r.p)
        V[k + "_null"] = f2(r.null_mean)
    V["gid_majority"] = f2(gid[gid.set == "labelled"].majority.iloc[0])
    tr = rd("transfer_by_encoding.csv").set_index("encoding")
    for k in tr.index:
        V[f"tr_{k}"] = f3(tr.loc[k, "auc_gene"])
        V[f"tr_{k}_w"] = f3(tr.loc[k, "within_gene"])
        V[f"tr_{k}_b"] = f3(tr.loc[k, "between_gene"])
        V[f"tr_{k}_var"] = f3(tr.loc[k, "auc_variant"])
    inc = rd("incremental_over_am.csv").set_index("model")
    keymap = {"AM alone": "am", "AM + position scores (published)": "am_pub",
              "AM + position scores (unrescaled)": "am_unres", "AM + conservation": "am_cons",
              "AM + position scores (within-protein percentile)": "am_pctp",
              "AM + position scores (within-gene percentile, labelled)": "am_pctc",
              "AM + position scores (within-gene z, labelled)": "am_zc",
              "AM + position scores (within-gene centred, labelled)": "am_ctrc",
              "AM + substitution descriptors": "am_sub", "AM + full feature set": "am_full"}
    for m, k in keymap.items():
        r = inc.loc[m]
        V[f"inc_{k}"] = f3(r.auc)
        V[f"inc_{k}_ap"] = f3(r.ap)
        V[f"inc_{k}_w"] = f3(r["within"])
        V[f"inc_{k}_b"] = f3(r["between"])
        if k != "am":
            V[f"inc_{k}_d"] = f"{r.delta:+.3f}"
            V[f"inc_{k}_pd"] = f3(r.delong_p)
            V[f"inc_{k}_pc"] = f3(r.clustered_p)
            V[f"inc_{k}_zd"] = f2(r.delong_z)
            V[f"inc_{k}_zc"] = f2(r.clustered_z)
            V[f"inc_{k}_ci"] = f"{r.delta_cluster_ci_lo:+.3f} to {r.delta_cluster_ci_hi:+.3f}"
    V["inc_min_pc"] = f2(inc.clustered_p.min())
    t2 = rd("TableS2_paired_tests.csv")
    _dis = int((((t2.delong_p < 0.05) != (t2.clustered_p < 0.05))).sum())
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 13: "thirteen"}
    V["n_test_disagree"] = words.get(_dis, str(_dis))
    V["n_tests"] = words.get(len(t2), str(len(t2)))
    V["n_clustered_smaller"] = words.get(int((t2.clustered_p < t2.delong_p).sum()),
                                         str(int((t2.clustered_p < t2.delong_p).sum())))
    rc = js("rank_convention.json")
    V["inc_am_pctc_ranknative"] = f3(rc["rank_over_n"]["am_plus_pct"])
    V["inc_am_pctc_ranknative_d"] = f"{rc['rank_over_n']['delta_vs_am']:+.3f}"

    wg = rd("within_gene_discrimination.csv").set_index("score")
    names = {"AlphaMissense (variant)": "amv", "AlphaMissense (residue mean)": "amr", "DCS": "dcs",
             "RCS": "rcs", "Normalised DCS": "ndcs", "Raw DCS (before rescaling)": "rawdcs",
             "Raw RCS (before rescaling)": "rawrcs", "Conservation": "cons",
             "Normalised mobility": "nd", "Mean RMSF": "rmsf", "|delta hydrophobicity|": "hyd",
             "|delta volume|": "vol", "|delta charge|": "chg"}
    for s, k in names.items():
        r = wg.loc[s]
        V[f"wg_{k}"] = f3(r.auc_within)
        V[f"wg_{k}_ci"] = r.within_ci.replace("-", "–")
        V[f"wg_{k}_p"] = pfmt(r.within_p_two_sided)
        V[f"wg_{k}_p1"] = pfmt(r.within_p_one_sided)
        V[f"wg_{k}_b"] = f3(r.auc_between)
        V[f"wg_{k}_auc"] = f3(r.auc)
        V[f"wg_{k}_aucci"] = r.auc_cluster_ci.replace("-", "–")
        V[f"wg_{k}_or"] = f2(r.clr_or_per_sd)
        V[f"wg_{k}_orci"] = r.clr_or_ci.replace("-", "–")
        V[f"wg_{k}_clrp"] = pfmt(r.clr_p)
        for gn in ("PGK1", "GHRHR", "KAT8"):
            if f"auc_{gn}" in r and pd.notna(r[f"auc_{gn}"]):
                V[f"wg_{k}_{gn}"] = f2(r[f"auc_{gn}"])
        if "clr_rank_p" in r and pd.notna(r.clr_rank_p):
            V[f"wg_{k}_clrrank_p"] = pfmt(r.clr_rank_p)
            V[f"wg_{k}_clrrank_or"] = f2(r.clr_rank_or_per_decile)

    # ── calibration of fold-specific models (step 11) ───────────────────
    mkey = {"Position-level scores": "pos", "Full 19-feature model": "full",
            "AlphaMissense only": "am", "Training prevalence only": "prev"}
    cal = rd("calibration_by_design.csv")
    for _, r in cal.iterrows():
        k = f"cal_{mkey[r.model]}_{r.split}"
        V[k + "_slope"] = f2(r.cal_slope)
        V[k + "_int"] = f2(r.cal_intercept)
        V[k + "_brier"] = f3(r.brier)
    cj = js("calibration.json")
    V["fold_r_train"] = f2(cj["fold_mean_pred_vs_train_prevalence_r"])
    V["fold_r_gene"] = f2(cj["fold_mean_pred_vs_gene_prevalence_r"])
    V["fold_r_train_gene"] = f2(cj["train_vs_gene_prevalence_r"])
    V["n_gene_pairs"] = f"{cj['n_gene_pairs']:,}"
    fc = rd("fold_corrections.csv")
    ckey = {"none": "none", "prior shift": "prior", "fold rank": "rank"}
    for _, r in fc.iterrows():
        k = f"fx_{mkey[r.model]}_{ckey[r.correction]}"
        V[k] = f3(r.auc)
        V[k + "_w"] = f3(r.auc_within)
        V[k + "_b"] = f3(r.auc_between)
        if r.cluster_ci != "—":
            V[k + "_ci"] = r.cluster_ci.replace("-", "–")
    fcal = rd("fold_calibration.csv")
    fcal = fcal[fcal.model == "Position-level scores"].set_index("gene")
    for gn in ("HMGCS2", "GBP1", "PGK1"):
        if gn in fcal.index:
            V[f"fold_{gn}_train"] = f2(fcal.loc[gn, "train_prevalence"])
            V[f"fold_{gn}_pred"] = f2(fcal.loc[gn, "mean_prediction"])
    V["fold_pred_lo"] = f2(fcal.mean_prediction.min())
    V["fold_pred_hi"] = f2(fcal.mean_prediction.max())
    # residual gap to 0.5 after the prior-shift correction, split by component
    wgt = d["pairs"]["w_within"]
    _pp = fc[(fc.model == "Position-level scores") & (fc.correction == "prior shift")].iloc[0]
    V["fx_pos_gap"] = f3(0.5 - _pp.auc)
    V["fx_pos_gap_w"] = f3(wgt * (0.5 - _pp.auc_within))
    V["fx_pos_gap_b"] = f3((1 - wgt) * (0.5 - _pp.auc_between))
    lp = rd("leave_two_genes_out.csv").set_index("model")
    for m, k in mkey.items():
        if m in lp.index:
            V[f"lpgo_{k}"] = f3(lp.loc[m, "auc_between_lpgo"])
            V[f"lpgo_{k}_logo"] = f3(lp.loc[m, "auc_between_logo"])
            V[f"lpgo_{k}_d"] = f"{lp.loc[m, 'difference']:+.3f}"
    V["lpgo_n_pairs"] = f"{int(lp.n_pairs.iloc[0]):,}"

    # ── influence, power and paired within-gene tests (step 12) ─────────
    rj = js("robustness.json")
    V["mdes_within"] = f2(rj["min_detectable_within_auc_80"])
    V["power_065"] = f2(rj["power_at_065"])
    V["power_070"] = f2(rj["power_at_070"])
    V["null_size"] = f3(rd("within_gene_power.csv").power.iloc[0])
    V["within_null_lo"] = f3(rj["within_null_95"][0])
    V["within_null_hi"] = f3(rj["within_null_95"][1])
    V["null_df_check"] = f3(rj["null_check_max_quantile_difference"])
    V["drop_pgk1_pairs"] = rj["n_within_pairs_drop_PGK1"]
    V["drop_pgk1_dcs"] = f2(rj["within_drop_PGK1"]["within_DCS"])
    V["drop_pgk1_am"] = f2(rj["within_drop_PGK1"]["within_AlphaMissense (variant)"])
    V["drop_pgk1_cons"] = f2(rj["within_drop_PGK1"]["within_Conservation"])
    V["pos_drop_lo"] = f3(rj["pos_auc_range_over_drops"][0])
    V["pos_drop_hi"] = f3(rj["pos_auc_range_over_drops"][1])
    V["pos_drop_worst"] = rj["pos_auc_most_influential"]
    gd = rd("gene_drop_influence.csv").set_index("dropped")
    V["drop_hmgcs2_auc"] = f3(gd.loc[V["pos_drop_worst"], "pos_auc"])
    V["drop_hmgcs2_b"] = f3(gd.loc[V["pos_drop_worst"], "pos_between"])
    V["drop_hmgcs2_n"] = int(gd.loc[V["pos_drop_worst"], "n_removed"])
    pw = rd("within_gene_paired.csv")
    pkey = {"DCS vs Conservation": "dcs_cons", "RCS vs Conservation": "rcs_cons",
            "DCS vs RCS": "dcs_rcs", "Normalised mobility vs Conservation": "nd_cons",
            "AlphaMissense (variant) vs DCS": "am_dcs",
            "AlphaMissense (variant) vs Conservation": "am_cons"}
    for _, r in pw.iterrows():
        k = f"pd_{pkey[r.comparison]}"
        V[k] = f"{r.difference:+.3f}"
        V[k + "_ci"] = f"{r.ci_lo:+.3f} to {r.ci_hi:+.3f}"
    V["n_paired_sig"] = int(pw.excludes_zero.sum())
    return V
