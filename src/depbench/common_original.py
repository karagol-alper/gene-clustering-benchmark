"""
depbench.common
Shared utilities for the cluster-leakage benchmark analysis.

Every script in src/ imports from here so that feature engineering,
cross-validation and inference are defined in exactly one place.
"""
from __future__ import annotations
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import PowerTransformer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, average_precision_score

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────── paths
# common.py lives at <repo>/src/depbench/common.py -> three levels up is <repo>
ROOT     = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.environ.get("DEPBENCH_DATA", os.path.join(ROOT, "data"))
RESULTS  = os.environ.get("DEPBENCH_RESULTS", os.path.join(ROOT, "results"))
FIGURES  = os.path.join(RESULTS, "figures")
for d in (RESULTS, FIGURES):
    os.makedirs(d, exist_ok=True)

SEED = 11

# ─────────────────────────────────────────────── ClinVar label mapping
PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN     = {"Benign", "Likely benign", "Benign/Likely benign"}

# ──────────────────────────────────── amino-acid physicochemical table
# (Kyte-Doolittle hydrophobicity, formal charge, relative side-chain
#  volume, polarity flag, aromaticity flag, Vihinen flexibility)
AA_PROPS = {
    "A": (1.8, 0, 1.0, 0, 0, 0.36), "R": (-4.5, 1, 6.0, 1, 0, 0.53),
    "N": (-3.5, 0, 4.0, 1, 0, 0.46), "D": (-3.5, -1, 4.0, 1, 0, 0.51),
    "C": (2.5, 0, 2.0, 0, 0, 0.35), "Q": (-3.5, 0, 5.0, 1, 0, 0.49),
    "E": (-3.5, -1, 5.0, 1, 0, 0.50), "G": (-0.4, 0, 0.5, 0, 0, 0.54),
    "H": (-3.2, 0.5, 5.0, 1, 1, 0.32), "I": (4.5, 0, 4.0, 0, 0, 0.46),
    "L": (3.8, 0, 4.0, 0, 0, 0.45), "K": (-3.9, 1, 6.0, 1, 0, 0.47),
    "M": (1.9, 0, 4.0, 0, 0, 0.33), "F": (2.8, 0, 5.0, 0, 1, 0.31),
    "P": (-1.6, 0, 3.0, 0, 0, 0.51), "S": (-0.8, 0, 2.0, 1, 0, 0.51),
    "T": (-0.7, 0, 3.0, 1, 0, 0.44), "W": (-0.9, 0, 7.0, 0, 1, 0.31),
    "Y": (-1.3, 0, 6.0, 1, 1, 0.42), "V": (4.2, 0, 3.0, 0, 0, 0.40),
    "X": (0.0, 0, 3.0, 0, 0, 0.40),
}
THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

DEP_COLS = ["DCS_Score", "RCS_Score", "Normalized_DCS"]
AA_COLS = [
    "dHydro", "dCharge", "dSize", "dPolarity", "dArom", "dFlex",
    "absHydro", "absCharge", "absSize", "ref_hydro", "ref_flex",
    "alt_hydro", "alt_flex", "ref_is_hinge", "alt_is_hinge", "hinge_disr",
]
ALL_COLS = DEP_COLS + AA_COLS


def parse_protein_change(pc: str) -> tuple[str, str]:
    """'p.Thr17Met' -> ('T', 'M'). Unknown residues map to 'X'."""
    s = str(pc)
    s = s[2:] if s.startswith("p.") else s
    ref3 = s[:3]
    alt3 = "".join(c for c in s[3:] if c.isalpha())[-3:]
    return THREE_TO_ONE.get(ref3, "X"), THREE_TO_ONE.get(alt3, "X")


def add_aa_features(df: pd.DataFrame, col: str = "Protein_Change") -> pd.DataFrame:
    """Append the 16 substitution-level descriptors."""
    rows = []
    for pc in df[col]:
        r, a = parse_protein_change(pc)
        R, A = AA_PROPS[r], AA_PROPS[a]
        rows.append(dict(
            dHydro=A[0] - R[0], dCharge=A[1] - R[1], dSize=A[2] - R[2],
            dPolarity=A[3] - R[3], dArom=A[4] - R[4], dFlex=A[5] - R[5],
            absHydro=abs(A[0] - R[0]), absCharge=abs(A[1] - R[1]),
            absSize=abs(A[2] - R[2]), ref_hydro=R[0], ref_flex=R[5],
            alt_hydro=A[0], alt_flex=A[5],
            ref_is_hinge=int(r in "GP"), alt_is_hinge=int(a in "GP"),
            hinge_disr=int(r in "GP" and a not in "GP"),
        ))
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def build_cohort(path: str | None = None) -> pd.DataFrame:
    """
    Read the ClinVar-mapped score table and return the labelled benchmark
    cohort: non-conflicting germline classifications, deduplicated on
    (gene, protein change), with substitution descriptors attached.
    """
    path = path or os.path.join(DATA_DIR, "ClinVar_All_DCS_A.csv")
    raw = pd.read_csv(path)
    raw = raw.drop_duplicates(subset=["Protein_ID", "Protein_Change", "ClinVar_Accession"])
    lab = raw[raw.ClinVar_Significance.isin(PATHOGENIC | BENIGN)].copy()
    lab["label"] = lab.ClinVar_Significance.isin(PATHOGENIC).astype(int)
    lab = lab.drop_duplicates(subset=["Gene_Name", "Protein_Change"]).reset_index(drop=True)
    return add_aa_features(lab)


# ───────────────────────────────────────────────── cross-validation
def make_model(kind: str, seed: int = 0):
    if kind == "lr":
        return LogisticRegression(class_weight="balanced", max_iter=2000)
    if kind == "rf":
        return RandomForestClassifier(n_estimators=200, min_samples_leaf=3,
                                      class_weight="balanced",
                                      random_state=seed, n_jobs=-1)
    if kind == "gb":
        return GradientBoostingClassifier(n_estimators=60, max_depth=2,
                                          random_state=seed)
    raise ValueError(f"unknown model kind: {kind}")


def _preprocessor():
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("scale", PowerTransformer())])


def cross_validate(X, y, groups, split="gene", model="lr", k=9,
                   leak_scale=False, leak_select=False):
    """
    Return (out-of-fold predictions, AUC, average precision).

    split       : 'variant' (leave-one-variant-out) or 'gene' (leave-one-gene-out)
    leak_scale  : fit the imputer/transformer on the full labelled set
    leak_select : fit univariate feature selection on the full labelled set

    Setting both leak flags False and split='gene' is the only fully
    controlled configuration (C8 in the manuscript).
    """
    Xdf = pd.DataFrame(np.asarray(X))
    y = np.asarray(y)
    preds = np.zeros(len(y))

    if leak_scale:
        Xdf = pd.DataFrame(_preprocessor().fit_transform(Xdf))
    if leak_select and k and k < Xdf.shape[1]:
        base = Xdf if leak_scale else pd.DataFrame(_preprocessor().fit_transform(Xdf))
        Xdf = base.loc[:, SelectKBest(f_classif, k=k).fit(base, y).get_support()]

    cv = LeaveOneOut() if split == "variant" else LeaveOneGroupOut()
    args = (Xdf, y, groups) if split == "gene" else (Xdf, y)

    for train, test in cv.split(*args):
        Xtr, Xte = Xdf.iloc[train].values, Xdf.iloc[test].values
        if not leak_scale:
            pre = _preprocessor().fit(Xtr)
            Xtr, Xte = pre.transform(Xtr), pre.transform(Xte)
        if not leak_select and k and k < Xtr.shape[1]:
            sel = SelectKBest(f_classif, k=k).fit(Xtr, y[train])
            Xtr, Xte = sel.transform(Xtr), sel.transform(Xte)
        m = make_model(model).fit(Xtr, y[train])
        preds[test] = m.predict_proba(Xte)[:, 1]

    return preds, roc_auc_score(y, preds), average_precision_score(y, preds)


def logo_predict(df, cols, k=None, model="lr"):
    """Leave-one-gene-out predictions for a named column subset."""
    preds, _, _ = cross_validate(df[cols].values, df.label.values,
                                 df.Gene_Name.values, split="gene",
                                 model=model, k=k)
    return preds


# ───────────────────────────────────────────────────────── inference
def delong_test(y, p1, p2):
    """
    Two-sided DeLong test for two correlated ROC AUCs.
    Midrank structural components, per DeLong et al. (1988) Biometrics 44:837.
    Returns (auc1, auc2, z, p).
    """
    y = np.asarray(y)
    pos, neg = y == 1, y == 0
    m, n = pos.sum(), neg.sum()

    def components(p):
        X, Y = p[pos], p[neg]
        v01 = np.array([(np.sum(Y < x) + 0.5 * np.sum(Y == x)) / n for x in X])
        v10 = np.array([(np.sum(X > b) + 0.5 * np.sum(X == b)) / m for b in Y])
        return v01, v10

    a1, b1 = components(np.asarray(p1))
    a2, b2 = components(np.asarray(p2))
    S = np.cov(np.vstack([a1, a2])) / m + np.cov(np.vstack([b1, b2])) / n
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    z = (a1.mean() - a2.mean()) / np.sqrt(var) if var > 0 else 0.0
    return a1.mean(), a2.mean(), z, 2 * (1 - stats.norm.cdf(abs(z)))


def _fast_auc(y, score):
    """Rank-based AUC (Mann-Whitney). Equivalent to roc_auc_score but avoids
    the validation overhead, which dominates in bootstrap loops."""
    order = np.argsort(score, kind="mergesort")
    s_sorted = score[order]
    ranks = np.empty(len(score), dtype=float)
    i = 0
    r = np.arange(1, len(score) + 1, dtype=float)
    while i < len(score):
        j = i
        while j + 1 < len(score) and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        ranks[order[i:j + 1]] = r[i:j + 1].mean()
        i = j + 1
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    return (ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def bootstrap_auc(y, score, groups=None, B=10000, cluster=True, seed=SEED):
    """
    Bootstrap CI for AUC.
    cluster=True resamples whole genes (Field & Welsh 2007); this is the
    interval that respects the dependency structure. cluster=False is the
    naive variant-level bootstrap, reported only for comparison.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y); score = np.asarray(score)
    aucs = []
    if cluster:
        gs = np.unique(groups)
        # index lists cached once; rebuilding them per replicate dominates runtime
        members = [np.where(groups == g)[0] for g in gs]
        for _ in range(B):
            pick = rng.integers(0, len(gs), len(gs))
            idx = np.concatenate([members[i] for i in pick])
            yb = y[idx]
            if yb.min() == yb.max():
                continue
            aucs.append(_fast_auc(yb, score[idx]))
    else:
        n = len(y)
        for _ in range(B):
            idx = rng.integers(0, n, n)
            yb = y[idx]
            if yb.min() == yb.max():
                continue
            aucs.append(_fast_auc(yb, score[idx]))
    a = np.array(aucs)
    return np.percentile(a, 2.5), np.percentile(a, 97.5), a.std()


def gene_prevalence(df, leave_one_out=False):
    """Baseline score: the pathogenic fraction of a variant's gene."""
    if not leave_one_out:
        return df.groupby("Gene_Name").label.transform("mean").values
    return np.array([
        df.drop(i).groupby("Gene_Name").label.mean().get(
            df.Gene_Name.iloc[i], df.label.mean())
        for i in range(len(df))
    ])


def out(name: str) -> str:
    """Path helper for result files."""
    return os.path.join(RESULTS, name)
