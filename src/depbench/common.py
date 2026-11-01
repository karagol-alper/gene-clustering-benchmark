"""
depbench.common
===============
Shared definitions for the cluster-leakage benchmark analysis.

Everything that could differ between analyses (cohort construction, feature
definitions, the cross-validation engine, inference) is defined once here so
that every script uses identical code paths.
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PowerTransformer

warnings.filterwarnings("ignore")

# ──────────────────────────────────────────────────────────────── paths
# common.py lives at <repo>/src/depbench/common.py; three levels up is <repo>
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.environ.get("DEPBENCH_DATA", os.path.join(ROOT, "data"))
RESULTS = os.environ.get("DEPBENCH_RESULTS", os.path.join(ROOT, "results"))
FIGURES = os.path.join(RESULTS, "figures")
for _d in (RESULTS, FIGURES):
    os.makedirs(_d, exist_ok=True)

SEED = 11
N_PERM = int(os.environ.get("DEPBENCH_NPERM", 2000))
N_BOOT = int(os.environ.get("DEPBENCH_NBOOT", 5000))
N_JOBS = int(os.environ.get("DEPBENCH_NJOBS", os.cpu_count() or 1))

# ─────────────────────────────────────────────── ClinVar label mapping
PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}

# ─────────────────────────────────── amino-acid physicochemical table
# (Kyte-Doolittle hydrophobicity, formal charge, relative side-chain volume,
#  polarity flag, aromaticity flag, flexibility index)
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
AM_COLS = ["AM_Specific_Score", "AM_Mean_Residue_Score"]


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
            ref_aa=r, alt_aa=a,
            dHydro=A[0] - R[0], dCharge=A[1] - R[1], dSize=A[2] - R[2],
            dPolarity=A[3] - R[3], dArom=A[4] - R[4], dFlex=A[5] - R[5],
            absHydro=abs(A[0] - R[0]), absCharge=abs(A[1] - R[1]),
            absSize=abs(A[2] - R[2]), ref_hydro=R[0], ref_flex=R[5],
            alt_hydro=A[0], alt_flex=A[5],
            ref_is_hinge=int(r in "GP"), alt_is_hinge=int(a in "GP"),
            hinge_disr=int(r in "GP" and a not in "GP"),
        ))
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def read_clinvar_scores(path: str | None = None) -> pd.DataFrame:
    """All ClinVar-mapped variants with position-level scores (labelled or not)."""
    path = path or os.path.join(DATA_DIR, "ClinVar_All_DCS_A.csv")
    raw = pd.read_csv(path)
    return raw.drop_duplicates(subset=["Protein_ID", "Protein_Change", "ClinVar_Accession"])


def build_cohort(path: str | None = None) -> pd.DataFrame:
    """
    Labelled benchmark cohort: P/LP versus B/LB, deduplicated on
    (gene, protein change), with substitution descriptors attached.
    Variants of uncertain significance and conflicting classifications are
    excluded.
    """
    raw = read_clinvar_scores(path)
    lab = raw[raw.ClinVar_Significance.isin(PATHOGENIC | BENIGN)].copy()
    lab["label"] = lab.ClinVar_Significance.isin(PATHOGENIC).astype(int)
    lab = lab.drop_duplicates(subset=["Gene_Name", "Protein_Change"]).reset_index(drop=True)
    return add_aa_features(lab)


# ─────────────────────────────────────────────────── cross-validation
def make_model(kind: str, seed: int = 0, **kw):
    if kind == "lr":
        return LogisticRegression(class_weight="balanced", max_iter=2000,
                                  C=kw.get("C", 1.0))
    if kind == "rf":
        return RandomForestClassifier(n_estimators=kw.get("n_estimators", 400),
                                      min_samples_leaf=kw.get("min_samples_leaf", 3),
                                      max_features=kw.get("max_features", "sqrt"),
                                      class_weight="balanced",
                                      random_state=seed, n_jobs=kw.get("n_jobs", N_JOBS))
    if kind == "gb":
        return GradientBoostingClassifier(n_estimators=kw.get("n_estimators", 60),
                                          max_depth=kw.get("max_depth", 2),
                                          learning_rate=kw.get("learning_rate", 0.1),
                                          random_state=seed)
    raise ValueError(f"unknown model kind: {kind}")


def preprocessor():
    """Median imputation followed by Yeo-Johnson power transformation."""
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("scale", PowerTransformer())])


class FoldCache:
    """Caches label-independent preprocessing (imputation + Yeo-Johnson)
    for each training/test split of a fixed feature matrix, so that
    permutation and tuning loops refit only the label-dependent steps."""

    def __init__(self, X):
        self.X = np.asarray(X, dtype=float)
        self._c = {}
        self._full = None

    def full(self):
        if self._full is None:
            self._full = preprocessor().fit_transform(self.X)
        return self._full

    def get(self, train, test):
        key = (hash(np.asarray(train).tobytes()), hash(np.asarray(test).tobytes()))
        if key not in self._c:
            pre = preprocessor().fit(self.X[train])
            self._c[key] = (pre.transform(self.X[train]), pre.transform(self.X[test]))
        return self._c[key]


def make_folds(y, groups, split):
    y = np.asarray(y)
    if split == "variant":
        return list(LeaveOneOut().split(np.zeros(len(y))))
    return list(LeaveOneGroupOut().split(np.zeros(len(y)), y, groups))


def cross_validate(X, y, groups, split="gene", model="lr", k=9,
                   leak_scale=False, leak_select=False, model_kw=None,
                   folds=None, cache=None, engine="sklearn"):
    """
    Out-of-fold predictions under one validation configuration.

    split       : 'variant' (leave-one-variant-out) or 'gene' (leave-one-gene-out)
    leak_scale  : fit imputation + Yeo-Johnson on the complete labelled set
    leak_select : fit univariate ANOVA-F selection of k features on the
                  complete labelled set
    folds       : optional explicit list of (train, test) index arrays, which
                  overrides `split`
    cache       : optional FoldCache for the same X (label-independent
                  preprocessing is then computed once per split)
    engine      : 'sklearn', or 'fast' for the vectorised Newton solver of the
                  identical penalised logistic objective (logistic model only;
                  used inside permutation loops)

    split='gene' with both leak flags False is the only fully controlled
    configuration (C8).  Returns (predictions, AUC, average precision).
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y)
    model_kw = model_kw or {}
    preds = np.full(len(y), np.nan)
    cache = cache if cache is not None else FoldCache(X)

    Xg = cache.full() if leak_scale else X
    sel_mask = None
    if leak_select and k and k < X.shape[1]:
        sel_mask = SelectKBest(f_classif, k=k).fit(cache.full(), y).get_support()

    if folds is None:
        folds = make_folds(y, groups, split)

    for train, test in folds:
        if leak_scale:
            Xtr, Xte = Xg[train], Xg[test]
        else:
            Xtr, Xte = cache.get(train, test)
        if sel_mask is not None:
            Xtr, Xte = Xtr[:, sel_mask], Xte[:, sel_mask]
        elif k and k < Xtr.shape[1]:
            sel = SelectKBest(f_classif, k=k).fit(Xtr, y[train])
            Xtr, Xte = sel.transform(Xtr), sel.transform(Xte)
        if engine == "fast" and model == "lr":
            from depbench.fastlr import fit_lr, predict_lr
            b = fit_lr(Xtr, y[train], C=model_kw.get("C", 1.0))
            preds[test] = predict_lr(b, Xte)
        else:
            m = make_model(model, **model_kw).fit(Xtr, y[train])
            preds[test] = m.predict_proba(Xte)[:, 1]

    ok = ~np.isnan(preds)
    return preds, roc_auc_score(y[ok], preds[ok]), average_precision_score(y[ok], preds[ok])


def logo_predict(df, cols, k=None, model="lr", gene_col="Gene_Name"):
    """Leave-one-gene-out predictions for a named column subset."""
    preds, _, _ = cross_validate(df[cols].values, df.label.values,
                                 df[gene_col].values, split="gene",
                                 model=model, k=k)
    return preds


def gene_prevalence(df, leave_one_out=False, gene_col="Gene_Name"):
    """Baseline score: pathogenic fraction of a variant's own gene.
    leave_one_out=True excludes the scored variant from its gene's fraction.
    Variants that are the only member of their gene receive the cohort base
    rate, a constant that carries no ranking information."""
    y = df.label.values.astype(float)
    g = df[gene_col].values
    if not leave_one_out:
        return df.groupby(gene_col).label.transform("mean").values
    out_ = np.empty(len(df))
    base = y.mean()
    for i in range(len(df)):
        same = (g == g[i])
        same[i] = False
        out_[i] = y[same].mean() if same.any() else base
    return out_


def gene_prevalence_logo(df, gene_col="Gene_Name"):
    """Gene-prevalence baseline under leave-one-gene-out: a held-out gene's
    prevalence is unknown, so every variant receives the pathogenic rate of
    the remaining genes.  This isolates the prior-shift component of pooled
    leave-one-gene-out AUC."""
    y = df.label.values.astype(float)
    g = df[gene_col].values
    return np.array([y[g != gi].mean() for gi in g])


def out(name: str) -> str:
    """Path helper for result files."""
    return os.path.join(RESULTS, name)
