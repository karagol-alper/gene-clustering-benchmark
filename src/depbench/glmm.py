"""
depbench.glmm
=============
Random-intercept logistic model fitted by maximum likelihood:

    logit P(y_ij = 1) = alpha + x_ij' beta + sigma z_j,   z_j ~ N(0, 1)

sigma is the between-gene standard deviation of the log-odds of
pathogenicity, the quantity varied in the simulations.  A profile-likelihood
interval is returned for sigma.

Each gene's integral over z is evaluated by the trapezoidal rule on a fine
uniform grid.  Gauss-Hermite rules, fixed or adaptive, are inaccurate in this
setting: when sigma is large, the integrand of a gene containing a single
class is a normal density truncated sharply at the point where the gene's
log-odds change sign, a shape that neither a fixed rule nor a rule centred at
the mode and scaled by the curvature represents well; both errors inflate
the maximised likelihood.  The trapezoidal rule is accurate for smooth,
rapidly decaying integrands provided the grid resolves the transition width
(of order 1/sigma); `check_grid` repeats the evaluation on a grid twice as
fine and twice as wide.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.special import log_expit, logsumexp

N_GRID = 96001
HALF_WIDTH = 12.0


class _Lik:
    def __init__(self, y, groups, X, n_grid=N_GRID, half_width=HALF_WIDTH):
        self.y = np.asarray(y, dtype=float)
        self.groups = np.asarray(groups)
        self.X = X
        self.z = np.linspace(-half_width, half_width, n_grid)
        dz = self.z[1] - self.z[0]
        w = np.full(n_grid, dz)
        w[0] = w[-1] = dz / 2                                  # trapezoidal weights
        self.logw = np.log(w) - 0.5 * self.z ** 2 - 0.5 * np.log(2 * np.pi)
        self.gidx = [np.where(self.groups == g)[0] for g in np.unique(self.groups)]
        if X.shape[1] == 0:                                    # clusters reduce to (k, n) patterns
            pats = {}
            for idx in self.gidx:
                key = (self.y[idx].sum(), len(idx))
                pats[key] = pats.get(key, 0) + 1
            self.patterns = [(k, n, c) for (k, n), c in pats.items()]
        else:
            self.patterns = None

    def negll(self, alpha, sigma, beta=None):
        if self.patterns is not None:
            e = alpha + sigma * self.z
            le, lne = log_expit(e), log_expit(-e)
            return -sum(c * logsumexp(k * le + (n - k) * lne + self.logw)
                        for k, n, c in self.patterns)
        eta = alpha + self.X @ beta
        tot = 0.0
        for idx in self.gidx:
            e = eta[idx][:, None] + sigma * self.z[None, :]
            yy = self.y[idx][:, None]
            tot += logsumexp((yy * log_expit(e) + (1 - yy) * log_expit(-e)).sum(0) + self.logw)
        return -tot


def fit_ri_logit(y, groups, X=None, sigma_max=100.0, n_grid=N_GRID, half_width=HALF_WIDTH):
    y = np.asarray(y, dtype=float)
    X = np.zeros((len(y), 0)) if X is None else np.asarray(X, dtype=float)
    p = X.shape[1]
    L = _Lik(y, groups, X, n_grid, half_width)
    f = lambda th: L.negll(th[0], np.exp(th[-1]), th[1:1 + p])

    best = None
    for s0 in (-1.0, 0.0, 1.0, 1.5, 2.0, 2.5, 3.0):
        th0 = np.r_[np.log(y.mean() / (1 - y.mean())) - 0.3 * np.exp(s0), np.zeros(p), s0]
        r = minimize(f, th0, method="Nelder-Mead" if p == 0 else "L-BFGS-B",
                     options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000} if p == 0 else {})
        if best is None or r.fun < best.fun:
            best = r
    th = best.x

    def prof(s):
        g = lambda a: L.negll(a[0], s, a[1:]) if p else L.negll(a[0], s)
        if p == 0:
            grid = np.linspace(-6 - 0.5 * s, 4, 60)
            a0 = grid[int(np.argmin([L.negll(a, s) for a in grid]))]
            r = minimize_scalar(lambda a: L.negll(a, s), bounds=(a0 - 0.5, a0 + 0.5),
                                method="bounded", options={"xatol": 1e-8})
            return float(r.fun)
        return float(minimize(g, th[:-1], method="L-BFGS-B").fun)

    r1 = minimize_scalar(lambda ls: prof(np.exp(ls)), bounds=(th[-1] - 1.0, th[-1] + 1.0),
                         method="bounded", options={"xatol": 1e-6})
    sigma = float(np.exp(r1.x))
    ll_max = float(min(r1.fun, best.fun))
    if p == 0:
        grid = np.linspace(-6 - 0.5 * sigma, 4, 60)
        a0 = grid[int(np.argmin([L.negll(a, sigma) for a in grid]))]
        alpha = float(minimize_scalar(lambda a: L.negll(a, sigma), bounds=(a0 - 0.5, a0 + 0.5),
                                      method="bounded", options={"xatol": 1e-8}).x)
        beta = []
    else:
        a = minimize(lambda a: L.negll(a[0], sigma, a[1:]), th[:-1], method="L-BFGS-B").x
        alpha, beta = float(a[0]), a[1:].tolist()

    crit = 3.841459 / 2
    gfun = lambda s: prof(s) - ll_max - crit
    lo = 0.0
    try:
        if gfun(1e-3) > 0:
            lo = brentq(gfun, 1e-3, sigma, xtol=1e-4)
    except ValueError:
        lo = 0.0
    try:
        hi = brentq(gfun, sigma, sigma_max, xtol=1e-3)
    except ValueError:
        hi = np.inf
    icc = lambda s: s ** 2 / (s ** 2 + np.pi ** 2 / 3) if np.isfinite(s) else 1.0
    return dict(alpha=alpha, beta=beta, sigma=sigma, sigma_lo=float(lo), sigma_hi=float(hi),
                icc_latent=float(icc(sigma)), icc_lo=float(icc(lo)), icc_hi=float(icc(hi)),
                loglik=float(-ll_max), integration=f"trapezoidal, {n_grid} points on "
                                                   f"[-{half_width:g}, {half_width:g}]")


def check_grid(y, groups, alpha, sigma):
    """Log-likelihood at (alpha, sigma) on the default grid and on a grid twice
    as fine and twice as wide (intercept-only model)."""
    X = np.zeros((len(y), 0))
    a = -_Lik(y, groups, X).negll(alpha, sigma)
    b = -_Lik(y, groups, X, 4 * (N_GRID - 1) + 1, 2 * HALF_WIDTH).negll(alpha, sigma)
    return float(a), float(b)
