"""Design-based estimators for NHANES: weighted quantiles, Woodruff confidence
limits, Taylor-linearised variance for domain means and proportions.

Conventions
-----------
* Weighted empirical CDF: F(t) = sum_{i: y_i <= t} w_i / sum_i w_i, over the domain.
* Quantile Q(p), default method "step": the smallest observed value v with
  F(v) >= p (the inverse of the weighted empirical CDF, no interpolation).
  This reproduced all 20 published percentiles of the CDC blood lead table we
  tested exactly (tests/validation.py), so it matches CDC practice.
  Optional method "interpolate" (the SAS/STAT SURVEYMEANS documentation,
  "Quantiles"): with distinct values v_1 < ... < v_m and F_k = F(v_k),
  if p < F_1, Q = v_1; if F_k <= p < F_{k+1}, Q = v_k + (p - F_k) / (F_{k+1} - F_k)
  * (v_{k+1} - v_k); if p >= F_m, Q = v_m.
* Variance: first-stage with-replacement Taylor linearisation using strata and
  PSUs (NHANES SDMVSTRA / SDMVPSU). For a domain estimator, out-of-domain rows
  contribute zero but the full design (all PSUs in every stratum) is kept, so
  the variance is the domain (subpopulation) variance, not a subset variance.
* Woodruff limits: SE of F at Q(p) from the linearised variance; limits on the
  probability scale p -/+ t_{df,0.975} * SE are mapped back through Q(.).
  Probability limits outside [0, 1] are clipped and the cell is flagged.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


def wquantile(y: np.ndarray, w: np.ndarray, ps) -> np.ndarray:
    """Weighted quantiles (see module docstring for the definition)."""
    ps = np.atleast_1d(np.asarray(ps, dtype=float))
    order = np.argsort(y, kind="mergesort")
    ys, ws = y[order], w[order]
    vals, start = np.unique(ys, return_index=True)
    wsum = np.add.reduceat(ws, start)
    F = np.cumsum(wsum) / wsum.sum()
    F[-1] = 1.0
    return _q_from_cdf(vals, F, ps)


QUANTILE_METHOD = "step"
_EPS = 1e-12


def _q_from_cdf(vals: np.ndarray, F: np.ndarray, ps: np.ndarray) -> np.ndarray:
    if QUANTILE_METHOD == "step":
        idx = np.searchsorted(F, np.asarray(ps) - _EPS, side="left")
        return vals[np.minimum(idx, len(vals) - 1)]
    m = len(vals)
    idx = np.searchsorted(F, ps, side="right")  # number of F_k <= p
    out = np.empty(len(ps))
    for j, (p, k) in enumerate(zip(ps, idx)):
        if k == 0:
            out[j] = vals[0]
        elif k >= m:
            out[j] = vals[-1]
        else:
            lo, hi = k - 1, k
            out[j] = vals[lo] + (p - F[lo]) / (F[hi] - F[lo]) * (vals[hi] - vals[lo])
    return out


@dataclass
class Design:
    """Strata/PSU structure for a set of rows.

    psu: integer PSU id per row, unique across strata (0..n_psu-1)
    psu_stratum: stratum index for each PSU id
    """
    psu: np.ndarray
    psu_stratum: np.ndarray
    n_psu: int
    n_strata: int

    @classmethod
    def from_arrays(cls, stratum_key: np.ndarray, psu_key: np.ndarray,
                    all_stratum: np.ndarray, all_psu: np.ndarray) -> "Design":
        """stratum_key/psu_key: per analysis row. all_*: every (stratum, psu)
        pair in the full design (e.g. from DEMO), so PSUs with no analysis
        rows still count."""
        pairs = sorted(set(zip(all_stratum.tolist(), all_psu.tolist())))
        pid = {p: i for i, p in enumerate(pairs)}
        strata = sorted({s for s, _ in pairs})
        sid = {s: i for i, s in enumerate(strata)}
        psu_stratum = np.array([sid[s] for s, _ in pairs])
        psu = np.array([pid[(s, p)] for s, p in zip(stratum_key.tolist(), psu_key.tolist())])
        return cls(psu=psu, psu_stratum=psu_stratum, n_psu=len(pairs), n_strata=len(strata))

    @property
    def df(self) -> int:
        return self.n_psu - self.n_strata

    def domain_df(self, mask: np.ndarray) -> int:
        psus = np.unique(self.psu[mask])
        return len(psus) - len(np.unique(self.psu_stratum[psus]))

    def var_total(self, z: np.ndarray, rows: np.ndarray | None = None) -> float:
        """With-replacement variance of a linearised total sum_i z_i.
        z is aligned to `rows` (indices into the design rows) or to all rows.
        V = sum_h n_h/(n_h-1) * sum_j (z_hj - mean_h)^2 over PSUs j in stratum h."""
        psu = self.psu if rows is None else self.psu[rows]
        zp = np.bincount(psu, weights=z, minlength=self.n_psu)
        nh = np.bincount(self.psu_stratum, minlength=self.n_strata).astype(float)
        s1 = np.bincount(self.psu_stratum, weights=zp, minlength=self.n_strata)
        s2 = np.bincount(self.psu_stratum, weights=zp * zp, minlength=self.n_strata)
        ok = nh >= 2
        ss = s2[ok] - s1[ok] ** 2 / nh[ok]
        return float(np.sum(nh[ok] / (nh[ok] - 1) * np.maximum(ss, 0.0)))


def var_mean(design: Design, rows: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """Domain mean (ratio estimator) and its linearised variance.
    rows/y/w: domain members only (out-of-domain rows contribute zero)."""
    W = w.sum()
    mu = float(np.sum(w * y) / W)
    z = w * (y - mu) / W
    return mu, design.var_total(z, rows)


def woodruff(design: Design, rows: np.ndarray, y: np.ndarray, w: np.ndarray, ps, alpha: float = 0.05,
             df: int | None = None):
    """Point estimates and nonsymmetric Woodruff limits for quantiles.

    Returns arrays (q, lcl, ucl, clipped) for each p.
    """
    ps = np.atleast_1d(np.asarray(ps, dtype=float))
    order = np.argsort(y, kind="mergesort")
    ys, ws, rs = y[order], w[order], rows[order]
    vals, start = np.unique(ys, return_index=True)
    wsum = np.add.reduceat(ws, start)
    W = ws.sum()
    F = np.cumsum(wsum) / W
    F[-1] = 1.0
    q = _q_from_cdf(vals, F, ps)
    t = stats.t.ppf(1 - alpha / 2, df if df is not None else design.df)
    lcl = np.full(len(ps), np.nan)
    ucl = np.full(len(ps), np.nan)
    clipped = np.zeros(len(ps), dtype=bool)
    for j, (p, qj) in enumerate(zip(ps, q)):
        ind = (ys <= qj).astype(float)
        Fq = float(np.sum(ws * ind) / W)
        z = ws * (ind - Fq) / W
        se = np.sqrt(design.var_total(z, rs))
        pl, pu = p - t * se, p + t * se
        if pl < 0 or pu > 1:
            clipped[j] = True
        pl, pu = max(pl, 0.0), min(pu, 1.0)
        lcl[j], ucl[j] = _q_from_cdf(vals, F, np.array([pl, pu]))
    return q, lcl, ucl, clipped


def proportion(design: Design, rows: np.ndarray, ind: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """Domain proportion and its standard error (linearised)."""
    p, v = var_mean(design, rows, ind.astype(float), w)
    return p, float(np.sqrt(v))


def effective_n(n: int, p: float, se: float) -> float:
    """n / design effect for a proportion (NCHS Series 2 No. 175 approach)."""
    if se <= 0 or p <= 0 or p >= 1:
        return float(n)
    deff = (se ** 2) / (p * (1 - p) / n)
    return n / deff if deff > 0 else float(n)
