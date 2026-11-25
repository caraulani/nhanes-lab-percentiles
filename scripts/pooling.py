"""Pooling of survey weights across NHANES cycles, by domain.

CDC's tutorial scales each cycle's weight by years(cycle) / total years when
cycles are combined. That is only right when every cycle covers the same
population. When an analyte was measured for some sex/age groups in only some
cycles (for example ferritin: females 12-49 in one cycle, everyone 12+ in
another), a single factor over-weights the groups measured in more cycles.

We therefore scale within units u = (sex, single year of age):

    w*_i = w_i * years(c_i) / sum_{c covers u_i} years(c)

where "c covers u" means cycle c has at least one measured participant in u.
Within every unit the scaled weights then sum to the years-averaged
population of the unit over the cycles that measured it, so cell totals
(`sum_weights`) estimate the population of the cell, and no group is
over-represented. When every cycle covers every unit this reduces to CDC's
years / total-years factor.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import CYCLE_YEARS


def pooled_weights(code: pd.Series, sex: pd.Series, age: pd.Series, w: pd.Series) -> tuple[np.ndarray, pd.Series]:
    """Return scaled weights and the set of covering cycles per (sex, age) unit."""
    unit = pd.MultiIndex.from_arrays([sex.astype(int), age.astype(int)])
    df = pd.DataFrame({"code": code.to_numpy()}, index=unit)
    cover = df.groupby(level=[0, 1]).code.agg(lambda s: frozenset(s))
    years = cover.map(lambda S: sum(CYCLE_YEARS[c] for c in S))
    denom = years.reindex(unit).to_numpy()
    num = code.map(CYCLE_YEARS).to_numpy()
    return w.to_numpy(float) * num / denom, cover


def unit_population(persons: dict, codes) -> pd.Series:
    """Years-averaged population (sum of MEC weights of examined persons) per (sex, age) unit."""
    parts = []
    tot = sum(CYCLE_YEARS[c] for c in codes)
    for c in codes:
        p = persons[c]
        p = p[(p.status == 2) & p.age.notna()]
        s = p.groupby([p.sex.astype(int), p.age.astype(int)]).w_mec.sum() * (CYCLE_YEARS[c] / tot)
        parts.append(s)
    return pd.concat(parts, axis=1).sum(axis=1)


def coverage(cover_units: set, pop: pd.Series, sex_values, age_lo: int, age_hi: int) -> tuple[float, str]:
    """Share of the labelled population (sex filter x age range) in units that
    were measured, and a plain description of the measured part."""
    idx = pop.index
    sel = idx.get_level_values(0).isin(sex_values) & (idx.get_level_values(1) >= age_lo) & \
        (idx.get_level_values(1) <= age_hi)
    lab = pop[sel]
    if lab.sum() <= 0:
        return 0.0, ""
    cov = [u in cover_units for u in lab.index]
    share = float(lab[cov].sum() / lab.sum())
    measured = [u for u, c in zip(lab.index, cov) if c]
    if not measured:
        return share, "none"
    sexes = sorted({s for s, _ in measured})
    ages = [a for _, a in measured]
    sx = "female" if sexes == [2] else "male" if sexes == [1] else "all"
    lo, hi = min(ages), max(ages)
    rng = f"{lo}+" if hi >= 80 else f"{lo}-{hi}"
    return share, f"{sx} {rng}"
