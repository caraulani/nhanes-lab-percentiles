"""Reproduce published CDC/NCHS numbers with this release's estimators.

Each target in tests/published_values.csv is recomputed from the raw NHANES
files using scripts/survey.py (the same code that builds the release).
Writes tests/validation_results.csv. Run:  uv run python -I tests/validation.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from common import CYCLE_LABEL, RAW, file_name, read_xpt  # noqa: E402
from survey import Design, proportion, var_mean, woodruff  # noqa: E402
import re  # noqa: E402
from pooling import pooled_weights  # noqa: E402


def demo(code):
    d, _ = read_xpt(RAW / CYCLE_LABEL[code] / (file_name("DEMO", code) + ".xpt"))
    return d


def lab(code, stem, cols):
    d, _ = read_xpt(RAW / CYCLE_LABEL[code] / (file_name(stem, code) + ".xpt"))
    return d[["SEQN"] + [c for c in cols if c in d.columns]]


def design_for(dm: pd.DataFrame, d: pd.DataFrame) -> Design:
    ex = dm[dm.RIDSTATR == 2]
    return Design.from_arrays(d.SDMVSTRA.astype(int).to_numpy(), d.SDMVPSU.astype(int).to_numpy(),
                              ex.SDMVSTRA.astype(int).to_numpy(), ex.SDMVPSU.astype(int).to_numpy())


# ---------------------------------------------------------------- blood lead

def lead_targets(t: pd.DataFrame) -> list[dict]:
    out = []
    for code in "GHI":
        dm = demo(code)
        pb = lab(code, "PBCD", ["LBXBPB", "WTSH2YR"])
        d = dm.merge(pb, on="SEQN")
        wcol = "WTSH2YR" if "WTSH2YR" in d.columns else "WTMEC2YR"
        d = d[d.LBXBPB.notna() & (d[wcol] > 0)].reset_index(drop=True)
        des = design_for(dm, d)
        rows = np.arange(len(d))
        y, w, age = d.LBXBPB.to_numpy(float), d[wcol].to_numpy(float), d.RIDAGEYR.to_numpy(float)
        groups = {"all ages 1+": age >= 1, "age 1-5": (age >= 1) & (age <= 5), "age 6-11": (age >= 6) & (age <= 11),
                  "age 12-19": (age >= 12) & (age <= 19), "age 20+": age >= 20}
        for r in t[(t.source_id == "cdc_lead_2011_2018") & (t.cycle == code)].itertuples():
            m = groups[r.subgroup]
            if r.statistic == "geometric_mean":
                mu, v = var_mean(des, rows[m], np.log(y[m]), w[m])
                tq = stats.t.ppf(0.975, des.df)
                est, lcl, ucl = np.exp(mu), np.exp(mu - tq * np.sqrt(v)), np.exp(mu + tq * np.sqrt(v))
            else:
                p = float(r.statistic[1:]) / 100
                q, lo, hi, _ = woodruff(des, rows[m], y[m], w[m], [p])
                est, lcl, ucl = q[0], lo[0], hi[0]
            out.append(dict(target_id=r.target_id, n=int(m.sum()), estimate=est, lcl=lcl, ucl=ucl, se=np.nan,
                            weight=wcol, design_df=des.df))
    return out


# ---------------------------------------------------------------- cholesterol (DB 515)

def chol_targets(t: pd.DataFrame) -> list[dict]:
    dm = demo("L")
    d = dm.merge(lab("L", "TCHOL", ["LBXTC", "WTPH2YR"]), on="SEQN", how="left")
    d = d.merge(lab("L", "HDL", ["LBDHDD"]), on="SEQN", how="left")
    # DB 515 sample: adults 20+ with total cholesterol and HDL-C and a positive phlebotomy weight
    d = d[(d.RIDAGEYR >= 20) & d.LBXTC.notna() & d.LBDHDD.notna() & (d.WTPH2YR > 0)].reset_index(drop=True)
    des = design_for(dm, d)
    rows, w = np.arange(len(d)), d.WTPH2YR.to_numpy(float)
    age, sex = d.RIDAGEYR.to_numpy(), d.RIAGENDR.to_numpy()
    out = []
    for r in t[t.source_id == "db515"].itertuples():
        ag, sx = r.subgroup.split(" ")
        lo, hi = {"20+": (20, 200), "20-39": (20, 39), "40-59": (40, 59), "60+": (60, 200)}[ag]
        m = (age >= lo) & (age <= hi)
        if sx == "men":
            m &= sex == 1
        if sx == "women":
            m &= sex == 2
        ind = (d.LBXTC >= 240) if r.target_id.startswith("tc240") else (d.LBDHDD < 40)
        p, se = proportion(des, rows[m], ind.to_numpy(float)[m], w[m])
        out.append(dict(target_id=r.target_id, n=int(m.sum()), estimate=100 * p, se=100 * se, lcl=np.nan,
                        ucl=np.nan, weight="WTPH2YR", design_df=des.df))
    return out


# ---------------------------------------------------------------- diabetes (DB 516)

def diabetes_targets(t: pd.DataFrame) -> list[dict]:
    dm = demo("L")
    d = dm.merge(lab("L", "GLU", ["LBXGLU", "WTSAF2YR"]), on="SEQN", how="left")
    d = d.merge(lab("L", "GHB", ["LBXGH"]), on="SEQN", how="left")
    d = d.merge(lab("L", "DIQ", ["DIQ010"]), on="SEQN", how="left")
    # DB 516: adults 20+, fasting subsample weight, pregnant women excluded
    d = d[(d.RIDAGEYR >= 20) & (d.WTSAF2YR > 0) & (d.RIDEXPRG != 1) & d.LBXGLU.notna() & d.LBXGH.notna()
          & d.DIQ010.isin([1, 2, 3])].reset_index(drop=True)
    des = design_for(dm, d)
    rows, w = np.arange(len(d)), d.WTSAF2YR.to_numpy(float)
    diag = (d.DIQ010 == 1).to_numpy()
    undiag = (~diag) & ((d.LBXGLU >= 126) | (d.LBXGH >= 6.5)).to_numpy()
    total = diag | undiag
    age, sex = d.RIDAGEYR.to_numpy(), d.RIAGENDR.to_numpy()
    out = []
    for r in t[t.source_id == "db516"].itertuples():
        ag, sx = r.subgroup.split(" ")
        lo, hi = {"20+": (20, 200), "20-39": (20, 39), "40-59": (40, 59), "60+": (60, 200)}[ag]
        m = (age >= lo) & (age <= hi)
        if sx == "men":
            m &= sex == 1
        if sx == "women":
            m &= sex == 2
        ind = {"total": total, "diag": diag, "undiag": undiag}[r.target_id.split("_")[1]]
        p, se = proportion(des, rows[m], ind.astype(float)[m], w[m])
        out.append(dict(target_id=r.target_id, n=int(m.sum()), estimate=100 * p, se=100 * se, lcl=np.nan,
                        ucl=np.nan, weight="WTSAF2YR", design_df=des.df))
    return out


# ---------------------------------------------------------------- CDC Nutrition Report (pooled cycles)

NR_FILES = {"LBXVIDMS": "VID", "LBXFER": "FERTIN", "LBDRFOSI": "FOLATE", "LBDFOTSI": "FOLFMS", "LBDB12SI": "VITB12"}
NR_YEARS = {"2011-2014": "GH", "2015-2018": "IJ"}
NR_COLS = ["mean", "p2_5", "p5", "p50", "p95", "p97_5"]


def parse_num(s):
    return float(s.replace(",", ""))


def nutrition_tables():
    """Parse tests/nutrition_report_tables.txt into one row per published number."""
    rows, cur = [], None
    for line in (HERE / "nutrition_report_tables.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        if line.startswith("@"):
            _, var, years, mtype, unit, *rest = line.split()
            cur = dict(var=var, years=years, mtype=mtype, unit=unit, females=bool(rest))
            continue
        parts = line.split("|")
        group, cells, n = parts[0], parts[1:7], parse_num(parts[7])
        for col, c in zip(NR_COLS, cells):
            m = re.match(r"([\d.,]+) \(([\d.,]+)–([\d.,]+)\)", c)
            stat = ("geometric_mean" if cur["mtype"] == "gm" else "mean") if col == "mean" else col
            rows.append(dict(cur, group=group, statistic=stat, published_value=parse_num(m.group(1)),
                             published_lcl=parse_num(m.group(2)), published_ucl=parse_num(m.group(3)),
                             published_n=int(n)))
    return pd.DataFrame(rows)


def nr_group_mask(group, age, sex, var, females):
    base = (age >= 20) if var == "LBDB12SI" else (age >= 1)
    if females:  # ferritin table: females 12-49
        base = (sex == 2) & (age >= 12) & (age <= 49)
    g = group.replace(" years", "y")
    if g.startswith("Total"):
        return base
    if g == "Male":
        return base & (sex == 1)
    if g == "Female":
        return base & (sex == 2)
    if "and older" in g:
        return base & (age >= int(re.match(r"(\d+)", g).group(1)))
    lo, hi = re.match(r"(\d+)–(\d+)y", g).groups()
    return base & (age >= int(lo)) & (age <= int(hi))


def nutrition_targets() -> pd.DataFrame:
    import weights as weights_mod
    fw = weights_mod.main().set_index(["code", "file"]).weight.to_dict()
    t = nutrition_tables()
    out = []
    for (var, years), g in t.groupby(["var", "years"], sort=False):
        parts, dms = [], []
        for code in NR_YEARS[years]:
            dm = demo(code)
            f = file_name(NR_FILES[var], code)
            x, _ = read_xpt(RAW / CYCLE_LABEL[code] / (f + ".xpt"))
            wname = fw[(code, f)]
            d = dm.merge(x[["SEQN", var] + ([wname] if wname in x.columns and wname not in dm.columns else [])], on="SEQN")
            d["w0"] = d[wname]
            d = d[(d.RIDSTATR == 2) & d[var].notna() & (d.w0 > 0)]
            d["code"] = code
            parts.append(d)
            ex = dm[dm.RIDSTATR == 2][["SDMVSTRA", "SDMVPSU"]].copy()
            ex["code"] = code
            dms.append(ex)
        d = pd.concat(parts, ignore_index=True)
        d["w"], _ = pooled_weights(d.code, d.RIAGENDR, d.RIDAGEYR, d.w0)
        dm = pd.concat(dms)
        des = Design.from_arrays((d.code + d.SDMVSTRA.astype(int).astype(str)).to_numpy(), d.SDMVPSU.astype(int).to_numpy(),
                                 (dm.code + dm.SDMVSTRA.astype(int).astype(str)).to_numpy(), dm.SDMVPSU.astype(int).to_numpy())
        rows = np.arange(len(d))
        y, w = d[var].to_numpy(float), d.w.to_numpy(float)
        age, sex = d.RIDAGEYR.to_numpy(float), d.RIAGENDR.to_numpy(float)
        tq = stats.t.ppf(0.975, des.df)
        for r in g.itertuples():
            m = nr_group_mask(r.group, age, sex, var, r.females)
            if r.statistic in ("mean", "geometric_mean"):
                yy = np.log(y[m]) if r.statistic == "geometric_mean" else y[m]
                mu, v = var_mean(des, rows[m], yy, w[m])
                lo_, hi_ = mu - tq * np.sqrt(v), mu + tq * np.sqrt(v)
                est, lcl, ucl = (np.exp(mu), np.exp(lo_), np.exp(hi_)) if r.statistic == "geometric_mean" else (mu, lo_, hi_)
            else:
                p = float(r.statistic[1:].replace("_", ".")) / 100
                q, lo, hi, _ = woodruff(des, rows[m], y[m], w[m], [p])
                est, lcl, ucl = q[0], lo[0], hi[0]
            out.append(dict(source_id="cdc_nutrition_report_2026", target_id=f"nr_{var}_{years}_{r.group}_{r.statistic}",
                            cycle=years, statistic=r.statistic, subgroup=r.group, published_value=r.published_value,
                            published_lcl=r.published_lcl, published_ucl=r.published_ucl, published_se=np.nan,
                            published_n=r.published_n, n=int(m.sum()), estimate=est, lcl=lcl, ucl=ucl, se=np.nan,
                            weight=";".join(sorted(set(fw[(c, file_name(NR_FILES[var], c))] for c in NR_YEARS[years]))),
                            design_df=des.df,
                            location=f"Table-B-{var}-{years}, row '{r.group}'"))
    res = pd.DataFrame(out)
    # pass rule (fixed before comparing): sample size identical and our estimate inside CDC's 95% CI
    res["n_match"] = res.n == res.published_n
    res["value_pass"] = (res.estimate >= res.published_lcl) & (res.estimate <= res.published_ucl)
    res["rule"] = "estimate inside the published 95% CI"
    res["difference"] = res.estimate - res.published_value
    res["rel_difference"] = res.difference / res.published_value
    return res


def tolerance(r) -> tuple[bool, str]:
    """Pass rules, fixed before comparing:
    - prevalence: |difference| <= 0.15 percentage points (published values are rounded to 0.1)
      and SE within 0.15 percentage points.
    - percentiles and geometric means: our estimate rounds to the published 3 significant
      digits within one unit in the last place, and our 95% limits are each within 10% of the
      published limits.
    - sample size: identical."""
    if r.statistic == "prevalence_pct":
        ok = abs(r.estimate - r.published_value) <= 0.15 and abs(r.se - r.published_se) <= 0.15
        return ok, "abs diff <= 0.15 pp and SE diff <= 0.15 pp"
    ulp = 10 ** (np.floor(np.log10(abs(r.published_value))) - 2)
    ok_pt = abs(float(f"{r.estimate:.3g}") - r.published_value) <= ulp + 1e-12
    ok_ci = (abs(r.lcl - r.published_lcl) <= 0.10 * r.published_lcl) and (abs(r.ucl - r.published_ucl) <= 0.10 * r.published_ucl)
    return ok_pt and ok_ci, "point within 1 unit of 3rd significant digit; CI limits within 10%"


def main() -> pd.DataFrame:
    t = pd.read_csv(HERE / "published_values.csv")
    res = pd.DataFrame(lead_targets(t) + chol_targets(t) + diabetes_targets(t))
    m = t.merge(res, on="target_id", how="left")
    m["n_match"] = m.n == m.published_n
    ok = [tolerance(r) for r in m.itertuples()]
    m["value_pass"] = [o for o, _ in ok]
    m["rule"] = [s for _, s in ok]
    m["difference"] = m.estimate - m.published_value
    m.to_csv(HERE / "validation_results.csv", index=False, float_format="%.5g")
    nr = nutrition_targets()
    nr.to_csv(HERE / "validation_nutrition_report.csv", index=False, float_format="%.5g")
    return m


if __name__ == "__main__":
    m = main()
    pd.set_option("display.width", 220)
    print(m[["target_id", "published_value", "estimate", "published_lcl", "lcl", "published_ucl", "ucl",
             "published_se", "se", "published_n", "n", "n_match", "value_pass"]].to_string(float_format=lambda x: f"{x:.4g}"))
    print("value pass:", int(m.value_pass.sum()), "/", len(m), " n match:", int(m.n_match.sum()), "/", len(m))
