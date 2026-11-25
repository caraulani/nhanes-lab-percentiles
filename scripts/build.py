"""Build every output of nhanes-lab-percentiles from the raw NHANES downloads.

    uv run python -I scripts/build.py            # full build
    uv run python -I scripts/build.py --only LBXTC,LBXGH   # a few variables (debug)

Steps
  1. catalog.py   -> build/catalog_raw.csv (every continuous lab variable)
  2. docscan.py   -> build/doc_notes.csv (method-change sentences from CDC docs)
  3. assemble analytes: curated rules (config/analyte_rules.csv) for clinical
     analytes, a default single-cycle rule for everything else
  4. person table per cycle (DEMO + reference-subset exclusion flags)
  5. survey-weighted percentiles per cell, Woodruff CIs, suppression
  6. write data/*.csv and web/*.json
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import catalog as catalog_mod  # noqa: E402
import docscan as docscan_mod  # noqa: E402
import weights as weights_mod  # noqa: E402
from pooling import coverage, pooled_weights, unit_population  # noqa: E402
from common import (BUILD, CONFIG, CYCLE_LABEL, CYCLE_YEARS, DATA, LIFEFROM, RAW,  # noqa: E402
                    RELEASE_DATE, VERSION, WEB, file_name, read_xpt)
from survey import Design, effective_n, proportion, woodruff, wquantile  # noqa: E402

PCTS = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
PCOL = {p: "p" + str(p).replace(".", "_") for p in PCTS}
CI_PCTS = [2.5, 50, 97.5]
AGE_BANDS = [("1-5", 1, 5), ("6-11", 6, 11), ("12-17", 12, 17), ("18-29", 18, 29), ("30-39", 30, 39),
             ("40-49", 40, 49), ("50-59", 50, 59), ("60-69", 60, 69), ("70-79", 70, 79), ("80+", 80, 200),
             ("18+", 18, 200)]
RACE = {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black",
        6: "Non-Hispanic Asian", 7: "Other race, including multiracial"}
SEXES = [("all", None), ("female", 2), ("male", 1)]

N_MIN_FULL = 120   # CLSI EP28: nonparametric 2.5th/97.5th need n >= 120
N_MIN_MEDIAN = 30  # median only, flagged
N_MIN_EXTREME = 300  # our rule for p1 and p99
DF_MIN = 8  # NCHS Series 2 No. 175 degrees-of-freedom flag
COVER_MIN = 0.95  # share of the labelled population that must lie in measured sex x age units

PRE_ORDER = ["P", "J", "I", "H", "G"]  # default: most recent pre-pandemic cycle with data
PRE_SET = "GHIJP"

# Files that are not part of the general NHANES sample design or cannot give
# person-level percentiles.
SMOKER_WEIGHTS = {"WTFSM", "WTALDS2Y", "WTSVS2YR"}
EXCLUDE_VARIABLES = {
    "LBXEMA": "Endomysial antibody measured only in TTG-positive samples (n < 50 per cycle).",
    "SSPHMA": "2017-2018 surplus-urine result; the same analyte is in the 2017-March 2020 file P_UVOC2 (URXPHMA).",
    "SSMUCA": "2017-2018 surplus-urine result; the same analyte is in the 2017-March 2020 file P_UVOC2 (URXMUCA).",
}
EXCLUDE_STEMS = {
    "SSTST": "2017-2018 surplus-serum sex-steroid results; CDC combined them into the 2017-March 2020 file P_TST, which is used.",
    "UCFLOW": "Urine flow rate describes the urine collection, not an analyte concentration.",
}
RENAME = {("COT", "LBXHCT"): "LBXHCOT", ("TST", "LBXTLG"): "LBXTLG", ("VITB12", "LBXB12"): "LBDB12"}
L_ALIAS = {"LBXTLG": "LBXTR"}  # 2021-2023 variable that continues an earlier analyte id
VALID_UNIT = re.compile(r"/|%|^fL$|^pg$|IU|U/L|mOsm|ratio|units", re.I)

SPECIMEN_BY_STEM = {
    "CBC": "whole blood", "GHB": "whole blood", "PBCD": "whole blood", "IHGEM": "whole blood", "CRCO": "whole blood",
    "FOLATE": "whole blood", "VOCWB": "whole blood", "VOCMWB": "whole blood", "FAR": "whole blood (washed RBC)",
    "FFMR": "whole blood (washed RBC)", "GLU": "plasma", "OGTT": "plasma", "FAS": "plasma",
}


# ----------------------------------------------------------------------------- persons

def load_persons(code: str) -> pd.DataFrame:
    label = CYCLE_LABEL[code]
    d = RAW / label
    demo, _ = read_xpt(d / (file_name("DEMO", code) + ".xpt"))
    wmec = "WTMECPRP" if code == "P" else "WTMEC2YR"
    p = pd.DataFrame({
        "SEQN": demo.SEQN.astype(int), "cycle": code, "sex": demo.RIAGENDR, "age": demo.RIDAGEYR,
        "race": demo.get("RIDRETH3"), "status": demo.RIDSTATR, "w_mec": demo[wmec],
        "stra": demo.SDMVSTRA, "psu": demo.SDMVPSU, "preg_exam": demo.get("RIDEXPRG"),
    })

    def merge(stem, cols):
        f = d / (file_name(stem, code) + ".xpt")
        if not f.exists():
            return
        df, _ = read_xpt(f)
        keep = [c for c in cols if c in df.columns]
        if keep:
            nonlocal p
            p = p.merge(df[["SEQN"] + keep].astype({"SEQN": int}), on="SEQN", how="left")

    merge("BMX", ["BMXBMI", "BMDBMIC"])
    merge("DIQ", ["DIQ010"])
    merge("KIQ_U", ["KIQ022", "KIQ025"])
    merge("MCQ", ["MCQ053", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ160M", "MCQ160L", "MCQ220"])
    merge("UCPREG", ["URXPREG"])
    merge("GHB", ["LBXGH"])
    merge("BIOPRO", ["LBXSCR"])
    merge("ALB_CR", ["URDACT"])
    merge("COT" if code != "G" else "COTNAL", ["LBXCOT"])
    if code == "L":  # person-level phlebotomy weight (CDC: use for all blood analytes)
        ph = None
        for f in sorted(d.glob("*.xpt")):
            m = read_xpt(f, metadataonly=True)[1]
            if "WTPH2YR" in m.column_names:
                df, _ = read_xpt(f)
                ph = df[["SEQN", "WTPH2YR"]] if ph is None else pd.concat([ph, df[["SEQN", "WTPH2YR"]]])
        ph = ph.dropna().drop_duplicates("SEQN")
        p = p.merge(ph.astype({"SEQN": int}), on="SEQN", how="left")
    for c in EXCL_ITEMS + ["BMDBMIC", "BMXBMI", "preg_exam"]:
        if c not in p.columns:
            p[c] = np.nan
    # serum creatinine on the 2017+ (Roche Cobas) scale: CDC forward equation (BIOPRO_J) for 2011-2016
    p["scr_bridged"] = p.LBXSCR * 0.9515 + 0.06608 if code in "GHI" else p.LBXSCR
    flags = exclusion_flags(p)
    p["excl"] = [";".join(k for k in REF_CRITERIA if r[k]) for r in flags.to_dict("records")]
    p["obese"] = flags["obesity"].to_numpy()
    for k in flags.columns:
        p["fl_" + k] = flags[k].to_numpy()
    return p


EXCL_ITEMS = ["DIQ010", "KIQ022", "KIQ025", "MCQ053", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F",
              "MCQ160M", "MCQ160L", "MCQ220", "URXPREG", "LBXGH", "LBXSCR", "URDACT", "LBXCOT"]
# Default reference-subset criteria: each uses items asked or measured the same way in every cycle
# (same variable, same age universe), so the subset means the same thing in every pooled cycle.
REF_CRITERIA = ["pregnant", "diabetes", "kidney", "smoking_cotinine", "anemia_treatment", "heart_or_stroke",
                "thyroid", "liver", "cancer"]
COTININE_CUT = 3.0  # ng/mL; Benowitz et al. 2009 (doi:10.1093/aje/kwn301): proposed overall cutpoint of 3 ng/mL


def egfr_ckd_epi_2021(scr, age, female):
    """Race-free CKD-EPI 2021 creatinine equation (Inker et al., NEJM 2021,
    doi:10.1056/NEJMoa2102953; coefficients as given by the National Kidney
    Foundation page cited in METHODS.md)."""
    k = np.where(female, 0.7, 0.9)
    a = np.where(female, -0.241, -0.302)
    r = scr / k
    e = 142 * np.minimum(r, 1) ** a * np.maximum(r, 1) ** -1.200 * 0.9938 ** age
    return np.where(female, e * 1.012, e)


def exclusion_flags(p: pd.DataFrame) -> pd.DataFrame:
    """One boolean column per criterion. A person is flagged only on positive evidence."""
    yes = lambda c: p[c] == 1  # noqa: E731
    age = p.age
    egfr = egfr_ckd_epi_2021(p.scr_bridged.astype(float), age.astype(float), p.sex == 2)
    f = pd.DataFrame({
        "pregnant": yes("preg_exam") | yes("URXPREG"),
        "diabetes": yes("DIQ010") | (p.LBXGH >= 6.5),
        "kidney": yes("KIQ022") | yes("KIQ025") | ((age >= 6) & (p.URDACT >= 30)) | ((age >= 18) & (egfr < 60)),
        "smoking_cotinine": (age >= 12) & (p.LBXCOT > COTININE_CUT),
        "anemia_treatment": yes("MCQ053"),
        "heart_or_stroke": yes("MCQ160B") | yes("MCQ160C") | yes("MCQ160D") | yes("MCQ160E") | yes("MCQ160F"),
        "thyroid": yes("MCQ160M"),
        "liver": yes("MCQ160L"),
        "cancer": yes("MCQ220"),
        # sensitivity only, not part of the default subset
        "obesity": ((age >= 18) & (p.BMXBMI >= 30)) | ((age >= 2) & (age < 18) & (p.BMDBMIC == 4)),
    })
    return f.fillna(False).astype(bool)


# ----------------------------------------------------------------------------- analytes

def lifefrom_markers() -> dict:
    m = pd.read_csv(LIFEFROM / "markers.csv")
    return dict(zip(m.marker_slug, m.canonical_name))


def kebab(s: str) -> str:
    s = s.lower().replace("µ", "u").replace("%", " percent ")
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return re.sub(r"-+", "-", s)


def assemble(cat: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    """Return (analyte table, per-cycle source table, exclusion log)."""
    excluded = []
    cat = cat[cat.variable.notna() & (cat.variable != "")].copy()
    # drop files outside the general design
    bad_stems = set()
    for (stem, code), g in cat.groupby(["stem", "code"]):
        ws = set(";".join(g.weights.fillna("").astype(str)).split(";")) - {""}
        if ws & SMOKER_WEIGHTS:
            bad_stems.add((stem, code))
            excluded.append(dict(item=f"{stem} ({CYCLE_LABEL[code]})", reason="Smokers' special sample with its own smoking weights; not part of the general population design."))
        elif stem.startswith("SS") and not ws:
            bad_stems.add((stem, code))
            excluded.append(dict(item=f"{stem} ({CYCLE_LABEL[code]})", reason="Surplus-specimen file without a sample weight; design not documented in the file."))
        elif code not in ("P",) and not g.file.iloc[0].endswith("_" + code):
            bad_stems.add((stem, code))
            excluded.append(dict(item=f"{g.file.iloc[0]} ({CYCLE_LABEL[code]})", reason="File from another cycle listed on this cycle's page."))
    cat = cat[[(s, c) not in bad_stems for s, c in zip(cat.stem, cat.code)]]
    cat["unit_ok"] = cat.unit.astype(str).str.contains(VALID_UNIT)
    for r in cat[~cat.unit_ok].drop_duplicates("variable").itertuples():
        excluded.append(dict(item=f"{r.variable} ({r.stem})", reason=f"Label has no measurement unit ('{r.unit}'): not a continuous lab result."))
    cat = cat[cat.unit_ok]
    cat = cat[~cat.variable.str.startswith("GTD")]
    curated_vars = set(pd.read_csv(CONFIG / "analyte_rules.csv").variable)
    maxd = cat.groupby("variable").n_distinct.transform("max")
    few = (maxd < 10) & ~cat.variable.isin(curated_vars)
    for r in cat[few].drop_duplicates("variable").itertuples():
        excluded.append(dict(item=f"{r.variable} ({r.stem})", reason="Fewer than 10 distinct values in every cycle: categorical or nearly constant (mostly below the LOD)."))
    cat = cat[~few]
    for s_, why in EXCLUDE_STEMS.items():
        if s_ in set(cat.stem):
            excluded.append(dict(item=f"{s_} (all cycles)", reason=why))
    cat = cat[~cat.stem.isin(EXCLUDE_STEMS)]
    cat["key"] = [RENAME.get((s, v), v) for s, v in zip(cat.stem, cat.variable)]
    cat["key"] = [L_ALIAS.get(k, k) if c == "L" else k for k, c in zip(cat.key, cat.code)]

    # SI / alternate-unit pairs within a file: same label stem, different unit
    alt_of = {}
    for (code, f), g in cat.groupby(["code", "file"], sort=False):
        for stem_lab, gg in g.groupby("label_stem", sort=False):
            if len(gg) < 2:
                continue
            names = list(gg.variable)
            nonsi = [n for n in names if not n.endswith("SI")]
            prim = nonsi[0] if nonsi else names[0]
            for n in names:
                if n != prim:
                    alt_of[(code, f, n)] = prim
        # name-based pairs: X / XSI and LBXname / LBDnameSI with differently worded labels
        vs = set(g.variable)
        for n in g.variable:
            if n.endswith("SI") and (code, f, n) not in alt_of:
                for cand in (n[:-2], "LBX" + n[3:-2], "URX" + n[3:-2]):
                    if cand in vs and cand != n and (code, f, cand) not in alt_of:
                        alt_of[(code, f, n)] = cand
                        break
    cat["is_alt"] = [(c, f, v) in alt_of for c, f, v in zip(cat.code, cat.file, cat.variable)]
    alts = cat[cat.is_alt].copy()
    alts["primary"] = [alt_of[(c, f, v)] for c, f, v in zip(alts.code, alts.file, alts.variable)]
    cat = cat[~cat.is_alt]
    for v, why in EXCLUDE_VARIABLES.items():
        if v in set(cat.variable):
            excluded.append(dict(item=v, reason=why))
    cat = cat[~cat.variable.isin(EXCLUDE_VARIABLES)]

    # one source file per (cycle, key): prefer the file with the most results
    cat = cat.sort_values("n", ascending=False).drop_duplicates(["code", "key"])
    rules = pd.read_csv(CONFIG / "analyte_rules.csv", dtype=str).fillna("")
    rules = rules.set_index("variable")
    lf = lifefrom_markers()

    rows = []
    for key, g in cat.groupby("key"):
        have = set(g.code)
        rule = rules.loc[key] if key in rules.index else None
        if rule is not None:
            pre = [c for c in rule.pre_cycles if c in have and c in PRE_SET]
            tier = "curated"
        else:
            pre = [next(c for c in PRE_ORDER if c in have)] if have & set(PRE_SET) else []
            tier = "default"
        latest = ["L"] if "L" in have else []
        first = g.sort_values("code", key=lambda s: s.map({c: i for i, c in enumerate("PJIHGL")})).iloc[0]
        stem = first.stem
        name = first.label_stem
        if rule is not None and rule.analyte_id:
            aid = rule.analyte_id
        else:
            aid = kebab(name)
        slug = rule.marker_slug if rule is not None else ""
        if slug and slug not in lf:
            raise SystemExit(f"unknown marker_slug {slug}")
        spec = rule.specimen if rule is not None and rule.specimen else (
            "urine" if key.startswith("UR") else SPECIMEN_BY_STEM.get(stem, "serum or plasma" if key.startswith("LB") else "blood"))
        rows.append(dict(key=key, analyte_id=aid, marker_slug=slug, analyte_name=name, unit=first.unit, stem=stem,
                         specimen=spec, tier=tier, pre_cycles="".join(pre), latest_cycles="".join(latest),
                         all_cycles="".join(sorted(have, key="GHIJPL".index)),
                         rule_note=rule.rule_note if rule is not None else "",
                         rule_source=rule.rule_source if rule is not None else "",
                         adjust_cycles=rule.adjust_cycles if rule is not None else "",
                         adjust_kind=rule.adjust_kind if rule is not None else "",
                         adjust_a=rule.adjust_a if rule is not None else "",
                         adjust_b=rule.adjust_b if rule is not None else ""))
    an = pd.DataFrame(rows)
    # analyte_id must be unique: disambiguate default-tier ids with the file stem
    dup = an.analyte_id.duplicated(keep=False)
    an.loc[dup & (an.tier == "default"), "analyte_id"] = [
        f"{a}-{s.lower().replace('_', '-')}" for a, s in zip(an.analyte_id[dup & (an.tier == "default")], an.stem[dup & (an.tier == "default")])]
    dup = an.analyte_id.duplicated(keep=False)
    an.loc[dup, "analyte_id"] = [f"{a}-{k.lower()}" for a, k in zip(an.analyte_id[dup], an.key[dup])]
    assert an.analyte_id.is_unique
    src = cat[["key", "code", "cycle", "file", "stem", "variable", "label", "unit", "n", "weights"]].copy()
    alts = alts.drop(columns=["key"]).merge(src[["code", "file", "variable", "key"]].rename(columns={"variable": "primary"}),
                      on=["code", "file", "primary"], how="inner")
    return an, src, excluded, alts


# ----------------------------------------------------------------------------- values

def read_values(code: str, file: str, var: str, weight: str) -> pd.DataFrame:
    """Read one analyte from one file. `weight` is the weight chosen from the
    CDC doc (scripts/weights.py). If it is not a column of the file it is
    taken from the person table (MEC exam weight or 2021-2023 WTPH2YR)."""
    df, meta = read_xpt(RAW / CYCLE_LABEL[code] / (file + ".xpt"))
    labels = dict(zip(meta.column_names, meta.column_labels))
    lc = find_lc(var, labels)
    out = pd.DataFrame({"SEQN": df.SEQN.astype(int), "value": df[var].astype(float)})
    out["lod"] = (df[lc] == 1) if lc else False
    out["w_file"] = df[weight].astype(float) if weight in df.columns else np.nan
    out.attrs["lc"] = lc or ""
    out.attrs["lod_doc"] = lod_in_doc(code, file, var)
    return out


def lod_in_doc(code: str, file: str, var: str) -> bool:
    """True if the CDC doc gives a lower limit of detection for this variable."""
    f = RAW / CYCLE_LABEL[code] / (file + ".htm")
    if not f.exists():
        return False
    t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", f.read_text(errors="replace")))
    for m in re.finditer(r"LLOD|limit of detection|detection limit", t, re.I):
        if var in t[m.start():m.start() + 3000]:
            return True
    return False


def find_lc(var: str, labels: dict) -> str | None:
    """Find the CDC comment-code variable (LOD flag) for an analyte."""
    cands = []
    if len(var) >= 3:
        base = var[:2] + "D" + var[3:]
        cands += [base + "LC", base[:-2] + "LC" if var.endswith("SI") else None, (base + "LC")[:8],
                  base[:7] + "C", base[:6] + "LC"]
    for c in cands:
        if c and c in labels and c != var and "comment" in (labels[c] or "").lower():
            return c
    # fallback: a comment variable whose label starts like the analyte label
    lab = re.sub(r"\(.*?\)\s*$", "", labels.get(var) or "").strip().lower()[:18]
    for c, l in labels.items():
        if c.endswith("LC") and l and "comment" in l.lower() and lab and l.lower().startswith(lab[:12]):
            return c
    return None


def adjust(x: np.ndarray, kind: str, a: float, b: float) -> np.ndarray:
    if kind == "linear":
        return a * x + b
    if kind == "log10linear":
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(x > 0, 10 ** (a * np.log10(x) + b), x)
    raise ValueError(kind)


# ----------------------------------------------------------------------------- estimation

def fmt(x):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return None
    return float(f"{x:.4g}")


def cell_stats(design: Design, rows: np.ndarray, y: np.ndarray, w: np.ndarray, lod: np.ndarray | None,
               mask: np.ndarray, lod_flag: str = "", cover_share: float = 1.0):
    n = int(mask.sum())
    res = dict(n=n, sum_weights=float(w[mask].sum()) if n else 0.0)
    flags = []
    if n == 0:
        return res, ["no_data"]
    yy, ww, rr = y[mask], w[mask], rows[mask]
    if lod is not None:
        ll = lod[mask]
        res["n_below_lod"] = int(ll.sum())
        # rounded as published, so the suppression rule matches the published share
        f_lod = fmt(100 * float(ww[ll].sum() / ww.sum())) / 100
        res["pct_below_lod"] = 100 * f_lod
    else:
        f_lod = 0.0
        flags.append(lod_flag)
    res["df"] = design.domain_df(rows[mask])
    if res["df"] < DF_MIN:
        flags.append("df<8")
    if cover_share < COVER_MIN:
        flags.append("population_not_sampled")
        return res, flags
    if n < N_MIN_MEDIAN:
        flags.append("n<30_suppressed")
        return res, flags
    full = n >= N_MIN_FULL
    ps = PCTS if full else [50]
    q, lo, hi, clip = woodruff(design, rr, yy, ww, [p / 100 for p in ps], df=design.df)
    for p, qq, l_, h_, c_ in zip(ps, q, lo, hi, clip):
        if p <= 100 * f_lod:
            flags.append(f"{PCOL[p]}<LOD")
            continue
        if p in (1, 99) and n < N_MIN_EXTREME:
            continue
        res[PCOL[p]] = qq
        if p in CI_PCTS:
            res[PCOL[p] + "_lcl"], res[PCOL[p] + "_ucl"] = l_, h_
            if c_:
                flags.append(f"{PCOL[p]}_ci_clipped")
    if not full:
        flags.append("median_only_n<120")
    elif n < N_MIN_EXTREME:
        flags.append("p1_p99_withheld_n<300")
    if f_lod <= 0.5 and full:
        res["mean"] = float(np.sum(ww * yy) / ww.sum())
    elif f_lod > 0.5:
        flags.append("over_50pct_below_lod")
    # effective sample size at the median (NCHS Series 2 No. 175 approach)
    if "p50" in res:
        pr, se = proportion(design, rr, (yy <= res["p50"]).astype(float), ww)
        ne = effective_n(n, pr, se)
        res["n_eff_median"] = ne
        if ne < 30:
            flags.append("n_eff<30")
    return res, flags


def estimate_block(args):
    """Estimate every cell for one analyte and one cycle set."""
    an, cycles, sources, persons_by_cycle, design_frames, adjust_spec = args
    parts = []
    weight_used, lcs, lod_doc = {}, {}, {}
    for code in cycles:
        s = sources[code]
        v = read_values(code, s["file"], s["variable"], s["weight"])
        lcs[code], lod_doc[code] = v.attrs["lc"], v.attrs["lod_doc"]
        if adjust_spec and code in adjust_spec["cycles"]:
            v["value"] = adjust(v.value.to_numpy(), adjust_spec["kind"], adjust_spec["a"], adjust_spec["b"])
        pers = persons_by_cycle[code]
        v = v.merge(pers, on="SEQN", how="inner")
        if s["weight"] in ("WTMEC2YR", "WTMECPRP"):
            v["w"] = v.w_mec
        elif s["weight"] == "WTPH2YR" and v.w_file.isna().all():
            v["w"] = v.WTPH2YR
        else:
            v["w"] = v.w_file
        weight_used[code] = s["weight"]
        v = v[(v.status == 2) & v.value.notna() & (v.w > 0) & v.age.notna() & (v.age >= 1)]
        v["code"] = code
        parts.append(v)
    d = pd.concat(parts, ignore_index=True)
    have = [c for c in cycles if (d.code == c).any()]
    d["w"], cover = pooled_weights(d.code, d.sex, d.age, d.w)
    cover_units = set(cover.index)
    pop = unit_population(persons_by_cycle, have)
    dem = pd.concat([design_frames[c] for c in have])
    design = Design.from_arrays(d.code + "_" + d.stra.astype(int).astype(str), d.psu.astype(int).to_numpy(),
                                (dem.code + "_" + dem.stra.astype(int).astype(str)).to_numpy(),
                                dem.psu.astype(int).to_numpy())
    rows = np.arange(len(d))
    y = d.value.to_numpy(float)
    w = d.w.to_numpy(float)
    # LOD: computed only when every cycle in the block carries a CDC comment (LC) variable
    coded = [c for c in have if lcs[c]]
    if len(coded) == len(have):
        lod, lod_flag, lod_status = d.lod.fillna(False).to_numpy(bool), "", "comment variable in every cycle"
    else:
        lod = None
        if coded:
            lod_flag, lod_status = "lod_coded_in_some_cycles_only", "comment variable in some cycles only"
        elif any(lod_doc[c] for c in have):
            lod_flag, lod_status = "lod_not_coded", "LLOD stated in the CDC doc; no comment variable"
        else:
            lod_flag, lod_status = "lod_not_documented", "none documented"
    age = d.age.to_numpy(float)
    sex = d.sex.to_numpy(float)
    race = d.race.to_numpy(float)
    ref = (d.excl == "").to_numpy()
    ref_bmi = ref & ~d.obese.fillna(False).to_numpy(bool)
    out = []
    label = ";".join(CYCLE_LABEL[c] for c in have)
    pops = (("everyone", np.ones(len(d), bool)), ("reference", ref), ("reference_bmi_lt30", ref_bmi))
    for pop_name, pm in pops:
        for sx, sv in SEXES:
            smask = pm if sv is None else pm & (sex == sv)
            sexvals = [1, 2] if sv is None else [sv]
            for band, lo_, hi_ in AGE_BANDS:
                m = smask & (age >= lo_) & (age <= hi_)
                if m.sum() == 0:
                    continue
                share, measured = coverage(cover_units, pop, sexvals, lo_, hi_)
                res, fl = cell_stats(design, rows, y, w, lod, m, lod_flag, share)
                cyc_cell = ";".join(CYCLE_LABEL[c] for c in have if ((d.code == c).to_numpy() & m).any())
                out.append(dict(table="sex_age", population=pop_name, sex=sx, age_band=band, race_ethnicity="",
                                cycles=cyc_cell, measured_population=measured, coverage_pct=100 * share,
                                flags=";".join(f for f in fl if f), **res))
            if pop_name == "reference_bmi_lt30":
                continue
            share, measured = coverage(cover_units, pop, sexvals, 18, 200)
            for rc, rname in RACE.items():
                m = smask & (age >= 18) & (race == rc)
                if m.sum() == 0:
                    continue
                res, fl = cell_stats(design, rows, y, w, lod, m, lod_flag, share)
                cyc_cell = ";".join(CYCLE_LABEL[c] for c in have if ((d.code == c).to_numpy() & m).any())
                out.append(dict(table="race", population=pop_name, sex=sx, age_band="18+", race_ethnicity=rname,
                                cycles=cyc_cell, measured_population=measured, coverage_pct=100 * share,
                                flags=";".join(f for f in fl if f), **res))
    info = dict(cycles=label, design_df=design.df, n_strata=design.n_strata, n_psu=design.n_psu,
                n_rows=len(d), weights=";".join(f"{CYCLE_LABEL[c]}:{weight_used[c]}" for c in have),
                lod_vars=";".join(sorted({lcs[c] for c in have} - {""})), lod_status=lod_status,
                lod_by_cycle=";".join(f"{CYCLE_LABEL[c]}:" + (f"comment {lcs[c]}" if lcs[c] else
                                      ("LLOD in doc, no comment variable" if lod_doc[c] else "none documented"))
                                      for c in have),
                variables=";".join(dict.fromkeys(sources[c]["variable"] for c in have)))
    return an["analyte_id"], label, out, info


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    print("catalog ...")
    cat = catalog_mod.main()
    print("doc scan ...")
    docscan_mod.main()
    print("weights from CDC docs ...")
    fw = weights_mod.main().set_index(["code", "file"]).weight.to_dict()
    an, src, excluded, alts = assemble(cat)
    if args.only:
        an = an[an.key.isin(args.only.split(","))]
    print(f"{len(an)} analytes")

    codes = sorted(set("".join(an.pre_cycles) + "".join(an.latest_cycles)), key="GHIJPL".index)
    persons = {c: load_persons(c) for c in codes}
    design_frames = {}
    for c, p in persons.items():
        e = p[(p.status == 2) & p.stra.notna()][["stra", "psu"]].drop_duplicates().copy()
        e["code"] = c
        design_frames[c] = e
    pcols = ["SEQN", "sex", "age", "race", "status", "w_mec", "stra", "psu", "excl", "obese", "WTPH2YR"]
    pslim = {c: p[[x for x in pcols if x in p.columns]] for c, p in persons.items()}

    jobs = []
    for a in an.itertuples(index=False):
        a = a._asdict()
        g = src[src.key == a["key"]].set_index("code")
        sources = {c: dict(file=g.loc[c, "file"], variable=g.loc[c, "variable"], weight=fw[(c, g.loc[c, "file"])])
                   for c in g.index}
        adj = None
        if a["adjust_kind"]:
            adj = dict(cycles=a["adjust_cycles"], kind=a["adjust_kind"], a=float(a["adjust_a"]), b=float(a["adjust_b"]))
        for cyc in (a["pre_cycles"], a["latest_cycles"]):
            if cyc:
                jobs.append((a, list(cyc), {c: sources[c] for c in cyc}, {c: pslim[c] for c in cyc},
                             {c: design_frames[c] for c in cyc}, adj if cyc != "L" else None))
    print(f"{len(jobs)} analyte x cycle-set blocks")
    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for i, r in enumerate(ex.map(estimate_block, jobs, chunksize=2)):
            results.append(r)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(jobs)}")
    write_outputs(an, src, alts, results, excluded, persons)
    import checksums
    checksums.main()


# ----------------------------------------------------------------------------- outputs

def alt_units(an: pd.DataFrame, alts: pd.DataFrame) -> dict:
    """alt unit and factor (alt = factor * primary) from CDC paired variables."""
    out = {}
    for key, g in alts.groupby("key"):
        r = g.sort_values("code", key=lambda s: s.map({c: i for i, c in enumerate("PJIHGL")})).iloc[0]
        a_ = an[an.key == key]
        if a_.empty:
            continue
        df, _ = read_xpt(RAW / r.cycle / (r.file + ".xpt"))
        x, z = df[r.primary].astype(float), df[r.variable].astype(float)
        ok = (x > 0) & (z > 0)
        if ok.sum() < 30:
            continue
        ratio = (z[ok] / x[ok])
        med = float(ratio.median())
        spread = float((ratio.quantile(0.95) - ratio.quantile(0.05)) / med)
        if spread <= 0.02:
            out[key] = dict(alt_unit=r.unit, alt_variable=r.variable, alt_factor=float(f"{med:.5g}"),
                            alt_factor_basis=f"median ratio {r.variable}/{r.primary} in {r.file} (5th-95th spread {spread:.2%})")
    return out


def write_outputs(an, src, alts, results, excluded, persons):
    DATA.mkdir(exist_ok=True)
    (WEB / "markers").mkdir(parents=True, exist_ok=True)
    au = alt_units(an, alts)
    rows, infos = [], []
    by_id = an.set_index("analyte_id")
    for aid, label, out, info in results:
        a = by_id.loc[aid]
        block = "2021-2023" if label == "2021-2023" else "2011-2020P"
        for r in out:
            rows.append(dict(analyte_id=aid, nhanes_variable=info["variables"],
                             marker_slug=a.marker_slug, analyte_name=a.analyte_name, specimen=a.specimen,
                             unit=a.unit, cycle_set=block, **r))
        infos.append(dict(analyte_id=aid, cycle_set=block, **info))
    res = pd.DataFrame(rows)
    pc = [PCOL[p] for p in PCTS]
    cic = [f"{PCOL[p]}_{s}" for p in CI_PCTS for s in ("lcl", "ucl")]
    cols = ["analyte_id", "nhanes_variable", "marker_slug", "analyte_name", "specimen", "unit", "cycle_set",
            "population", "sex", "age_band", "measured_population", "coverage_pct", "cycles", "n", "sum_weights", "n_below_lod", "pct_below_lod", "mean"] + pc + cic + [
            "df", "n_eff_median", "flags"]
    for c in cols:
        if c not in res.columns:
            res[c] = np.nan
    for c in ["mean", "pct_below_lod", "n_eff_median", "coverage_pct"] + pc + cic:
        res[c] = res[c].map(lambda x: fmt(x) if pd.notna(x) else None)
    res["sum_weights"] = res.sum_weights.round(0)
    sa = res[res.table == "sex_age"][cols]
    ra = res[(res.table == "race")][["race_ethnicity"] + cols]
    ra = ra[["analyte_id", "nhanes_variable", "marker_slug", "analyte_name", "specimen", "unit", "cycle_set",
             "population", "sex", "race_ethnicity"] + [c for c in cols if c not in (
                 "analyte_id", "nhanes_variable", "marker_slug", "analyte_name", "specimen", "unit", "cycle_set",
                 "population", "sex", "age_band")]]
    sa.to_csv(DATA / "percentiles_sex_age.csv", index=False)
    ra.to_csv(DATA / "percentiles_race_ethnicity.csv", index=False)

    info = pd.DataFrame(infos)
    # analytes.csv
    cat_rows = []
    for a in an.itertuples(index=False):
        g = src[src.key == a.key]
        docs = ";".join(f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{ {'G':'2011','H':'2013','I':'2015','J':'2017','P':'2017','L':'2021'}[c] }/DataFiles/{f}.htm"
                        for c, f in zip(g.code, g.file))
        i_pre = info[(info.analyte_id == a.analyte_id) & (info.cycle_set == "2011-2020P")]
        i_l = info[(info.analyte_id == a.analyte_id) & (info.cycle_set == "2021-2023")]
        alt = au.get(a.key, {})
        cat_rows.append(dict(
            analyte_id=a.analyte_id, nhanes_variable=a.key, marker_slug=a.marker_slug, analyte_name=a.analyte_name,
            specimen=a.specimen, unit=a.unit, alt_unit=alt.get("alt_unit", ""), alt_variable=alt.get("alt_variable", ""),
            alt_factor=alt.get("alt_factor", ""), alt_factor_basis=alt.get("alt_factor_basis", ""),
            tier=a.tier, cdc_files=";".join(f"{CYCLE_LABEL[c]}:{f}" for c, f in zip(g.code, g.file)),
            cycles_available=";".join(CYCLE_LABEL[c] for c in a.all_cycles),
            cycles_used_2011_2020P=i_pre.cycles.iloc[0] if len(i_pre) else "",
            cycles_used_2021_2023=i_l.cycles.iloc[0] if len(i_l) else "",
            weights=";".join(i_pre.weights.tolist() + i_l.weights.tolist()),
            design_df_2011_2020P=int(i_pre.design_df.iloc[0]) if len(i_pre) else "",
            lod_comment_variables=";".join(sorted(set(";".join(info[info.analyte_id == a.analyte_id].lod_vars).split(";")) - {""})),
            lod_status_2011_2020P=i_pre.lod_status.iloc[0] if len(i_pre) else "",
            lod_status_2021_2023=i_l.lod_status.iloc[0] if len(i_l) else "",
            lod_by_cycle=";".join(i_pre.lod_by_cycle.tolist() + i_l.lod_by_cycle.tolist()),
            method_adjustment=(f"{a.adjust_kind} a={a.adjust_a} b={a.adjust_b} applied to {';'.join(CYCLE_LABEL[c] for c in a.adjust_cycles)}" if a.adjust_kind else ""),
            method_note=a.rule_note or ("Default rule: most recent pre-pandemic cycle only; method history across cycles not reviewed."
                                        if a.tier == "default" else ""),
            method_source=a.rule_source, doc_urls=docs))
    cat_df = pd.DataFrame(cat_rows).sort_values("analyte_id")
    cat_df.to_csv(DATA / "analytes.csv", index=False)

    # coverage.csv: n per analyte x cycle (all cycles available), and use
    cov = []
    for a in an.itertuples(index=False):
        for r in src[src.key == a.key].itertuples():
            used = r.code in a.pre_cycles or r.code in a.latest_cycles
            cov.append(dict(analyte_id=a.analyte_id, nhanes_variable=r.variable, cycle=CYCLE_LABEL[r.code],
                            cdc_file=r.file, n_results=r.n, used=used,
                            cycle_set=("2021-2023" if r.code == "L" else "2011-2020P") if used else "",
                            reason_not_used="" if used else (
                                "2017-2018 is contained in 2017-March 2020 pre-pandemic files" if r.code == "J" and "P" in a.all_cycles else
                                "method change or rule; see analytes.csv method_note" if a.tier == "curated" else
                                "default rule uses only the most recent pre-pandemic cycle")))
    pd.DataFrame(cov).to_csv(DATA / "coverage.csv", index=False)

    # exclusions log and design info (build/, not shipped)
    BUILD.mkdir(exist_ok=True)
    pd.DataFrame(excluded).drop_duplicates().to_csv(BUILD / "excluded_items.csv", index=False)
    info.to_csv(BUILD / "block_info.csv", index=False)
    excl_counts = []
    items = {"pregnant": ["preg_exam", "URXPREG"], "diabetes": ["DIQ010", "LBXGH"],
             "kidney": ["KIQ022", "URDACT", "LBXSCR"], "smoking_cotinine": ["LBXCOT"], "anemia_treatment": ["MCQ053"],
             "heart_or_stroke": ["MCQ160B"], "thyroid": ["MCQ160M"], "liver": ["MCQ160L"], "cancer": ["MCQ220"],
             "obesity": ["BMXBMI", "BMDBMIC"]}
    for c, p in persons.items():
        ex = p.status == 2
        for crit in REF_CRITERIA + ["obesity"]:
            assessable = p.loc[ex, items[crit]].notna().any(axis=1)
            excl_counts.append(dict(cycle=CYCLE_LABEL[c], criterion=crit,
                                    role="sensitivity only" if crit == "obesity" else "default",
                                    n_examined=int(ex.sum()), n_assessable=int(assessable.sum()),
                                    n_flagged=int(p.loc[ex, "fl_" + crit].sum())))
        excl_counts.append(dict(cycle=CYCLE_LABEL[c], criterion="any default criterion", role="default",
                                n_examined=int(ex.sum()), n_assessable=int(ex.sum()),
                                n_flagged=int((p.loc[ex, "excl"] != "").sum())))
        excl_counts.append(dict(cycle=CYCLE_LABEL[c], criterion="any default criterion or obesity",
                                role="sensitivity only", n_examined=int(ex.sum()), n_assessable=int(ex.sum()),
                                n_flagged=int(((p.loc[ex, "excl"] != "") | p.loc[ex, "obese"]).sum())))
    pd.DataFrame(excl_counts).to_csv(DATA / "reference_exclusions.csv", index=False)

    write_web(cat_df, sa, an)
    print("wrote", len(sa), "sex/age rows,", len(ra), "race rows,", len(cat_df), "analytes")


def write_web(cat_df, sa, an):
    global LF_NAMES
    LF_NAMES = lifefrom_markers()
    def cells(df):
        pops = {}
        for r in df.itertuples(index=False):
            p = {str(x).rstrip("0").rstrip(".") if isinstance(x, float) else str(x): getattr(r, PCOL[x])
                 for x in PCTS}
            p = {("2.5" if k == "2.5" else "97.5" if k == "97.5" else k): v for k, v in p.items()}
            cell = dict(n=int(r.n), p={k: (None if v is None or (isinstance(v, float) and math.isnan(v)) else v)
                                       for k, v in p.items()},
                        mean=None if r.mean is None or (isinstance(r.mean, float) and math.isnan(r.mean)) else r.mean,
                        pct_below_lod=None if r.pct_below_lod is None or (isinstance(r.pct_below_lod, float) and math.isnan(r.pct_below_lod)) else r.pct_below_lod,
                        flags=[f for f in str(r.flags).split(";") if f and f != "nan"],
                        measured_population=r.measured_population if isinstance(r.measured_population, str) else None)
            if r.population == "reference_bmi_lt30":
                continue
            pops.setdefault(r.population, {}).setdefault(r.sex, {})[r.age_band] = cell
        return pops

    index = []
    (WEB / "markers").mkdir(parents=True, exist_ok=True)
    for f in (WEB / "markers").glob("*.json"):
        f.unlink()
    for a in cat_df.itertuples(index=False):
        d = sa[sa.analyte_id == a.analyte_id]
        pre, lat = d[d.cycle_set == "2011-2020P"], d[d.cycle_set == "2021-2023"]
        main = pre if len(pre) else lat
        if main.empty:
            continue
        name = LF_NAMES.get(a.marker_slug, a.analyte_name) if a.marker_slug and a.analyte_id == a.marker_slug else a.analyte_name
        doc = dict(analyte_id=a.analyte_id, marker_slug=a.marker_slug or None, name=name, unit=a.unit,
                   specimen=a.specimen,
                   cycles=(a.cycles_used_2011_2020P if len(pre) else a.cycles_used_2021_2023).split(";"),
                   method_note=a.method_note, populations=cells(main))
        if a.alt_unit:
            doc["alt_unit"], doc["alt_factor"] = a.alt_unit, a.alt_factor
        if len(pre) and len(lat):
            doc["cycles_2021_2023"] = dict(cycles=["2021-2023"], nhanes_variable=lat.nhanes_variable.iloc[0],
                                           populations=cells(lat))
        (WEB / "markers" / f"{a.analyte_id}.json").write_text(json.dumps(doc, separators=(",", ":"), allow_nan=False))
        ev = main[(main.population == "everyone") & (main.sex == "all")]
        n_total = int(ev[ev.age_band.isin(["1-5", "6-11", "12-17", "18+"])].n.sum())
        has_ref = bool(((main.population == "reference") & main.p50.notna()).any())
        index.append(dict(slug=a.marker_slug or None, analyte_id=a.analyte_id, name=name, unit=a.unit,
                          specimen=a.specimen, n_total=n_total, has_reference_subset=has_ref))
    idx = dict(version=VERSION, generated=RELEASE_DATE,
               source="CDC/NCHS National Health and Nutrition Examination Survey (NHANES), laboratory, demographic, examination and questionnaire public-use files",
               cycles=["2011-2012", "2013-2014", "2015-2016", "2017-2018", "2017-March 2020 pre-pandemic", "August 2021-August 2023"],
               markers=index)
    (WEB / "index.json").write_text(json.dumps(idx, indent=1, allow_nan=False))


if __name__ == "__main__":
    main()
