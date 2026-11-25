"""Consistency checks on the built release.

Most tests read only the shipped files. Tests marked ``needs_raw`` also read the
NHANES downloads in ../_raw/nhanes (see scripts/download.py) and are skipped
when those files are absent."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from common import CYCLE_LABEL, CYCLE_YEARS, LIFEFROM, RAW, file_name, read_xpt  # noqa: E402

needs_raw = pytest.mark.skipif(
    not (RAW / "2011-2012").is_dir(),
    reason=f"raw NHANES files not found in {RAW}; run scripts/download.py to enable this test")

PC = ["p1", "p2_5", "p5", "p10", "p25", "p50", "p75", "p90", "p95", "p97_5", "p99"]
POP_FILES = {  # CDC ACS control totals (https://wwwn.cdc.gov/nchs/nhanes/ResponseRates.aspx)
    "G": "acs_totals_1112.xls", "H": "2013-2014_ACS_Control_Total.xlsx", "I": "2015-2016_acs_control_total.xlsx",
    "P": "ACS-Population-Totals-For-2017-March2020.xlsx", "L": "ACS-Population-Totals-For-August-2021-August-2023.xlsx",
}


@pytest.fixture(scope="module")
def sa():
    return pd.read_csv(ROOT / "data/percentiles_sex_age.csv")


@pytest.fixture(scope="module")
def ra():
    return pd.read_csv(ROOT / "data/percentiles_race_ethnicity.csv")


@pytest.fixture(scope="module")
def analytes():
    return pd.read_csv(ROOT / "data/analytes.csv")


def test_percentiles_monotonic(sa, ra):
    for df in (sa, ra):
        v = df[PC].to_numpy(float)
        for i in range(v.shape[0]):
            row = v[i][~np.isnan(v[i])]
            assert np.all(np.diff(row) >= 0), df.iloc[i][["analyte_id", "population", "sex"]].to_dict()


def test_cells_meet_sample_size_rules(sa, ra):
    for df in (sa, ra):
        non_median = [c for c in PC if c != "p50"]
        assert (df.loc[df[non_median].notna().any(axis=1), "n"] >= 120).all()
        assert (df.loc[df.p50.notna(), "n"] >= 30).all()
        assert (df.loc[df.p1.notna() | df.p99.notna(), "n"] >= 300).all()
        assert (df.loc[df["mean"].notna(), "n"] >= 120).all()
        assert df.loc[df.n < 30, PC + ["mean"]].isna().all().all()


def test_lod_rule(sa):
    # no published percentile at or below the weighted share below LOD
    for p, c in zip([1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99], PC):
        bad = sa[sa[c].notna() & (sa.pct_below_lod >= p)]  # NaN share (LOD unknown) never matches
        assert bad.empty, bad[["analyte_id", "age_band", "pct_below_lod", c]].head()
    assert sa.loc[sa.pct_below_lod > 50, PC[:6] + ["mean"]].isna().all().all()


def test_confidence_limits_bracket_estimates(sa):
    for p in ["p2_5", "p50", "p97_5"]:
        d = sa[sa[p].notna() & sa[p + "_lcl"].notna()]
        assert (d[p + "_lcl"] <= d[p] + 1e-9).all() and (d[p + "_ucl"] >= d[p] - 1e-9).all()


def test_sex_cells_add_up(sa):
    k = ["analyte_id", "cycle_set", "population", "age_band"]
    w = sa.pivot_table(index=k, columns="sex", values="sum_weights", aggfunc="first").dropna()
    assert np.allclose(w["all"], w["female"] + w["male"], rtol=1e-6, atol=2)
    n = sa.pivot_table(index=k, columns="sex", values="n", aggfunc="first").dropna()
    assert (n["all"] == n["female"] + n["male"]).all()


def test_reference_is_subset(sa):
    k = ["analyte_id", "cycle_set", "sex", "age_band"]
    n = sa.pivot_table(index=k, columns="population", values="n", aggfunc="first").dropna()
    assert (n["reference"] <= n["everyone"]).all()


@needs_raw
@pytest.mark.parametrize("code", list(POP_FILES))
def test_mec_weights_sum_to_population(code):
    """Sum of MEC exam weights over examined participants equals CDC's ACS
    control total for the civilian noninstitutionalized population."""
    d, _ = read_xpt(RAW / CYCLE_LABEL[code] / (file_name("DEMO", code) + ".xpt"))
    w = "WTMECPRP" if code == "P" else "WTMEC2YR"
    tot = d.loc[d.RIDSTATR == 2, w].sum()
    x = pd.read_excel(RAW / "_poptotals" / POP_FILES[code], sheet_name=0, header=None)
    row = x[x.iloc[:, 1].astype(str).str.strip() == "All Ages"].iloc[0]
    acs = float(row.iloc[2])
    assert abs(tot - acs) / acs < 0.005, (code, tot, acs)


def test_pooled_weights_are_year_averaged(sa):
    """For a full-sample MEC-weighted analyte pooled over 2011-2020P, the sum of
    scaled weights for 'all ages' cells is the years-weighted average population."""
    d = sa[(sa.analyte_id == "urine-creatinine") & (sa.cycle_set == "2011-2020P") & (sa.population == "everyone")
           & (sa.sex == "all") & sa.age_band.isin(["1-5", "6-11", "12-17", "18+"])]
    assert d.sum_weights.sum() < 330e6 and d.sum_weights.sum() > 250e6


def test_web_index_and_markers(analytes):
    idx = json.loads((ROOT / "web/index.json").read_text())
    assert set(idx) >= {"version", "generated", "source", "cycles", "markers"}
    ids = [m["analyte_id"] for m in idx["markers"]]
    assert len(ids) == len(set(ids)) == len(analytes)
    for m in idx["markers"]:
        doc = json.loads((ROOT / "web/markers" / f"{m['analyte_id']}.json").read_text())
        for k in ["analyte_id", "marker_slug", "name", "unit", "specimen", "cycles", "method_note", "populations"]:
            assert k in doc
        for pop in doc["populations"].values():
            for sex in pop.values():
                for cell in sex.values():
                    assert set(cell) == {"n", "p", "mean", "pct_below_lod", "flags", "measured_population"}
                    assert set(cell["p"]) == {"1", "2.5", "5", "10", "25", "50", "75", "90", "95", "97.5", "99"}


def test_marker_slugs_exist(analytes):
    if not (LIFEFROM / "markers.csv").exists():
        pytest.skip(f"LifeFrom release 1 markers.csv not found in {LIFEFROM}")
    m = pd.read_csv(LIFEFROM / "markers.csv")
    used = set(analytes.marker_slug.dropna())
    assert used <= set(m.marker_slug)


@needs_raw
def test_adjustment_equations_appear_in_cdc_docs():
    """Every bridging equation coefficient in config/analyte_rules.csv appears in
    the CDC document named as its source."""
    import html
    import re
    docs = {"BIOPRO_J doc bridging table": "2017-2018/BIOPRO_J.htm", "APOB_I doc (Analytic Notes)": "2015-2016/APOB_I.htm",
            "GLU_I doc": "2015-2016/GLU_I.htm", "INS_H doc": "2013-2014/INS_H.htm", "OGTT_I doc": "2015-2016/OGTT_I.htm"}
    r = pd.read_csv(ROOT / "config/analyte_rules.csv", dtype=str).fillna("")
    for x in r[r.adjust_kind != ""].itertuples():
        t = (RAW / docs[x.rule_source]).read_text(errors="replace")
        t = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))
        seg = t
        if x.rule_source.startswith("BIOPRO"):
            i = t.find(x.variable + " (")
            i = i if i >= 0 else t.find(x.variable + "1 ")
            assert i >= 0, x.variable
            seg = t[i:i + 250]
        assert x.adjust_a in seg and x.adjust_b.lstrip("-") in seg, x.variable


def test_curated_weights_match_config(analytes):
    """The weight chosen for every curated analyte and cycle (from the CDC docs,
    scripts/weights.py) equals the expected weight recorded in config/analyte_rules.csv."""
    rules = pd.read_csv(ROOT / "config/analyte_rules.csv", dtype=str).fillna("")
    got = analytes.set_index("nhanes_variable").weights.fillna("").to_dict()
    for r in rules.itertuples():
        if r.variable not in got:
            continue
        assert got[r.variable] == r.expected_weights, (r.variable, got[r.variable], r.expected_weights)


def test_ferritin_female_adult_pooled_median(sa):
    """Ferritin 2011-2020P was measured in females 12-49 only in 2015-2016
    but in everyone 12+ in 2017-March 2020. With domain-wise pooled weights the
    female 18+ median is about 65.7 ng/mL (independent estimate); a single
    cycle-level factor gave 59.1."""
    c = sa[(sa.analyte_id == "ferritin") & (sa.cycle_set == "2011-2020P") & (sa.population == "everyone")
           & (sa.sex == "female")].set_index("age_band")
    assert abs(c.loc["18+", "p50"] - 65.7) / 65.7 < 0.02
    bands = ["18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"]
    assert abs(c.loc[bands, "sum_weights"].sum() - c.loc["18+", "sum_weights"]) / c.loc["18+", "sum_weights"] < 1e-6


def test_band_weights_add_up_to_adult_cell(sa):
    """Domain-wise pooling: the adult (18+) cell's population equals the sum of its bands."""
    bands = ["18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"]
    d = sa[sa.population == "everyone"]
    k = ["analyte_id", "cycle_set", "sex"]
    b = d[d.age_band.isin(bands)].groupby(k).sum_weights.sum()
    a = d[d.age_band == "18+"].set_index(k).sum_weights
    j = pd.concat([a, b], axis=1, keys=["adult", "bands"]).dropna()
    assert np.allclose(j.adult, j.bands, rtol=1e-6, atol=10)


@needs_raw
def test_pooled_sum_weights_not_understated():
    """Blood lead ages 1-5 are not in the public 2017-March 2020 file. The pooled
    1-5 cell must estimate the years-averaged population measured in 2011-2016,
    not 6/9.2 of it."""
    sa = pd.read_csv(ROOT / "data/percentiles_sex_age.csv")
    cell = sa[(sa.analyte_id == "lead-blood") & (sa.cycle_set == "2011-2020P") & (sa.population == "everyone")
              & (sa.sex == "all") & (sa.age_band == "1-5")].iloc[0]
    tots = []
    for code in "GHI":
        d, _ = read_xpt(RAW / CYCLE_LABEL[code] / (file_name("DEMO", code) + ".xpt"))
        pb, _ = read_xpt(RAW / CYCLE_LABEL[code] / (file_name("PBCD", code) + ".xpt"))
        x = d.merge(pb, on="SEQN")
        w = "WTSH2YR" if "WTSH2YR" in x.columns else "WTMEC2YR"
        x = x[(x.RIDSTATR == 2) & x.LBXBPB.notna() & (x[w] > 0) & x.RIDAGEYR.between(1, 5)]
        tots.append(x[w].sum())
    assert abs(cell.sum_weights - np.mean(tots)) / np.mean(tots) < 0.01, (cell.sum_weights, np.mean(tots))


def test_cells_labelled_by_population_measured(sa, ra):
    """A cell whose label covers groups that were not sampled carries no percentiles,
    is flagged, and says which population was measured."""
    for df in (sa, ra):
        low = df[df.coverage_pct < 95]
        assert low["flags"].str.contains("population_not_sampled").all()
        assert low[PC + ["mean"]].isna().all().all()
    f = sa[(sa.analyte_id == "ferritin") & (sa.cycle_set == "2021-2023") & (sa.population == "everyone")
           & (sa.sex == "all") & (sa.age_band == "18+")].iloc[0]
    assert f.measured_population == "female 18-49" and pd.isna(f.p50)
    assert "population_not_sampled" in f["flags"]


def test_lod_never_silently_zero(sa, analytes):
    """pct_below_lod is published only when every cycle of the block has a CDC comment variable."""
    st = analytes.set_index("analyte_id")
    for block, col in (("2011-2020P", "lod_status_2011_2020P"), ("2021-2023", "lod_status_2021_2023")):
        d = sa[sa.cycle_set == block]
        coded = d.analyte_id.map(st[col]) == "comment variable in every cycle"
        assert d.loc[~coded, "pct_below_lod"].isna().all()
        fl = d.loc[~coded & (d.n > 0), "flags"].fillna("")
        assert fl.str.contains("lod_not_coded|lod_not_documented|lod_coded_in_some_cycles_only").all()


def test_reference_criteria_are_cycle_consistent():
    r = pd.read_csv(ROOT / "data/reference_exclusions.csv")
    default = r[(r.role == "default") & (r.criterion != "any default criterion")]
    assert set(default.criterion) == {"pregnant", "diabetes", "kidney", "smoking_cotinine", "anemia_treatment",
                                      "heart_or_stroke", "thyroid", "liver", "cancer"}
    # every criterion is assessable in every cycle
    assert (default.groupby("criterion").n_assessable.min() > 0).all()
    assert default.groupby("criterion").cycle.nunique().nunique() == 1


def test_readme_flags_match_data(sa, ra):
    """The flag table in README.md lists exactly the flags that occur in the data."""
    import re
    txt = (ROOT / "README.md").read_text()
    sec = txt[txt.index("### Flags"):]
    sec = sec[:sec.index("\n### ", 5)]
    documented = set(re.findall(r"^\| `([^`]+)` \|", sec, re.M))
    found = set()
    for df in (sa, ra):
        for x in df["flags"].dropna():
            found.update(re.sub(r"^p\d+(_5)?(?=[_<])", "pX", f) if not f.startswith("p1_p99") else f
                         for f in x.split(";"))
    assert found == documented, (found - documented, documented - found)
