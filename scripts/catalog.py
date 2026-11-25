"""Scan every downloaded NHANES laboratory XPT file and list candidate analyte
variables: continuous results with a unit in the SAS label.

Writes build/catalog_raw.csv (one row per cycle x file x variable).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pyreadstat

import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from common import CYCLES, NON_LAB_STEMS, RAW, BUILD, file_stem, read_xpt

UNIT_RE = re.compile(r"\(([^()]*)\)\s*$")
TRAIL_UNIT_RE = re.compile(
    r"(ug/dL|ug/L|ng/mL|ng/L|pg/mL|mg/dL|mg/L|umol/L|nmol/L|pmol/L|mmol/L|ug/g|ng/g|pg/g|fg/g|uM/L|%|IU/L|U/L|mIU/L|uIU/mL|mg/g)\s*$",
    re.I)
SKIP_LABEL = re.compile(r"comment|weight|respondent|sequence|pool|number of samples|status|"
                        r"time of venipuncture|session|hours since|fasting status|fasting hours|age group|gender|ethnicity|"
                        r"dilution|instrument|batch|run date|flag", re.I)
SKIP_NAME = re.compile(r"^(SEQN|WT|SDMV|RID|RIA|SAMPLEID|PHA|PHD|PHQ|LBDB?[A-Z0-9]*LC$)")


def unit_of(label: str) -> str | None:
    m = UNIT_RE.search(label)
    if m and re.search(r"[a-zA-Z%/]", m.group(1)) and len(m.group(1)) <= 24:
        return m.group(1).strip()
    m = TRAIL_UNIT_RE.search(label)
    return m.group(1) if m else None


def stem_of(label: str) -> str:
    s = UNIT_RE.sub("", label)
    s = TRAIL_UNIT_RE.sub("", s)
    return re.sub(r"\s+", " ", s).strip(" ,-")


def main() -> pd.DataFrame:
    rows = []
    for code, (label, _, _) in CYCLES.items():
        d = RAW / label
        if not d.exists():
            continue
        for f in sorted(d.glob("*.xpt")):
            stem = file_stem(f.name, code)
            if stem in NON_LAB_STEMS:
                continue
            df, meta = read_xpt(f)
            if "SEQN" not in df.columns:
                rows.append(dict(cycle=label, code=code, file=f.stem, stem=stem, variable="", label="",
                                 unit="", n=0, weights="", note="no SEQN (pooled samples)"))
                continue
            weights = [c for c in df.columns if c.startswith("WT")]
            labels = dict(zip(meta.column_names, meta.column_labels))
            for v in df.columns:
                lab = labels.get(v) or ""
                if SKIP_NAME.search(v) or SKIP_LABEL.search(lab):
                    continue
                u = unit_of(lab)
                if not u:
                    continue
                s = df[v]
                if not pd.api.types.is_numeric_dtype(s):
                    continue
                n = int(s.notna().sum())
                if n == 0:
                    continue
                rows.append(dict(cycle=label, code=code, file=f.stem, stem=stem, variable=v, label=lab,
                                 unit=u, label_stem=stem_of(lab), n=n, weights=";".join(weights),
                                 n_distinct=int(s.nunique()), note=""))
    out = pd.DataFrame(rows)
    BUILD.mkdir(exist_ok=True)
    out.to_csv(BUILD / "catalog_raw.csv", index=False)
    return out


if __name__ == "__main__":
    o = main()
    print(o.groupby("cycle").variable.count())
