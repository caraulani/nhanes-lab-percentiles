"""Shared constants and paths."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Raw downloads live outside the repository, in ../_raw/ next to the repository folder.
RAW = Path(os.environ.get("NHANES_RAW", ROOT.parent / "_raw" / "nhanes")).expanduser()
LIFEFROM = Path(os.environ.get("LIFEFROM_RAW", ROOT.parent / "_raw" / "lifefrom-release-1")).expanduser()
BUILD = ROOT / "build"
DATA = ROOT / "data"
WEB = ROOT / "web"
CONFIG = ROOT / "config"

VERSION = "1.0.0"
RELEASE_DATE = "2026-11-25"

# code -> (label, page query, data folder year)
CYCLES = {
    "G": ("2011-2012", "CycleBeginYear=2011", "2011"),
    "H": ("2013-2014", "CycleBeginYear=2013", "2013"),
    "I": ("2015-2016", "CycleBeginYear=2015", "2015"),
    "J": ("2017-2018", "CycleBeginYear=2017", "2017"),
    "P": ("2017-2020P", "Cycle=2017-2020", "2017"),
    "L": ("2021-2023", "CycleBeginYear=2021", "2021"),
}
# Years of data collection each cycle's weights represent (CDC weighting tutorial:
# 2017-March 2020 pre-pandemic files represent 3.2 years).
CYCLE_YEARS = {"G": 2.0, "H": 2.0, "I": 2.0, "J": 2.0, "P": 3.2, "L": 2.0}
CYCLE_LABEL = {k: v[0] for k, v in CYCLES.items()}

NON_LAB_STEMS = {"DEMO", "BMX", "DIQ", "KIQ_U", "MCQ", "SMQ", "RHQ", "BPQ", "FASTQX"}


def file_stem(name: str, code: str) -> str:
    """'ALB_CR_G.xpt' -> 'ALB_CR'; 'P_ALB_CR.xpt' -> 'ALB_CR'."""
    s = name.rsplit(".", 1)[0]
    if code == "P" and s.startswith("P_"):
        return s[2:]
    suf = "_" + code
    return s[: -len(suf)] if s.endswith(suf) else s


def file_name(stem: str, code: str) -> str:
    return f"P_{stem}" if code == "P" else f"{stem}_{code}"


def read_xpt(path, metadataonly: bool = False):
    """Read a SAS transport file. Labels are decoded as Latin-1 because a few
    NHANES files carry non-UTF-8 bytes in labels."""
    import pyreadstat
    return pyreadstat.read_xport(str(path), encoding="latin1", metadataonly=metadataonly)
