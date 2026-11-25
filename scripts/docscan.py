"""Pull method-change sentences out of the NHANES .htm documentation files.

For every laboratory doc, extract sentences that mention a change of method,
instrument or laboratory, a crossover study, or a regression adjustment.
Writes build/doc_notes.csv. These sentences were read by hand when deciding
which cycles to pool (see config/cycle_rules.csv).
"""
from __future__ import annotations

import html
import re

import pandas as pd

import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from common import BUILD, CYCLES, NON_LAB_STEMS, RAW, file_stem

KEY = re.compile(
    r"crossover|cross-over|method change|changed|change in|new (instrument|method|assay|laborator)|"
    r"switched|replaced|regression|adjust|comparab|not comparable|bridg|discontinu|recalibrat|"
    r"standardi[sz]ed|different (method|instrument|laborator)|moved to|transition",
    re.I,
)


def doc_text(path) -> str:
    t = path.read_text(errors="replace")
    t = re.sub(r"<script.*?</script>|<style.*?</style>", " ", t, flags=re.S | re.I)
    t = re.sub(r"<(br|p|li|h\d|tr|div)[^>]*>", "\n", t, flags=re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    t = html.unescape(t)
    t = re.sub(r"[ \t\r\f\v]+", " ", t)
    return t


def sections(t: str) -> dict:
    """Split on the standard NHANES doc headings."""
    heads = ["Component Description", "Eligible Sample", "Description of Laboratory Methodology",
             "Laboratory Method Files", "Laboratory Quality Assurance and Monitoring",
             "Data Processing and Editing", "Analytic Notes", "Codebook and Frequencies"]
    pos = sorted((t.find(h), h) for h in heads if t.find(h) >= 0)
    out = {}
    for i, (p, h) in enumerate(pos):
        end = pos[i + 1][0] if i + 1 < len(pos) else len(t)
        out[h] = t[p:end]
    return out


def main() -> pd.DataFrame:
    rows = []
    for code, (label, _, _) in CYCLES.items():
        d = RAW / label
        if not d.exists():
            continue
        for f in sorted(d.glob("*.htm")):
            stem = file_stem(f.name, code)
            if stem in NON_LAB_STEMS:
                continue
            t = doc_text(f)
            secs = sections(t)
            body = " ".join(v for k, v in secs.items() if k != "Codebook and Frequencies") or t[:20000]
            body = re.sub(r"\s+", " ", body)
            for s in re.split(r"(?<=[.!?])\s+(?=[A-Z])", body):
                if KEY.search(s) and 20 < len(s) < 1500:
                    rows.append(dict(cycle=label, code=code, file=f.stem, stem=stem, sentence=s.strip()))
    out = pd.DataFrame(rows)
    BUILD.mkdir(exist_ok=True)
    out.to_csv(BUILD / "doc_notes.csv", index=False)
    return out


if __name__ == "__main__":
    o = main()
    print(len(o), "sentences from", o.file.nunique(), "docs")
