"""Choose the survey weight for every laboratory file from what its CDC doc says.

For each laboratory file we read the documentation text outside the codebook,
collect the sentences that mention weights, and decide:

1. The doc names a weight variable that is in the file -> that variable
   ("doc names <WT> and it is in the file").
2. The doc names a weight variable that is not in the file (for example the
   2021-2023 phlebotomy weight WTPH2YR, or WTMEC2YR) -> that variable, taken
   from DEMO or another file of the same cycle.
3. The doc names no variable but says the MEC exam / examined sample weights
   should be used -> WTMEC2YR (WTMECPRP for 2017-March 2020).
4. The doc says nothing usable and the file has exactly one weight column ->
   that column (basis recorded as such).
5. Otherwise -> MEC exam weight, basis "default: doc names no weight".

Writes build/file_weights.csv: cycle, file, weight, basis, quote.
"""
from __future__ import annotations

import html
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BUILD, CYCLE_LABEL, CYCLES, NON_LAB_STEMS, RAW, file_stem, read_xpt  # noqa: E402

WT = re.compile(r"\bWT[A-Z0-9]{3,8}\b")


def doc_body(path: Path) -> str:
    t = path.read_text(errors="replace")
    t = re.sub(r"<script.*?</script>|<style.*?</style>", " ", t, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    t = re.sub(r"\s+", " ", html.unescape(t))
    i = t.find("Codebook and Frequencies", 3000)
    return t[:i] if i > 0 else t


def choose(code: str, cols: list[str], body: str, blood: bool = False) -> tuple[str, str, str]:
    mec = "WTMECPRP" if code == "P" else "WTMEC2YR"
    sents = [s for s in re.split(r"(?<=[.!?])\s+", body) if re.search(r"weight", s, re.I)]
    # drop the codebook-style list at the top of the doc ("SEQN - ... WTxxx - ...")
    sents = [s for s in sents if "Respondent sequence number" not in s and " - " not in s[:60]]
    named = []
    for s in sents:
        for w in WT.findall(s):
            if w not in named:
                named.append(w)
    quote = " ".join(sents)[:400]
    in_file = [w for w in named if w in cols]
    if in_file:
        sub = [w for w in in_file if w != "WTPH2YR"]
        return (sub or in_file)[0], f"doc names {(sub or in_file)[0]} and it is in the file", quote
    file_specific = re.search(r"(weights? (are |is )?(included )?in this (data )?file|special sampl\w* weights in this file|"
                              r"specific sample weights for this subsample)", body, re.I)
    subsample_cols = [c for c in cols if c != "WTPH2YR"]
    misnamed = [w for w in named if not w.startswith(("WTMEC", "WTPH", "WTINT"))]
    if misnamed and len(subsample_cols) == 1:
        import difflib
        r = difflib.SequenceMatcher(None, misnamed[0], subsample_cols[0]).ratio()
        if r >= 0.6:
            return subsample_cols[0], (f"doc names {misnamed[0]}, which is not in the file (doc typo); "
                                       f"using the file's subsample weight {subsample_cols[0]}"), quote
        return subsample_cols[0], (f"CONFLICT: doc names {misnamed[0]} but the file carries {subsample_cols[0]}; "
                                   f"file weight used"), quote
    if file_specific and len(subsample_cols) == 1:
        return subsample_cols[0], "doc says to use the special sample weight in this file", quote
    if named:
        w = named[0]
        if w.startswith("WTMEC") and code == "P":
            w = "WTMECPRP"
        return w, f"doc names {w}; taken from DEMO or another file of the cycle", quote
    if re.search(r"phlebotomy weight", body, re.I) and code == "L":
        return "WTPH2YR", "doc says to use the phlebotomy weight", quote
    if code == "L" and blood and not subsample_cols:
        return "WTPH2YR", ("doc names no weight; CDC weighting tutorial: 2021-2023 blood analytes use the "
                           "phlebotomy weight WTPH2YR"), quote
    if re.search(r"(MEC exam(ination)? (sample )?weight|examined sample weight|MEC weight)", body, re.I):
        return mec, "doc says to use MEC exam weights", quote
    if len(subsample_cols) == 1:
        return subsample_cols[0], "only weight column in the file; doc names none", quote
    return mec, "default: doc names no weight; MEC exam weight", quote


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
            m = read_xpt(f, metadataonly=True)[1]
            cols = [c for c in m.column_names if c.startswith("WT")]
            doc = f.with_suffix(".htm")
            body = doc_body(doc) if doc.exists() else ""
            blood = any(c.startswith("LB") for c in m.column_names)
            w, basis, quote = choose(code, cols, body, blood)
            rows.append(dict(cycle=label, code=code, file=f.stem, stem=stem, weight_columns=";".join(cols),
                             weight=w, basis=basis, quote=quote))
    out = pd.DataFrame(rows)
    BUILD.mkdir(exist_ok=True)
    out.to_csv(BUILD / "file_weights.csv", index=False)
    return out


if __name__ == "__main__":
    o = main()
    print(o.basis.str.replace(r"WT\w+", "WT*", regex=True).value_counts())
