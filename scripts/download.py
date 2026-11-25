"""Download NHANES laboratory files, DEMO and the covariate files used for the
reference subset, with their .htm documentation, into _raw/nhanes/<cycle>/.

Usage:  uv run python -I scripts/download.py [--raw <folder>]   (default: ../_raw/nhanes next to the repository)

The list of laboratory files is read from the NHANES "Data, Documentation, Codebooks"
pages (Component=Laboratory) for each cycle. Every .xpt and its .htm doc is fetched.
A manifest (manifest.csv) records URL, local path, bytes, SHA-256 and fetch time.
Files already present are not fetched again.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

BASE = "https://wwwn.cdc.gov"
UA = {"User-Agent": "Mozilla/5.0 (compatible; nhanes-lab-percentiles build script)"}

# cycle code -> (label, page query, data folder year, file naming)
CYCLES = {
    "G": ("2011-2012", "CycleBeginYear=2011", "2011"),
    "H": ("2013-2014", "CycleBeginYear=2013", "2013"),
    "I": ("2015-2016", "CycleBeginYear=2015", "2015"),
    "J": ("2017-2018", "CycleBeginYear=2017", "2017"),
    "P": ("2017-2020P", "Cycle=2017-2020", "2017"),
    "L": ("2021-2023", "CycleBeginYear=2021", "2021"),
}

# Non-laboratory files needed for demographics, weights and the reference subset.
COVARIATES = ["DEMO", "BMX", "DIQ", "KIQ_U", "MCQ", "SMQ", "RHQ", "BPQ"]


def fname(cycle: str, stem: str) -> str:
    return f"P_{stem}" if cycle == "P" else f"{stem}_{cycle}"


def get(url: str, tries: int = 4) -> requests.Response | None:
    for i in range(tries):
        try:
            r = requests.get(url, headers=UA, timeout=120)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r
        except requests.RequestException as e:  # noqa: PERF203
            print(f"  retry {i+1} {url}: {e}", file=sys.stderr)
            time.sleep(3 * (i + 1))
    return None


def lab_links(page_html: str) -> list[str]:
    links = re.findall(r'href="(/Nchs/Data/Nhanes/Public/[^"]+\.xpt)"', page_html, re.I)
    seen, out = set(), []
    for l in links:
        if l not in seen:
            seen.add(l)
            out.append(l)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=str(Path(__file__).resolve().parents[1].parent / "_raw" / "nhanes"))
    ap.add_argument("--cycles", default="".join(CYCLES))
    args = ap.parse_args()
    raw = Path(args.raw).expanduser()
    (raw / "_pages").mkdir(parents=True, exist_ok=True)
    manifest_path = raw / "manifest.csv"
    rows = []
    if manifest_path.exists():
        rows = list(csv.DictReader(manifest_path.open()))
    have = {r["url"] for r in rows}

    def fetch(url: str, dest: Path, kind: str, cycle: str) -> None:
        if dest.exists() and dest.stat().st_size > 0:
            if url not in have:
                rows.append(dict(url=url, path=str(dest.relative_to(raw)), kind=kind, cycle=cycle,
                                 bytes=dest.stat().st_size,
                                 sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                                 fetched=datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc).isoformat()))
                have.add(url)
            return
        r = get(url)
        if r is None:
            print(f"  missing {url}")
            return
        body = r.content
        if kind == "xpt" and body[:6] != b"HEADER":
            print(f"  not an XPT file: {url}")
            return
        dest.write_bytes(body)
        rows.append(dict(url=url, path=str(dest.relative_to(raw)), kind=kind, cycle=cycle, bytes=len(body),
                         sha256=hashlib.sha256(body).hexdigest(),
                         fetched=datetime.now(timezone.utc).isoformat()))
        have.add(url)
        time.sleep(0.3)

    for code in args.cycles:
        label, query, year = CYCLES[code]
        out = raw / label
        out.mkdir(exist_ok=True)
        page = raw / "_pages" / f"lab_{label}.html"
        if not page.exists():
            r = get(f"{BASE}/nchs/nhanes/search/datapage.aspx?Component=Laboratory&{query}")
            page.write_text(r.text)
        links = lab_links(page.read_text())
        print(f"{label}: {len(links)} laboratory files")
        for link in links:
            stem = link.rsplit("/", 1)[1]
            fetch(BASE + link, out / stem, "xpt", label)
            doc = link[:-4] + ".htm"
            fetch(BASE + doc, out / (stem[:-4] + ".htm"), "doc", label)
        for cov in COVARIATES:
            f = fname(code, cov)
            fetch(f"{BASE}/Nchs/Data/Nhanes/Public/{year}/DataFiles/{f}.xpt", out / f"{f}.xpt", "xpt", label)
            fetch(f"{BASE}/Nchs/Data/Nhanes/Public/{year}/DataFiles/{f}.htm", out / f"{f}.htm", "doc", label)
        with manifest_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["url", "path", "kind", "cycle", "bytes", "sha256", "fetched"])
            w.writeheader()
            w.writerows(rows)
    print(f"manifest: {manifest_path} ({len(rows)} files)")


if __name__ == "__main__":
    main()
