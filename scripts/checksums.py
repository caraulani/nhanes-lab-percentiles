"""Write SHA256SUMS for every released file (data/, web/, config/, scripts/, tests/, docs)."""
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = ["data/*.csv", "web/index.json", "web/markers/*.json", "config/*.csv", "scripts/*.py", "tests/*.py",
            "tests/*.csv", "tests/*.txt", ".gitignore", "README.md", "METHODS.md", "CITATION.cff", ".zenodo.json", "LICENSE",
            "CHANGELOG.md", "pyproject.toml", "uv.lock"]


def main():
    files = sorted({p for pat in PATTERNS for p in ROOT.glob(pat) if p.is_file()})
    lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT).as_posix()}" for p in files]
    (ROOT / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    print(f"SHA256SUMS: {len(lines)} files")


if __name__ == "__main__":
    main()
