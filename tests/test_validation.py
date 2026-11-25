"""Reproduce published CDC/NCHS numbers (see tests/published_values.csv and
tests/sources.csv). Tolerances are fixed in tests/validation.py:tolerance()."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validation  # noqa: E402

# Every test here recomputes published numbers from the raw NHANES files.
pytestmark = pytest.mark.skipif(
    not (validation.RAW / "2011-2012").is_dir(),
    reason=f"raw NHANES files not found in {validation.RAW}; run scripts/download.py to enable validation")


@pytest.fixture(scope="module")
def results():
    return validation.main()


def test_sample_sizes_match_exactly(results):
    bad = results[~results.n_match]
    assert bad.empty, bad[["target_id", "published_n", "n"]].to_string()


def test_values_within_tolerance(results):
    bad = results[~results.value_pass]
    assert bad.empty, bad[["target_id", "published_value", "estimate", "published_lcl", "lcl",
                           "published_ucl", "ucl", "published_se", "se"]].to_string()


def test_each_source_covered(results):
    assert set(results.source_id) == {"cdc_lead_2011_2018", "db515", "db516"}


@pytest.fixture(scope="module")
def nutrition():
    return validation.nutrition_targets()


def test_nutrition_report_pooled_sample_sizes(nutrition):
    """Pooled 2011-2014 / 2015-2018 samples match CDC Nutrition Report table n exactly."""
    bad = nutrition[~nutrition.n_match]
    assert bad.empty, bad[["target_id", "published_n", "n"]].to_string()


def test_nutrition_report_pooled_values_inside_published_ci(nutrition):
    bad = nutrition[~nutrition.value_pass]
    assert bad.empty, bad[["target_id", "published_value", "estimate", "published_lcl", "published_ucl"]].to_string()


def test_nutrition_report_covers_requested_biomarkers(nutrition):
    assert {v.split("_")[1] for v in nutrition.target_id} >= {"LBXVIDMS", "LBXFER", "LBDRFOSI", "LBDFOTSI", "LBDB12SI"}
