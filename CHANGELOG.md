# Changelog

## 1.0.0 (2026-11-25)

First release.

- 479 analytes: 110 curated clinical analytes, with pre-pandemic cycles 2011-2012 to 2017-March 2020 pooled within sex and single year of age after CDC method-bridging equations, plus 369 analytes taken from a single cycle without method review. August 2021-August 2023 is reported separately.
- 55,395 sex/age rows and 21,903 race/ethnicity rows, for the `everyone` and `reference` populations and the sensitivity population `reference_bmi_lt30`.
- Survey weights chosen per file from the CDC documentation (`scripts/weights.py`); curated analytes' weights pinned by a test.
- Suppression and flags for small cells, results below the limit of detection, cells covering unsampled groups (`population_not_sampled`) and NCHS presentation standards.
- Validation against 320 published CDC/NCHS numbers, including 258 pooled-cycle values from the CDC Nutrition Report.
- Web JSON for lifefrom.ai.
