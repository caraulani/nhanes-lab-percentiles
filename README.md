---
license: cc-by-4.0
pretty_name: NHANES Lab Percentiles
language:
- en
tags:
- nhanes
- laboratory
- biomarkers
- percentiles
- survey-statistics
- public-health
size_categories:
- 10K<n<100K
configs:
- config_name: percentiles_sex_age
  data_files: data/percentiles_sex_age.csv
- config_name: percentiles_race_ethnicity
  data_files: data/percentiles_race_ethnicity.csv
- config_name: analytes
  data_files: data/analytes.csv
- config_name: coverage
  data_files: data/coverage.csv
---

# NHANES Lab Percentiles

Interactive version: https://lifefrom.ai/research/open-data/population-percentiles

Survey-weighted percentiles for 479 blood and urine laboratory analytes measured in the US National Health and Nutrition Examination Survey (NHANES). Results are given by sex and age band, and for adults also by race and Hispanic origin. The tables describe the US civilian noninstitutionalized population.

Version 1.0.0, released 2026-11-25, by LifeFrom (https://lifefrom.ai).

**Source data: National Center for Health Statistics (NCHS), Centers for Disease Control and Prevention. National Health and Nutrition Examination Survey, public-use laboratory, demographic, examination and questionnaire files, 2011-2012 to August 2021-August 2023.** NHANES data are in the public domain. NCHS and CDC are not involved in this release.

## What it is for

The tables answer one question: "Where does a value sit in the US population of the same sex and age?" For example:

- a researcher can place a cohort's measurements against a national distribution;
- an educator can show how wide normal variation is;
- an app can show a percentile next to a result.

Each value comes with its sample size, sum of weights, share below the limit of detection, and flags. For the 2.5th, 50th and 97.5th percentiles, design-based 95% confidence limits are also given.

## Populations

- **`everyone`**: all examined participants with a result, including people with diabetes, kidney disease, smokers and pregnant women.
- **`reference`**: examined participants without the common conditions listed under "Reference subset" below (definition in section 7 of METHODS.md).
- **`reference_bmi_lt30`**: the reference subset further restricted to BMI below 30 (sensitivity population).

A value's percentile depends on the laboratory method. Results from another lab or assay may not be comparable with NHANES results. Method notes are in `analytes.csv`.

## Files

| File | Rows | Content |
|---|---|---|
| `data/percentiles_sex_age.csv` | 55,395 | Long format: one row per analyte × cycle set × population × sex × age band. Populations: `everyone`, `reference`, and the sensitivity population `reference_bmi_lt30` |
| `data/percentiles_race_ethnicity.csv` | 21,903 | Adults 18+: one row per analyte × cycle set × population (`everyone`, `reference`) × sex × race/Hispanic origin. For research description only; see the note below |
| `data/analytes.csv` | 479 | Catalog: CDC files, variables, units, specimen, cycles used, weight per cycle, LOD status per cycle, method adjustments and notes, doc URLs |
| `data/coverage.csv` | 1,503 | Every analyte × cycle in the raw data, with n, and whether and why it was used |
| `data/reference_exclusions.csv` | 72 | Per cycle and criterion: examined participants, the number assessable, and the number flagged |
| `web/index.json`, `web/markers/<analyte_id>.json` | 479 | Compact JSON for the lifefrom.ai website |
| `config/analyte_rules.csv` | 110 | Hand-set pooling and method-adjustment rules for clinical analytes, with the CDC doc each rule comes from |
| `tests/published_values.csv`, `tests/sources.csv`, `tests/nutrition_report_tables.txt` | 62 + 258 | Published CDC/NCHS numbers we reproduce |
| `scripts/` | | `download.py` fetches the raw files; `build.py` rebuilds everything |

### Columns of `percentiles_sex_age.csv`

| Column | Meaning |
|---|---|
| `analyte_id` | Stable kebab-case id. Equal to `marker_slug` when the analyte is the main match for a LifeFrom marker |
| `nhanes_variable` | NHANES variable name(s) used, `;`-separated if the name changed across cycles |
| `marker_slug` | LifeFrom marker id (`markers.csv` of LifeFrom release 1); empty if there is no match |
| `analyte_name` | CDC's SAS label without the unit |
| `specimen` | serum, plasma, whole blood, urine, or a broader term where not curated |
| `unit` | Unit as in the CDC label (the conventional unit when CDC gives two) |
| `cycle_set` | `2011-2020P` (pooled pre-pandemic cycles) or `2021-2023` (August 2021-August 2023, never pooled) |
| `population` | `everyone` or `reference` |
| `sex` | `all`, `female`, `male` |
| `age_band` | `1-5`, `6-11`, `12-17`, `18-29` … `70-79`, `80+` (NHANES top-codes age at 80), `18+` |
| `measured_population` | The part of the labelled sex × age group that was actually measured in at least one cycle, for example `female 18-49` |
| `coverage_pct` | Share of the labelled group's population (years-averaged MEC weights) that lies in measured sex × single-year-of-age units. Below 95%, the cell has no percentiles |
| `cycles` | Cycles that contributed participants to this cell |
| `n` | Unweighted number of participants |
| `sum_weights` | Sum of the (pooled) survey weights: the estimated population size of the cell |
| `n_below_lod`, `pct_below_lod` | Unweighted count and weighted % of results below the limit of detection, as coded by CDC. Empty when any cycle in the block lacks a CDC comment variable; the flags say why |
| `mean` | Weighted mean (CDC fill values used below the LOD) |
| `p1` … `p99` | Weighted percentiles 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99 |
| `p2_5_lcl`, `p2_5_ucl`, `p50_lcl`, `p50_ucl`, `p97_5_lcl`, `p97_5_ucl` | Woodruff 95% confidence limits |
| `df` | Design degrees of freedom of the cell (PSUs minus strata containing its participants) |
| `n_eff_median` | Effective sample size at the median (n / design effect) |
| `flags` | `;`-separated. See the flag table below |

`percentiles_race_ethnicity.csv` has the same columns. In place of `age_band` it has `race_ethnicity`, which uses the NHANES `RIDRETH3` groups. All its rows are ages 18+.

An empty value means the value was suppressed, not that it is zero. The flags say why.

### Flags

These are the only flag values that occur in the data. `pX` stands for a percentile column such as `p2_5`.

| Flag | Meaning |
|---|---|
| `n<30_suppressed` | n < 30: only n and sum of weights published |
| `median_only_n<120` | 30 ≤ n < 120: median and its CI only |
| `p1_p99_withheld_n<300` | n < 300: 1st and 99th percentiles withheld |
| `pX<LOD` | That percentile is at or below the weighted share of results below the LOD, so it is withheld |
| `over_50pct_below_lod` | More than half the results are below the LOD: the mean is also withheld |
| `population_not_sampled` | Less than 95% of the labelled group was measured in any cycle: percentiles and mean withheld. See `measured_population` |
| `lod_not_coded` | The CDC doc states an LLOD, but the file has no comment variable, so below-LOD results cannot be identified |
| `lod_not_documented` | No comment variable and no LLOD in the CDC doc |
| `lod_coded_in_some_cycles_only` | A comment variable exists in some pooled cycles only, so `pct_below_lod` is not computed |
| `df<8` | Design degrees of freedom of the cell below 8 (NCHS presentation standard) |
| `n_eff<30` | Effective sample size at the median below 30 (NCHS presentation standard) |
| `pX_ci_clipped` | The Woodruff probability limit fell outside [0, 1] and was clipped |

### Race and Hispanic origin table

`percentiles_race_ethnicity.csv` describes the distribution of values across the self-reported NHANES race and Hispanic origin groups (`RIDRETH3`), which are social categories; differences between groups can reflect diet, sun exposure, smoking, income, access to care and other exposures.

### Web JSON

`web/index.json`:

```
{version, generated, source, cycles, markers:[{slug, analyte_id, name, unit, specimen, n_total, has_reference_subset}]}
```

`web/markers/<analyte_id>.json`:

```
{analyte_id, marker_slug, name, unit, alt_unit?, alt_factor?, specimen, cycles, method_note,
 populations:{everyone:{all|female|male:{"<age band>":{n, p:{"1":…,"2.5":…,…,"99":…}, mean, pct_below_lod, flags, measured_population}}}, reference:{…}},
 cycles_2021_2023?:{cycles, nhanes_variable, populations}}
```

`populations` holds the 2011-2020P block when it exists. Otherwise it holds 2021-2023. `alt_factor` converts to the second unit CDC provides: value_alt = alt_factor × value.

## How it was built (summary)

Full detail is in METHODS.md.

- **Cycles.**
  - The `2011-2020P` block pools 2011-2012, 2013-2014, 2015-2016 and 2017-March 2020 pre-pandemic. 2017-2018 is used only for analytes missing from the 2017-2020 files.
  - Weights are scaled by years, as in the CDC weighting tutorial, but within each sex × single year of age. Each group is averaged only over the cycles that measured it. This keeps groups measured in only some cycles (for example ferritin in older women) from being under-weighted, and makes `sum_weights` estimate the cell's population.
  - August 2021-August 2023 is reported separately, because CDC advises against combining it with earlier cycles.
- **Method changes.** For 110 clinical analytes we read every cycle's CDC documentation. There are three outcomes:
  - CDC published a bridging equation: we applied it to put older cycles on the current method, then pooled.
  - CDC published no usable equation: we used only the most recent consistent span.
  - No change is documented: we pooled all cycles.

  `analytes.csv` records the rule and the doc it came from. For the other 369 analytes, mostly environmental chemicals, we did not review method history. Their pre-pandemic block uses only the most recent cycle.
- **Weights.** The weight for each file is the one its CDC doc names: fasting, special subsample, blood-metal, surplus-specimen or MEC exam weight. For 2021-2023 blood analytes it is the phlebotomy weight. `scripts/weights.py` reads the docs. For the 110 curated analytes we reviewed every weight and a test pins it.
- **Percentiles.** The smallest value whose weighted cumulative share reaches p, with no interpolation. This definition reproduced all 20 CDC blood-lead percentiles we tested exactly.
- **Confidence limits.** Woodruff limits use Taylor-linearised variance with NHANES strata and PSUs. Subgroups are domains, not subsets.
- **Suppression.**
  - n < 30: nothing published.
  - 30 ≤ n < 120: median only.
  - n < 300: no 1st or 99th percentile.
  - No percentile at or below the share of results below the limit of detection.
  - No percentiles for a cell whose label covers groups that were never measured (`population_not_sampled`).
  - Flags for design df < 8 and effective n < 30 (NCHS presentation standards for proportions, applied to the median).
- **Reference subset.** Every default criterion is measured the same way in every cycle:
  - pregnancy;
  - diabetes (self-report or HbA1c ≥ 6.5%);
  - kidney markers (self-report, dialysis, ACR ≥ 30 mg/g, or race-free CKD-EPI 2021 eGFR < 60 on bridged creatinine; thresholds from KDIGO 2024);
  - serum cotinine > 3 ng/mL (Benowitz et al. 2009);
  - anemia treatment;
  - heart disease or stroke;
  - thyroid, liver or cancer history.

  CRP is not used because it was not measured in 2011-2014. Obesity is a separate sensitivity population (`reference_bmi_lt30`), not part of the default subset.

## Validation

`tests/validation.py` recomputes 320 published CDC/NCHS numbers from the raw files with this release's code:

- **Blood lead.** 20 percentiles and 3 geometric means from CDC's 2011-2018 blood lead table. All 20 percentiles match exactly. All 3 geometric means match to the published digits. All 40 confidence limits fall within 6%.
- **Cholesterol.** 24 prevalences from NCHS Data Brief 515 (2021-2023). The largest difference is 0.05 percentage points.
- **Diabetes.** 15 prevalences from NCHS Data Brief 516 (2021-2023). The largest difference is 0.04 percentage points.
- **Pooled cycles: CDC Nutrition Report reference-interval tables.** These are the 2026 online tables. We checked 258 numbers (means and the 2.5th, 5th, 50th, 95th and 97.5th percentiles), each recomputed on the same pooled cycles:
  - vitamin D 2011-2014 and 2015-2018;
  - ferritin, females 12-49, 2015-2018;
  - RBC folate, serum folate and vitamin B12, 2011-2014.

  All 258 estimates fall inside CDC's 95% CI. The median relative difference is 0.13% and the largest is 2.5%.

All 320 sample sizes match exactly. Consistency tests check that:

- percentiles never decrease within a row;
- every published cell meets the sample-size rules;
- no percentile is published at or below the share of results below the LOD;
- MEC weights sum to CDC's population control totals within 0.5%;
- the curated analytes' weights match the reviewed configuration;
- age bands add up to the adult cell;
- the pooled ferritin female 18+ median (65.3 ng/mL) is within 2% of an independent estimate (about 65.7);
- `pct_below_lod` is never filled in when the LOD coding is unknown.

Run `uv run python -I -m pytest -q tests`. Tests that read only the shipped files always run. The validation tests and four checks in `test_outputs.py` recompute values from the raw NHANES files; they are skipped, with a message, when those files are not present.

## Getting the raw data and rebuilding

The raw NHANES files are not redistributed here. `scripts/download.py` fetches them from the CDC website:

1. For each cycle it reads the NHANES "Data, Documentation, Codebooks" page for the Laboratory component (`https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&...`).
2. It downloads every laboratory `.xpt` file listed there with its `.htm` documentation, plus the DEMO, BMX, DIQ, KIQ_U, MCQ, SMQ, RHQ and BPQ files used for the reference subset.
3. Files go to `../_raw/nhanes/<cycle>/`, a folder next to the repository folder (change it with `--raw` or the `NHANES_RAW` environment variable). Files already present are skipped. `manifest.csv` records each file's URL, size, SHA-256 and fetch time.

The CDC population control totals used by one test are saved by hand from https://wwwn.cdc.gov/nchs/nhanes/ResponseRates.aspx into `../_raw/nhanes/_poptotals/`.

```
uv sync
uv run python -I scripts/download.py     # about 300 MB
uv run python -I scripts/build.py        # rebuilds data/ and web/
uv run python -I -m pytest -q tests
```

`build.py` also reads the LifeFrom release 1 marker list (`../_raw/lifefrom-release-1/markers.csv`, or `LIFEFROM_RAW`) to fill `marker_slug`.

## Limitations

- **Pooled years.** Pooled estimates describe 2011-2020 as a whole, not any single year. Trends, such as falling blood lead, are averaged.
- **Groups differ by cycle.** Some cycles measured an analyte only in some sexes or ages. For example, blood lead for ages 1-5 is not in the public 2017-2020 file, and ferritin in 2015-2016 covers females 12-49 and children only. Weights are pooled within sex × age, so each group is averaged over the cycles that measured it. That means different groups in one cell can represent different years. The `cycles` column shows which cycles fed each cell.
- **Method changes in the default tier.** For the 369 analytes we did not review, method changes between cycles are possible. Using a single cycle avoids pooling across them, but values may still differ from other years.
- **Applied CDC equations add error.** They are estimates from bridging studies of a few hundred samples.
- **Below the LOD.** CDC fill values are LOD/√2. Percentiles at or below the share below the LOD are not published. The mean uses the fill values. For analytes without a CDC comment variable (many clinical chemistry tests), the share below the LOD is unknown and left empty. Their lowest percentiles may include values at or below the assay's LLOD.
- **The reference subset relies on self-report and single measurements.** Pregnancy status is released only for females 20-44. Conditions in the medical-conditions questionnaire are asked only from age 20. A single elevated ACR or low eGFR is not a CKD diagnosis, which KDIGO defines by persistence over 3 months. For analytes used to define the subset (HbA1c, cotinine, urine albumin/creatinine ratio, creatinine), the reference distribution is cut off by design.
- **"All sexes" cells** mix two distributions for sex-dependent analytes such as testosterone or ferritin. Use the female and male cells.
- **Race and Hispanic origin** are self-reported census-style categories.

## How to cite

LifeFrom; Caraulani, J. (2026). NHANES Lab Percentiles: survey-weighted percentiles of US laboratory analytes by sex and age (Version 1.0.0) [Data set]. https://github.com/caraulani/nhanes-lab-percentiles

Please also cite the source: National Center for Health Statistics. National Health and Nutrition Examination Survey Data. Hyattsville, MD: U.S. Department of Health and Human Services, Centers for Disease Control and Prevention. https://www.cdc.gov/nchs/nhanes/

See `CITATION.cff`.

## Licence

- **Data** (`data/`, `web/`, `config/`, `tests/*.csv`): CC BY 4.0.
- **Code** (`scripts/`, `tests/*.py`): MIT.

The underlying NHANES data, produced by NCHS/CDC, are in the public domain (U.S. Government work). See `LICENSE`.
