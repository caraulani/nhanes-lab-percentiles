# Methods

This file gives enough detail to rebuild every number in the release. All code is in `scripts/`. Every source we read was accessed on 2026-11-25.

## 1. Source data

Source: the National Health and Nutrition Examination Survey (NHANES), run by the National Center for Health Statistics (NCHS), part of the US Centers for Disease Control and Prevention (CDC). We used the public-use files only.

`scripts/download.py` reads the NHANES "Data, Documentation, Codebooks" page for the Laboratory component of each cycle at `https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&...`. For each cycle it then downloads:

- every laboratory `.xpt` file listed on that page;
- the matching `.htm` documentation for each of those files;
- DEMO, BMX, DIQ, KIQ_U, MCQ, SMQ, RHQ and BPQ, with their docs.

| Code | Cycle | Weight years | Notes |
|---|---|---|---|
| G | 2011-2012 | 2.0 | |
| H | 2013-2014 | 2.0 | |
| I | 2015-2016 | 2.0 | |
| J | 2017-2018 | 2.0 | Used only when an analyte is not in the 2017-March 2020 files |
| P | 2017-March 2020 pre-pandemic (`P_` files) | 3.2 | Contains the 2017-2018 participants |
| L | August 2021-August 2023 | 2.0 | Kept separate (see 3.2) |

In total the script downloaded 446 `.xpt` files and 446 `.htm` docs (892 files). `../_raw/nhanes/manifest.csv` (the raw folder sits next to the repository folder) records the URL, SHA-256 and fetch time of each one. The raw files are not part of the release.

Files are read with `pyreadstat`, with SAS labels decoded as Latin-1. A few NHANES labels contain non-UTF-8 bytes.

## 2. Which variables count as analytes

`scripts/catalog.py` scans every laboratory file. It keeps a numeric variable when all of the following hold:

1. Its SAS label ends in a unit, such as `(mg/dL)`. The unit must look like a measurement unit: it contains `/`, `%`, `fL`, `pg`, `IU`, `U/L`, `mOsm`, "ratio" or "units".
2. It is not a weight, comment code, sequence number, status or timing variable.
3. It has at least 10 distinct values. This drops categorical serology results. The rule is waived for the curated analytes in `config/analyte_rules.csv`. Basophil and eosinophil counts are reported to 0.1 × 1000 cells/µL, so they have few distinct values.

`scripts/build.py` then removes the following:

- **Pooled-sample files** (BFRPOL, DOXPOL, PCBPOL, PSTPOL). These have no participant ID, so person-level percentiles cannot be computed from them.
- **Smokers' special samples.** These files carry their own smoking weights (WTFSM, WTALDS2Y, WTSVS2YR) and are outside the general population design.
- **Surplus-specimen files with no weight variable.**
- **SSTST (2017-2018).** CDC folded these results into `P_TST`.
- **UCFLOW.** Urine flow rate describes the collection, not an analyte concentration.
- **A small number of single variables**, each logged with its reason.

`build/excluded_items.csv` (written by the build, not shipped) lists every exclusion and its reason.

**Unit pairs.** CDC often publishes one result twice: once in a conventional unit and once in SI. Two variables in the same file are treated as one analyte when their labels match apart from the unit, or when their names follow the `X`/`XSI` or `LBX…`/`LBD…SI` pattern. The conventional variable is the primary one. The release reports percentiles in the primary unit.

**Second unit.** `alt_unit` and `alt_factor` are given only when both variables are in the CDC file. `alt_factor` is the median ratio of the two variables. It is published only when the 5th to 95th percentile spread of that ratio is 2% or less of the median, which shows the relation is linear. We never convert units from outside sources.

**Analyte key.** An analyte is keyed by its NHANES variable name. Four exceptions:

- Hydroxycotinine is `LBXHCT` in the COT files of 2013-2018. In the CBC files, `LBXHCT` means hematocrit. The cotinine variable is therefore keyed as `LBXHCOT`, the name CDC uses from 2021.
- The 2021-2023 triglyceride variable `LBXTLG` is shown under the id `triglycerides`. The CSV keeps `LBXTLG` in `nhanes_variable`, and the method note says the method changed.
- Vitamin B12 is `LBXB12` in 2011-2012 and `LBDB12` in 2013-2014. Both are keyed as `LBDB12` (see 3.2).
- When one variable appears in several files of the same cycle, the file with the most results is used. For example, urine creatinine is taken from ALB_CR.

## 3. Cycles, pooling and method changes

### 3.1 Two cycle sets

Each analyte gets up to two blocks of estimates:

- **`2011-2020P`**: pooled from G, H, I and P. J replaces P only when the analyte is not in the P files.
- **`2021-2023`**: L on its own.

The CDC weighting tutorial says: "It is generally not recommended to combine the August 2021-August 2023 cycle with other cycles given the 1.5-year gap between this cycle and the 2017-March 2020 cycle." (https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx). We follow this, so 2021-2023 is never pooled with earlier cycles.

### 3.2 Curated analytes (`config/analyte_rules.csv`, 110 variables)

These are the analytes with routine clinical use. For each one we read the CDC documentation of every cycle. `scripts/docscan.py` pulls out every sentence that mentions a method, instrument or laboratory change, a crossover or bridging study, or a regression adjustment, and writes them to `build/doc_notes.csv`. We read those sentences and set a rule for each analyte:

- **No documented method change:** pool all available cycles.
- **Documented change with a CDC-published forward equation:** apply the equation to the earlier cycles so that they are on the newest method's scale, then pool. The equation, the cycles it was applied to and the source doc are in `analytes.csv` under `method_adjustment` and `method_source`. Equations used:
  - **Biochemistry profile** (2017 move from Beckman DxC to Roche Cobas 6000; `BIOPRO_J` bridging table): applied to 2011-2016 for albumin, ALP (log-Deming), AST, ALT, BUN, cholesterol, CK, creatinine, GGT, iron, LDH, triglycerides and uric acid. CDC marked "adjustment not recommended" for bicarbonate, calcium, chloride, glucose, potassium, sodium, phosphorus, total bilirubin and total protein, so those are pooled unadjusted.
  - **Fasting glucose and 2-hour OGTT glucose** (2015 move from Roche C501 to C311; `GLU_I`, `OGTT_I`): applied to 2011-2014. CDC had already adjusted the 2015-2016 C501 results.
  - **Insulin** (2013 move from Roche to Tosoh; `INS_H`): the forward Deming equation on the log10 scale, applied to 2011-2012.
  - **Apolipoprotein B** (2015 move from Siemens ProSpec to Roche Cobas; `APOB_I`): applied to 2011-2014.
- **Documented change with no usable published adjustment:** use only the most recent consistent span. The reason is in `method_note`. Examples:
  - **CBC:** CDC states the 2012 change from Coulter HMX to DxH could not be bridged, so CBC uses 2013-2020P.
  - **hs-CRP:** the 2017 equation is valid only up to 23 mg/L, so hs-CRP uses 2017-2020P only.
  - **Testosterone:** the 2013 bridging equation's units are ambiguous in the doc, so testosterone uses 2013-2016.
  - **Globulin and calculated osmolality:** derived from components that need adjustment, so they use 2017-2020P only.
  - **Estradiol and SHBG:** comparability of the 2017-2020 panel with earlier cycles is not documented, so they use 2017-2020P only.
- **Vitamin B12** pools 2011-2012 (`LBXB12`) with 2013-2014 (`LBDB12`). CDC re-derived the 2013-2014 values under a new name with a QC-based Deming adjustment (VITB12_H doc). Both cycles used the Roche E170 assay, and the CDC Nutrition Report pools the two cycles in its 2011-2014 reference-interval table. Our pooled 2011-2014 B12 percentiles reproduce that table (section 8).
- **Target group change:** in the P file, transferrin receptor was measured in 2019-2020 only in ages 1-5 and in females 12-49. Transferrin receptor therefore uses 2015-2016 and 2017-2018.

### 3.3 Other analytes (default rule, 369 analytes)

We did not review the method history of the environmental chemicals, fatty acids, folate forms and other analytes. For each of these, the `2011-2020P` block uses only the most recent pre-pandemic cycle that has the analyte, in the order P, J, I, H, G. No pooling is done, so no method change can be pooled across. This is stated in `method_note` and `tier = default`.

### 3.4 Pooled weights, by domain

The CDC tutorial scales each cycle's weight by its share of the years pooled. Its example for 2015-2016 plus 2017-March 2020 is `WTMEC52YR = 3.2/5.2 * WTMECPRP`, with two-year cycles getting `2/5.2`. That factor is only right when every cycle measured the same population.

Some analytes were measured in different groups in different cycles. Ferritin, for example, was measured in females 12-49 and children in 2015-2016, but in everyone 12 and older in 2017-March 2020. With one factor per cycle, women 18-49 were weighted by both cycles (2/5.2 + 3.2/5.2), and women 50 and older by one (3.2/5.2). That under-represented older women. It also understated `sum_weights` wherever a group was measured in only some cycles. Version 1.0.0 drafts had this error (`scripts/pooling.py` documents it).

We now scale within units of sex × single year of age:

`w*_i = w_i × years(cycle_i) / (sum of years over the cycles that measured anyone in unit_i)`

Cycle years are 2.0 for G, H, I and J, and 3.2 for P. When every cycle covers every unit, this reduces to the CDC factor. Within each unit, the scaled weights sum to the years-averaged population of that unit over the cycles that measured it. So every cell's `sum_weights` estimates the population of the cell, and no group is over-weighted.

**Before and after for ferritin**, everyone, 2011-2020P (both 2015-2016 and 2017-March 2020):

| Cell | n | Median before | Median after (95% CI) | Sum of weights before | after |
|---|---|---|---|---|---|
| Females 18+ | 5,805 | 59.1 | 65.3 (62.4-68.6) | 98.9 million | 121.3 million |
| All 18+ | 9,863 | 89.1 | 103 (98.3-107) | 168.7 million | 234.7 million |
| Males 18+ | 4,058 | 156 | 156 | 69.8 million | 113.4 million |

The male 18+ median is unchanged because males were measured in only one cycle. Its weights were simply scaled too low before.

An independent calculation put the corrected female 18+ median at about 65.7. A test pins that case (`test_ferritin_female_adult_pooled_median`). Across all `everyone` cells, 27 medians moved by more than 0.5%. They belong to vitamin B12 (16, now pooled over two cycles), ferritin (6), soluble transferrin receptor (4) and methylmalonic acid (1).

**Cells whose label covers unsampled groups.** For each cell we compute `coverage_pct`. This is the share of the labelled population (sex filter × age range), weighted by years-averaged MEC exam weights, that falls in sex × age units where at least one cycle measured the analyte. `measured_population` describes the measured part, for example `female 18-49`.

If coverage is below 95%, the cell's percentiles and mean are withheld and the cell is flagged `population_not_sampled`. Example: in 2021-2023, ferritin was measured only in children 1-5 and in females 12-49. Its "all sexes, 18+" cell therefore says `measured_population = female 18-49` and `coverage_pct = 26.8`, and has no percentiles. The female 18-29, 30-39 and 40-49 cells are published.

## 4. Weights

`scripts/weights.py` chooses the weight for every laboratory file from what the CDC doc for that file says. It reads the doc text outside the codebook and collects the sentences that mention weights. The chosen weight and its basis are written to `build/file_weights.csv`. The rules, in order:

1. **The doc names a weight variable that is in the file** (for example `WTSH2YR` in PBCD_H, `WTFOLPRP` in P_FOLATE, `WTTSTPP` in P_TST): use it. If the file also carries `WTPH2YR`, the named subsample weight wins.
2. **The doc names a weight that is not in the file, but the name is a clear typo of the file's only subsample weight.** GLU_G, GLU_H and INS_H name `WTSFA2YR` while the file has `WTSAF2YR`. TRIGLY_L names `WTSAFPRP` while the file has `WTSAF2YR`. In these cases we use the file's weight and record the typo. Name similarity must be at least 0.6.
3. **The doc says to use "the special sampling weights in this file"** (fasting, OGTT and subsample A/B/C files): use the file's subsample weight.
4. **The doc names a weight that lives elsewhere** (`WTPH2YR` for COT_L): take it from another file of the same cycle.
5. **The doc says to use MEC exam weights or examined sample weights:** use `WTMEC2YR`, or `WTMECPRP` for P.
6. **2021-2023 blood files whose doc names no weight** (BCHE_L): use `WTPH2YR`. The CDC weighting tutorial says 2021-2023 blood analytes use the phlebotomy weight.
7. **Otherwise:** use the file's only weight column if there is one, or the MEC exam weight. The basis is recorded as such.

For the 110 curated analytes, we reviewed the chosen weight for every cycle and recorded it in `config/analyte_rules.csv` (`expected_weights`). The test `test_curated_weights_match_config` fails if the build picks anything else.

Two doc/file conflicts were found, in PAHS_I and PERNTS_I. Those files belong to the smokers' special sample, which is excluded anyway.

Rows used in a cell must meet all of these conditions:

- the participant was examined (`RIDSTATR = 2`);
- the result is not missing;
- the chosen weight is positive;
- age is 1 year or more.

The weight used for each cycle is listed in `analytes.csv` (`weights`).

## 5. Variance and confidence intervals

Variance uses first-stage with-replacement Taylor linearisation, with strata `SDMVSTRA` and PSUs `SDMVPSU`, as NCHS recommends (variance estimation tutorial). Strata are keyed by cycle, so pooled cycles never share a stratum.

Cells are domains, not subsets. Out-of-domain rows contribute zero, and every PSU in every stratum of the cycles used is kept. The design comes from all examined participants in DEMO.

The variance of a linearised total is `V = sum_h n_h/(n_h - 1) * sum_j (z_hj - mean_h(z))^2`. Here `z_hj` is the sum over PSU `j` in stratum `h`.

**Degrees of freedom.** The design degrees of freedom are the number of PSUs minus the number of strata over the cycles used. Examples: 72 for G+H+I+P (140 PSUs, 68 strata), 62 for G+H+I+J, 15 for L. These df give the t multiplier. Each cell also reports the domain df: PSUs minus strata among those containing the cell's participants.

**Percentiles.** The weighted empirical CDF is `F(t) = sum(w_i : y_i <= t) / sum(w_i)`. The percentile `Q(p)` is the smallest observed value `v` with `F(v) >= p`, with no interpolation.

We first coded the SAS SURVEYMEANS interpolated definition. It missed 5 of the 20 published CDC blood-lead percentiles in the third significant digit. The step definition reproduced all 20 exactly (section 8), so it is the one used. The interpolated definition is kept in `scripts/survey.py` (`QUANTILE_METHOD`) for comparison.

**Woodruff limits** are given for the 2.5th, 50th and 97.5th percentiles in every published cell:

1. Compute `F(Q(p))` and its linearised standard error. The residuals are `w_i (I(y_i <= Q(p)) - F(Q(p))) / W`.
2. Form `p ± t(0.975, df) × SE`.
3. Clip to [0, 1]. If clipping occurs, the cell gets the flag `<pct>_ci_clipped`.
4. Map both limits back through `Q(·)`.

The limits are not symmetric.

**Mean** is the weighted arithmetic mean. It is computed with fill values for results below the limit of detection (LOD), as CDC releases them.

## 6. Cells, suppression and flags

**Groups.** Sex is all, female or male (`RIAGENDR`). Age bands use `RIDAGEYR`: 1-5, 6-11, 12-17, 18-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+, and 18+. NHANES top-codes age at 80. A band is output only if the analyte was measured in that age range.

A second table splits adults (18+) by `RIDRETH3` and sex. The categories are Mexican American, Other Hispanic, Non-Hispanic White, Non-Hispanic Black, Non-Hispanic Asian, and Other race including multiracial.

**Rules applied to every cell:**

| Rule | Effect | Flag |
|---|---|---|
| Probability limit of a Woodruff interval falls outside [0, 1] and is clipped | Published, flagged | `pX_ci_clipped` |
| n < 30 | Only n and sum of weights are published | `n<30_suppressed` |
| 30 ≤ n < 120 | Median and its CI only | `median_only_n<120` |
| n < 300 | 1st and 99th percentiles not published | `p1_p99_withheld_n<300` |
| Weighted % below LOD ≥ p | That percentile is not published; the share is the rounded value published in `pct_below_lod` | `pX<LOD` |
| Weighted % below LOD > 50% | All percentiles at or below the median, and the mean, are not published | `over_50pct_below_lod` |
| Labelled population < 95% covered by measured sex × age units | Percentiles and mean not published | `population_not_sampled` |
| A CDC comment (LOD) variable is missing in some or all cycles of the block | `pct_below_lod` left empty; LOD rules cannot be applied | `lod_not_coded`, `lod_not_documented` or `lod_coded_in_some_cycles_only` |
| Domain df < 8 | Published, flagged | `df<8` |
| Effective n at the median < 30 | Published, flagged | `n_eff<30` |

Where the rules come from:

- **n ≥ 120** is the CLSI EP28 minimum for nonparametric 2.5th and 97.5th percentiles, as set in the task specification. We did not fetch CLSI EP28 itself.
- **n ≥ 300 for p1 and p99** is our own rule, not from a published standard.
- **df < 8 and effective n < 30** come from the NCHS Data Presentation Standards for Proportions: Parker JD, Talih M, Malec DJ, et al. Vital Health Stat 2(175), 2017, text table on page 2. Those standards are written for proportions. We apply them to percentiles by treating the median as the proportion `F(Q(0.5))`. Effective n is `n / deff`, where `deff = SE^2 / (p(1 - p)/n)`.

**Limit of detection.** We use CDC's own LOD coding. Results below the LOD are released by CDC as LOD/√2 fill values. The CDC comment variable (`...LC`, value 1 = below LOD) marks them. `build.py` finds the comment variable by the CDC naming pattern, falling back to label matching.

For every analyte and cycle we also check the CDC doc for a stated lower limit of detection: an "LLOD" or "limit of detection" passage that names the variable. `analytes.csv` records the result per cycle in `lod_by_cycle`, as one of:

- `comment <variable>`;
- `LLOD in doc, no comment variable`;
- `none documented`.

It also records one status per block in `lod_status_2011_2020P` and `lod_status_2021_2023`.

`pct_below_lod` (weighted) and `n_below_lod` (unweighted) are published only when every cycle in the block has a comment variable. Otherwise both are empty and every cell carries one of these flags:

- `lod_not_coded`: the LLOD is stated in a doc, but there is no comment variable;
- `lod_not_documented`: neither a comment variable nor a stated LLOD;
- `lod_coded_in_some_cycles_only`.

The value is never set to 0 by default.

Results for the 110 curated analytes (2011-2020P block):

| Status | Analytes |
|---|---|
| Comment variable in every cycle | 29 |
| LLOD stated in the doc, no comment variable | 54 |
| None documented | 18 |
| Comment variable in some cycles only | 7 |

The remaining 2 curated analytes are measured only in 2021-2023.

## 7. Populations

**`everyone`**: all examined participants with a result, weighted.

**`reference`**: `everyone` minus participants with any of the default exclusions below.

Every default criterion uses the same variable, asked or measured with the same age universe, in every cycle (2011-2012 to 2021-2023). So the subset means the same thing in every pooled cycle. We checked availability by cycle. `data/reference_exclusions.csv` reports, per cycle and criterion:

- the number of examined participants;
- the number for whom the criterion could be assessed;
- the number flagged.

A participant is excluded only on positive evidence. A missing answer or lab value does not exclude.

| Criterion | Definition | Age universe (same in every cycle) | Source |
|---|---|---|---|
| Pregnant | `RIDEXPRG = 1` or urine pregnancy test `URXPREG = 1` | Females 20-44 (CDC releases pregnancy status only for these ages) | Fulgoni et al. 2019, Table 1 ("Pregnancy") |
| Diabetes | Told by a doctor (`DIQ010 = 1`), or HbA1c ≥ 6.5% | DIQ010: ages 1+. HbA1c: 12+ | Thresholds as in NCHS Data Brief 516. Fasting glucose is not used because it exists only for the fasting subsample. Insulin and pill items are not used because their skip pattern differs in the 2017-March 2020 and 2021-2023 files |
| Kidney | Told of weak or failing kidneys (`KIQ022 = 1`), dialysis in the past 12 months (`KIQ025 = 1`), urine albumin/creatinine ≥ 30 mg/g (ages 6+), or eGFR < 60 mL/min/1.73 m² (ages 18+) | KIQ: 20+. ACR: 6+ (2011-2014 measured ACR from age 6, later cycles from 3). Creatinine: 12+ | KDIGO 2024 (below) |
| Smoking | Serum cotinine > 3 ng/mL | Ages 12+ | Benowitz et al. 2009 (below) |
| Anemia treatment | `MCQ053 = 1` (treatment for anemia in the past 3 months) | Ages 1+ | Fulgoni et al. 2019, Table 1 ("Anemia") |
| Heart disease or stroke | `MCQ160B` to `MCQ160F` = 1 (heart failure, coronary heart disease, angina, heart attack, stroke) | Ages 20+ | Fulgoni et al. 2019, Table 1 |
| Thyroid condition | `MCQ160M = 1` | Ages 20+ | Fulgoni et al. 2019, Table 1 |
| Liver condition | `MCQ160L = 1` | Ages 20+ (`MCQ500`, used only from 2017, is not used) | Fulgoni et al. 2019, Table 1 |
| Cancer | `MCQ220 = 1` | Ages 20+ | Fulgoni et al. 2019, Table 1 |

**Dropped criteria:**

- **hs-CRP > 10 mg/L.** CRP was not measured in 2011-2014, so this criterion cannot be applied consistently across the pooled cycles. We dropped it rather than apply it to some cycles only.
- **Fasting glucose ≥ 126 mg/dL.** It exists only for the fasting subsample. Diabetes is defined by self-report and HbA1c instead.

**Obesity (BMI ≥ 30 for 18+, CDC BMI category `BMDBMIC = 4` for ages 2-17) is a sensitivity flag only.** The reference-interval papers we read do not use BMI as an exclusion: Fulgoni et al. 2019, Table 1, and Hollowell et al. 2002 (cited below). So it is not part of the default subset. The population `reference_bmi_lt30`, in `percentiles_sex_age.csv` only and not in the web JSON, applies the default criteria plus the obesity flag, for sensitivity analysis.

**Sources for the criteria:**

- **Cotinine cut-off.** Benowitz NL, Bernert JT, Caraballo RS, Holiday DB, Wang J. Optimal serum cotinine levels for distinguishing cigarette smokers and nonsmokers within different racial/ethnic groups in the United States between 1999 and 2004. Am J Epidemiol 2009;169(2):236-248. doi:10.1093/aje/kwn301. We verified the DOI on Crossref and read the abstract through Europe PMC (PMID 19019851). The abstract reports optimal cutpoints of 3.08 ng/mL for adults and 2.99 ng/mL for adolescents, and proposes "a new overall cutpoint of 3 ng/mL". We use > 3 ng/mL for ages 12+, the ages the cutpoints were derived for.
- **Kidney thresholds.** Levin A, Ahmed SB, Carrero JJ, et al. Executive summary of the KDIGO 2024 Clinical Practice Guideline for the Evaluation and Management of Chronic Kidney Disease. Kidney Int 2024;105:684-701. doi:10.1016/j.kint.2023.10.016. We read the PDF from kdigo.org. It uses eGFR < 60 mL/min/1.73 m² for GFR categories G3a-G5 (Recommendation 3.15.1.1) and ACR ≥ 30 mg/g as the threshold for albuminuria to confirm (practice points in section 1.3, and Practice Point 1.4.3). The full guideline is: KDIGO 2024 Clinical Practice Guideline for the Evaluation and Management of Chronic Kidney Disease. Kidney Int 2024;105(4S):S117-S314. doi:10.1016/j.kint.2023.10.018. We verified it on Crossref but did not read it in full.

  KDIGO defines CKD by markers persisting for more than 3 months. NHANES measures once, so our flag marks a single elevated measurement, not CKD.
- **eGFR equation.** eGFR uses the race-free CKD-EPI 2021 creatinine equation: Inker LA, Eneanya ND, Coresh J, et al. New creatinine- and cystatin C-based equations to estimate GFR without race. N Engl J Med 2021;385(19):1737-1749. doi:10.1056/NEJMoa2102953, verified on Crossref. We took the coefficients from the National Kidney Foundation page https://www.kidney.org/ckd-epi-creatinine-equation-2021-0.
- **Creatinine bridging.** Serum creatinine for 2011-2016 is first put on the 2017+ Roche Cobas scale with the CDC forward equation from the BIOPRO_J bridging table (0.9515 × Scr + 0.06608). For 2021-2023, CDC's BIOPRO_L bridging table says "Adjustment Not Recommended" for creatinine, so those values are used as released.
- **Fulgoni et al. 2019.** Fulgoni VL, Agarwal S, Kellogg MD, Lieberman HR. Establishing pediatric and adult RBC reference intervals with NHANES data using piecewise regression. Am J Clin Pathol 2019;151(2):128-142. doi:10.1093/ajcp/aqy116.

  Its Table 1 lists:

  - pregnancy and lactation;
  - anemia;
  - coronary heart failure, coronary heart disease, angina and acute myocardial infarction;
  - cancer, gout, emphysema, liver condition and celiac disease;
  - stroke and thyroid condition;
  - hepatitis C and HIV;
  - daily aspirin use;
  - hematologic prescription drugs;
  - donated blood.

  We mirror the items asked the same way in every cycle. We do not mirror:

  - lactation, gout, emphysema/COPD and celiac disease, because the items differ across cycles;
  - hepatitis C, HIV, aspirin, hematologic drugs and blood donation, because they are not in the files we downloaded.
- **Hollowell et al. 2002.** Hollowell JG, Staehling NW, Flanders WD, et al. Serum TSH, T4, and thyroid antibodies in the United States population (1988 to 1994): NHANES III. J Clin Endocrinol Metab 2002;87(2):489-499. doi:10.1210/jcem.87.2.8182. We verified the title and DOI on Crossref but did not read the paper. We cite it as related NHANES reference-population work, not for any criterion.

**Examined participants flagged by the default criteria** (any criterion): 3,587 in 2011-2012 (38%), 3,889 in 2013-2014 (40%), 3,704 in 2015-2016 (39%), 6,079 in 2017-March 2020 (43%), 4,148 in 2021-2023 (47%) and 3,710 in 2017-2018 (43%). Adding the obesity flag raises these to 4,799, 5,215, 5,054, 8,196, 5,446 and 4,487. Counts by criterion are in `data/reference_exclusions.csv`.

Some analytes are themselves used to define the exclusions: HbA1c, serum cotinine, urine albumin/creatinine ratio and serum creatinine. For these, the reference distribution is cut off by design. Their `reference` tails must not be read as natural limits.

The reference subset describes US residents without the conditions above, as recorded in NHANES.

## 8. Validation (`tests/validation.py`, `tests/test_validation.py`)

We fetched three CDC/NCHS publications (`tests/sources.csv`) and copied 62 published numbers into `tests/published_values.csv`, with table and row locations. We then recomputed each one from the raw files with `scripts/survey.py`, the estimator code that builds the release. The pass rules were fixed in `tolerance()` before comparing.

1. **CDC Early Release Blood Lead Table, NHANES 2011-2018** (National Exposure Report). This gives 20 weighted percentiles and 3 geometric means, with 95% CIs, for 2011-12, 2013-14 and 2015-16. Weights are MEC for 2011-12 and the blood-metal subsample weight for 2013-16.
   - All 23 sample sizes match exactly.
   - All 20 percentile point estimates match exactly.
   - All 3 geometric means match to the published 3 digits.
   - 19 of 40 percentile confidence limits match exactly. All 40 are within 5.9% of the published limit.
2. **NCHS Data Brief 515** (Carroll et al. 2024, doi:10.15620/cdc/165796). This gives the crude prevalence of total cholesterol ≥ 240 mg/dL and HDL-C < 40 mg/dL, 2021-2023, adults 20+, by sex and age, with phlebotomy weights. That is 24 estimates.
   - All 24 sample sizes match exactly.
   - The largest difference is 0.046 percentage points.
   - The largest difference in standard error is 0.047 percentage points.
3. **NCHS Data Brief 516** (Gwira et al. 2024, doi:10.15620/cdc/165794). This gives the crude prevalence of total, diagnosed and undiagnosed diabetes, 2021-2023, adults 20+, fasting weights, pregnant women excluded. That is 15 estimates.
   - All 15 sample sizes match exactly.
   - The largest difference is 0.044 percentage points.

Result: 62 of 62 targets pass. These checks test the estimator, the weights, the design variables and the sample selection on single cycles.

4. **CDC National Report on Biochemical Indicators of Diet and Nutrition (2026 release), "Reference intervals" tables.** These test pooled estimates on matching cycles. The tables pool 2011-2014 or 2015-2018. We fetched them on 2026-11-25 from `https://www.cdc.gov/nutritionreport/report/biomarkers//2011-2018/<VAR>/Table-B-<VAR>-<YEARS>.html`, reached from https://www.cdc.gov/nutritionreport/data-tables.html. The values are copied in `tests/nutrition_report_tables.txt`. Tables used:
   - serum 25-hydroxyvitamin D (`LBXVIDMS`), 2011-2014 and 2015-2018;
   - serum ferritin, females 12-49 (`LBXFER`), 2015-2018;
   - RBC folate (`LBDRFOSI`), 2011-2014;
   - serum total folate (`LBDFOTSI`), 2011-2014;
   - serum vitamin B12, ages 20+ (`LBDB12SI`), 2011-2014.

   For every published row (total, age groups, male, female) we recomputed the arithmetic or geometric mean and the 2.5th, 5th, 50th, 95th and 97.5th percentiles. We used the release's weight choice, its domain-wise pooled weights (`scripts/pooling.py`) and its estimators, on the same cycles. That is 258 numbers in all.

   The pass rule, fixed before comparing, is: identical sample size, and our estimate inside CDC's published 95% CI. Results:
   - **258 of 258 pass.**
   - All 258 sample sizes match exactly.
   - 135 of 258 estimates equal the published value at its printed precision.
   - The median absolute relative difference is 0.13%. The largest is 2.5%, a ferritin 2.5th percentile.

   | Table | Numbers | Pass | n match | Largest relative difference |
   |---|---|---|---|---|
   | 25(OH)D 2011-2014 | 54 | 54 | 54 | 1.0% |
   | 25(OH)D 2015-2018 | 54 | 54 | 54 | 1.8% |
   | Ferritin, females 12-49, 2015-2018 | 24 | 24 | 24 | 2.5% |
   | RBC folate 2011-2014 | 54 | 54 | 54 | 1.3% |
   | Serum folate 2011-2014 | 36 | 36 | 36 | 0.9% |
   | Vitamin B12 2011-2014 | 36 | 36 | 36 | 0.6% |

   Results are in `tests/validation_nutrition_report.csv`.

   Our release pools different cycle spans from these tables (for example 2011-2020P for vitamin D). So this checks the pooling machinery on matching spans, not the release's own cells.

## 9. Consistency checks (`tests/test_outputs.py`)

- Percentiles are non-decreasing in every row.
- Every published value meets the sample-size rules in section 6.
- No percentile is published at or below the share of results below LOD.
- Confidence limits bracket their estimates.
- The female and male cells add up to the all-sex cell, for both n and weights.
- Reference n is never larger than everyone n.
- **MEC weights sum to the population.** For G, H, I, P and L, the sum of MEC weights over examined participants equals CDC's ACS control total for the civilian noninstitutionalized population within 0.5%. The control totals are the "Population Totals" files linked from https://wwwn.cdc.gov/nchs/nhanes/ResponseRates.aspx, saved by hand in `../_raw/nhanes/_poptotals/`.
- Every bridging-equation coefficient in `config/analyte_rules.csv` appears in the CDC document named as its source.
- Each curated analyte's weight per cycle equals `expected_weights` in `config/analyte_rules.csv`.
- **Ferritin pooled median:** the female 18+ median is within 2% of 65.7 ng/mL, and band populations add up to the 18+ population for every analyte.
- **Population totals are not understated:** the pooled blood-lead 1-5 cell's `sum_weights` equals the mean of the 2011-2016 single-cycle totals within 1%.
- Cells with coverage below 95% carry no percentiles and are flagged. The 2021-2023 ferritin "all 18+" cell reports `female 18-49`.
- `pct_below_lod` is empty, with a flag, whenever any cycle lacks a CDC comment variable.
- Every default reference criterion is assessable in every cycle.
- The web JSON follows the schema.
- Every `marker_slug` exists in LifeFrom's release 1 `markers.csv` (skipped when that file is absent).

Tests that need the raw NHANES files (section 8, the MEC weight totals, the bridging-equation check and the blood-lead population check) are skipped with a message when `../_raw/nhanes` is absent. All other checks read only the shipped files.

## 10. Rebuild

```
uv sync
uv run python -I scripts/download.py     # about 300 MB into ../_raw/nhanes
uv run python -I scripts/build.py        # about 30 s on 8 cores
uv run python -I -m pytest -q tests      # validation and checks
```

Raw files go to `../_raw/`, a folder next to the repository folder. Set `NHANES_RAW` or `LIFEFROM_RAW` to point at other raw folders. `LIFEFROM_RAW` holds the LifeFrom release 1 `markers.csv` used to fill `marker_slug`. CDC sometimes revises files. Compare `../_raw/nhanes/manifest.csv` with a fresh download to see whether a file changed.
