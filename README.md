# CRC Survival Analysis: SEER

A survival analysis pipeline for colorectal cancer. Walks through the three core methods used in clinical oncology research to analyze time-to-event data.

---

## What This Project Does

It answers questions like:
- "What percentage of stage III CRC patients are still alive at 5 years?"
- "Does Lev+5FU actually improve survival over observation, and by how much?"
- "After controlling for age, sex, and tumor grade, what independently predicts death?"

Those are the same questions that appear in every CRC clinical trial publication. This notebook builds the machinery to answer them.

---

## The Three Methods, Explained Simply

### 1. Kaplan-Meier (KM) Curve
**What it is:** A step-function curve that tracks what fraction of patients are still alive at each point in time.

**The key problem it solves is censoring:**  
Not everyone in a study dies before the study ends. Some patients are still alive at the last follow-up. You can't just throw them out (that biases the estimate toward worse outcomes) and you can't count them as dead (that's wrong). KM handles this correctly: when a patient is censored, they're removed from the denominator going forward without stepping the curve down.

**How to read the output:**  
- Y-axis: probability of still being alive (1.0 = 100%)
- X-axis: time (years here)
- The curve steps down at each death
- The shaded band is the 95% confidence interval, which widens over time as fewer patients remain at risk
- The dashed lines crossing at 50% give you **median survival**, the time at which half the cohort has died

### 2. Log-Rank Test
**What it is:** A statistical test that compares two or more KM curves and asks: "Is the difference between these groups real, or could it be random chance?"

**Output:** A p-value. If p < 0.05, the survival difference between groups is statistically significant.

**Important limitation:** The log-rank test only tells you *whether* groups differ, not *by how much* or *why*. For that, you need Cox regression.

### 3. Cox Proportional Hazards (PH) Regression
**What it is:** A regression model for survival data that estimates each variable's *independent* effect on survival, controlling for everything else simultaneously.

**The key output, the Hazard Ratio (HR):**
- HR = 1.0 → no effect
- HR = 2.0 → twice the risk of death at any given moment ("twice the hazard")
- HR = 0.5 → half the risk → protective

Example from this data: if `nodes` has HR = 1.09, that means each additional positive lymph node is associated with a 9% higher hazard of death, after controlling for age, sex, treatment, and tumor grade.

**The proportional hazards assumption:**  
Cox PH assumes the hazard ratio between groups stays *constant over time* (the survival curves never cross). If two curves cross, this assumption is violated for that variable and the results are unreliable. The notebook checks this with a statistical test (Schoenfeld residuals) after fitting.

**The forest plot:**  
Displays each variable's HR and confidence interval as a horizontal line. The vertical reference line at HR = 1.0 is "no effect." If the confidence interval crosses 1.0, the result is not statistically significant. Variables to the right of 1.0 increase hazard; variables to the left are protective.

---

## The Dataset (Phase 1)

**Source:** NCCTG Colon Cancer Trial (Moertel et al., *New England Journal of Medicine*, 1990)  
**N:** 929 patients with stage III colon cancer  
**Three treatment arms:** Observation (no adjuvant), Levamisole alone, Levamisole + 5-Fluorouracil  
**This trial established Lev+5FU as the standard of care for stage III colon cancer in the 1990s.**

The dataset ships with the `lifelines` Python package and loads from a public URL, so there are no files to download.

**Key columns used:**
| Column | What it means |
|--------|--------------|
| `time` | Days from registration to death or last follow-up |
| `status` / `event` | 1 = died, 0 = censored (still alive) |
| `rx` / `treatment` | Treatment arm: Obs, Lev, or Lev+5FU |
| `age` | Age at diagnosis |
| `sex` | 1 = male, 0 = female |
| `nodes` | Number of positive (cancerous) lymph nodes found at surgery |
| `obstruct` | Bowel obstruction at presentation (1 = yes), worse prognosis |
| `differ` | Tumor differentiation: 1 = well, 2 = moderate, 3 = poor |
| `extent` | How far tumor has grown: 1 = submucosa → 4 = into adjacent structures |
| `etype` | Event type: 1 = recurrence, 2 = death. **We filter to etype == 2** for overall survival |

---

## Notebook Walkthrough

**Cell 1: Imports**  
Loads pandas, matplotlib, seaborn, and the lifelines survival analysis functions.

**Cell 2: Load data**  
Downloads the colon trial dataset from a public URL and prints the raw shape.

**Cell 3: Filter and prep**  
Filters to `etype == 2` (death events only, not recurrence). Converts days → years. Reports cohort size, death rate, and median follow-up.

**Cell 4: Cohort summary**  
Basic descriptive statistics: age, sex, nodes, obstruction rate, grade distribution.

**Cell 5: KM overall survival**  
Single KM curve for the full cohort. Shows median OS with a reference line. Saves as `km_overall.png`.

**Cell 6: KM by treatment arm**  
Three KM curves (Obs / Lev / Lev+5FU) on the same plot. Multivariate log-rank p-value in the corner. Saves as `km_by_treatment.png`.

**Cell 7: Bin lymph nodes into N-stage**  
Creates a new column grouping node counts into N0 (0), N1 (1–3), N2a (4–6), N2b (7+). These mirror the AJCC staging categories used clinically.

**Cell 8: KM by N-stage**  
Four KM curves showing the survival gradient across nodal burden. Saves as `km_by_nodes.png`.

**Cell 9: KM by tumor grade**  
Three KM curves (well / moderate / poorly differentiated). Saves as `km_by_grade.png`.

**Cell 10: Build Cox model matrix**  
Prepares the data for Cox regression: one-hot encodes the treatment column (dropping the reference category "Obs" to avoid multicollinearity), drops rows with missing values, converts any boolean columns to integers.

**Cell 11: Fit Cox model**  
Fits the Cox PH model with a small L2 penalty for numerical stability. Prints the full summary table with HR, 95% CI, and p-value for every predictor.

**Cell 12: Forest plot**  
Visual display of the Cox summary. Each predictor shown as a horizontal line; the dot is the HR estimate, the line is the 95% CI. Saves as `cox_forest.png`.

**Cell 13: PH assumption check**  
Tests whether the proportional hazards assumption holds for each variable. If p < 0.05 for a variable, its effect changes over time and Cox PH is technically not valid for that variable without modification.

**Cell 14: Univariate Cox screen**  
Runs Cox regression separately for each predictor. Standard practice before multivariable modeling, since it shows which variables individually predict survival. Results sorted by p-value.

**Cell 15: Cox-predicted survival curves**  
Uses the fitted Cox model to predict survival for two hypothetical patient profiles: a 65-year-old male with N1 disease on observation vs. the same patient on Lev+5FU. This is the clinical utility of Cox: you can predict for specific patients, not just compare group averages.

**Cell 16: Table 1**  
Standard cohort characteristics table by treatment arm. The format used in every clinical paper: N, death rate, age (median, IQR), sex, nodes, obstruction rate, poor differentiation rate.

---

## Notebook 2: Tree-Based Prediction of 5-Year Mortality (`crc_tree_models.ipynb`)

Reframes the same cohort as a prognostic classification problem: which stage III patients die within 5 years of surgery? Decision tree, random forest, gradient boosting (scikit-learn), and XGBoost are compared against a Cox proportional hazards reference model on a held-out test set.

**Design**
- Outcome: 5-year landmark label. Patients censored alive before 5 years are excluded from the classifiers (14 of 929); the Cox reference model uses full time-to-event data.
- Features: age, sex, positive nodes, differentiation, local extent, obstruction, perforation, adherence, surgery-to-registration interval, treatment arm, a nodes-missingness indicator, and two derived features (AJCC N2 indicator, any-local-complication composite) that are evaluated rather than assumed useful.
- Protocol: stratified 75/25 split; all selection by 5-fold cross-validation or a validation slice inside the training set; test-set AUC with 1000-resample bootstrap intervals, log loss, Brier score, and calibration.

**Results**
- The root split is nodal burden at more than 4 positive nodes, matching the trial's own `node4` covariate and the AJCC N1/N2 boundary within one node.
- Grid search selects depth-1 boosting: an additive model with unconstrained per-feature response and no interactions.
- Random forest, gradient boosting, and XGBoost reach test AUC 0.68 to 0.69 versus 0.65 for Cox, with better log loss and Brier score. Bootstrap intervals overlap, so the advantage is consistent but not statistically established on 229 test patients.
- Partial dependence shows a plateau in predicted mortality from 5 to 10 positive nodes that the proportional hazards form cannot represent.
- Permutation importance on the test set rejects the local-complication composite despite non-trivial impurity importance.

**Appendices:** a minimal L2 boosting implementation reproducing the library models' train/test dynamics, and a check that tree splits are invariant to monotone transforms of the inputs.

**Figures** are saved to `assets/` with the prefixes `tree_`, `boosting_`, `xgboost_`, `roc_`, `calibration_`, `permutation_`, `partial_dependence_`.

**Environment:** Python 3.14 with scikit-learn 1.8, xgboost 3.4, lifelines 0.30, pandas, matplotlib. Data is cached in `data/colon_ncctg.csv` and downloaded once if missing.

---

## Phase 2: Real SEER Data

This notebook is built to swap in real SEER data with minimal changes.

**Steps:**
1. Apply for SEER access at [seer.cancer.gov/data/access.html](https://seer.cancer.gov/data/access.html). Fill out the DUA form. Approval is typically same-day.
2. Install SEER*Stat (Windows). Export a CRC cohort as CSV.
3. Replace the load cell (Cell 2) with the SEER loader documented at the bottom of the notebook.

**Why SEER over this trial dataset?**
- 100,000+ patients vs 929
- Includes MSI status (from 2018 onward), critical for immunotherapy era analysis
- Right-sided vs left-sided colon, colon vs rectum, different tumor biology
- Racial/ethnic disparities with enough N to actually analyze
- Community practice patients, not just academic trial centers

---

## Connection to the Colo Project

The Colo project runs two AI agents that debate CRC literature and produce evidence-graded hypotheses. When they discuss survival outcomes, they're referencing published trial data.

This SEER analysis is the real-world empirical counterpart:
- Trial data = efficacy under controlled conditions
- SEER data = effectiveness in the actual population
- The gap between them is clinically meaningful

The planned Epidemiology Agent C in Colo would run survival analyses like this dynamically and inject population-level survival statistics into the agent debate, grounding theoretical arguments in real numbers.
