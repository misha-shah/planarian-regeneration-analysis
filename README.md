# PFNA & Planarian Regeneration — Statistical Analysis

Python reproduction of the statistical analysis from our published research:

> Casto, E. G., Shah, M. C., & Lorenz, C. M. (2025). *The effects of
> perfluorononanoic acid on pluripotent stem cell regeneration in planaria.*
> Rock Canyon High School, Department of Science.

The study tested whether PFNA (a "forever chemical" found in firefighting
foam and increasingly in Colorado drinking water) inhibits regeneration in
planaria — flatworms whose pluripotent neoblast stem cells let them regrow
a severed head. Fifteen planaria per group were decapitated and regeneration
length was tracked over 5 days across three conditions: a spring-water
control, a low PFNA concentration (29.5 mg/L), and a high PFNA concentration
(59 mg/L). The original analysis used a one-way ANOVA followed by a Fisher's
Least Significant Difference (LSD) post-hoc test — this script rebuilds that
exact pipeline in Python.

## What it does

1. Loads regeneration measurements for the three groups and prints summary
   statistics (n, mean, standard deviation)
2. Runs a **one-way ANOVA** to test whether regeneration differs across the
   three groups overall, and prints the full ANOVA summary table (SS, df,
   MS, F, p) in the same format as the paper's Table 2
3. Runs a **Fisher's LSD post-hoc test** on every pair of groups
   (control vs. low, control vs. high, low vs. high), reporting the mean
   difference, 95% confidence interval, p-value, and whether it's
   significant — matching the paper's Table 3
4. Calculates **effect size (Cohen's d)** for each pair
5. Reproduces the paper's "control was X% greater than group Y" figures
6. Generates two charts: a boxplot (matching Figure 14) and a bar chart with
   SD error bars (matching Figure 15)

## Project structure

```
planarian-regeneration-analysis/
├── analysis.py           # main script — run this
├── data/
│   └── sample_data.csv   # reconstructed data (see note below — important)
├── output/
│   ├── regeneration_boxplot.png
│   └── regeneration_bar_chart.png
├── requirements.txt
└── README.md
```

## Setup

```bash
git clone <your-repo-url>
cd planarian-regeneration-analysis
pip install -r requirements.txt
```

## Run it

```bash
python analysis.py
```

## ⚠️ About the data — read this before showing it to anyone

**`data/sample_data.csv` is reconstructed, not the raw underlying
measurements.** The paper reports only summary statistics per group (n,
mean, SD, quartiles — its Table 4), not the individual per-planarian
readings. This script generates synthetic per-planarian values that are
mathematically forced to match the paper's reported n, mean, and standard
deviation *exactly* for each group. Rerunning the sanity check below shows
those reconstructed numbers reproduce the paper's own ANOVA table (F ≈
8.69 vs. the paper's 8.673; identical 95% confidence intervals in the LSD
table) — so the *statistical conclusions* are a faithful reproduction, but
the individual data points themselves are simulated, not the actual
per-trial measurements, and won't match the exact shape of the true
distribution (e.g., the real data may have been more skewed, since some
planaria disintegrated during the high/low trials).

Be upfront about this if you show it to anyone: it's "I rebuilt my
published analysis pipeline and validated it against my own reported
results," not "here is our raw lab data."

If you still have the real per-planarian measurements from the lab
notebook, replace `data/sample_data.csv` with them directly — the two
required columns are:

| column              | meaning                                   |
|---------------------|--------------------------------------------|
| `group`             | `control`, `low`, or `high`                 |
| `regeneration_mm`   | change in length from day 0 to day 5 (mm)  |

## Why ANOVA + Fisher's LSD (not a simple t-test)

A t-test compares exactly two groups. With three groups (control, low,
high), running three separate t-tests inflates the false-positive rate.
ANOVA tests all three at once for *any* difference, and Fisher's LSD then
locates *which* pairs differ — using the pooled error variance from the
ANOVA (`MS_within`) rather than re-estimating noise separately for each
pair, which is more powerful when sample sizes are small, as they are here.

## Possible extensions

- Add a Tukey HSD comparison alongside Fisher's LSD (Tukey is more
  conservative and corrects for multiple comparisons more strictly)
- Model regeneration as a function of day (1–5) instead of a single Day-5
  endpoint, using a repeated-measures or mixed-effects approach
- Add a script that scrapes/loads EPA PFAS drinking-water monitoring data
  and compares it against the concentrations tested here
