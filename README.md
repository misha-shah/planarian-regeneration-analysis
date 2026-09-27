# PFNA & Planarian Regeneration — Statistical Analysis

Python reproduction of the statistical analysis from our published research:

> Casto, E. G., Shah, M. C., & Lorenz, C. M. (2025). *The effects of
> perfluorononanoic acid on pluripotent stem cell regeneration in planaria.*
> Rock Canyon High School, Department of Science.

The study tested whether PFNA (a "forever chemical" found in firefighting
foam and increasingly in Colorado drinking water) inhibits regeneration in
planaria. Planaria are flatworms whose pluripotent neoblast stem cells allow them regenerative capabilities. Fifteen planaria per group were decapitated and regeneration
length was tracked over the course of 5 days across three conditions: a spring-water
control, a low PFNA concentration (29.5 mg/L), and a high PFNA concentration
(59 mg/L). The original analysis used a one-way ANOVA followed by a Fisher's
Least Significant Difference (LSD) post-hoc test. This script rebuilds that
exact pipeline in Python.

## What it does

1. Loads regeneration measurements for the three groups and prints summary
   statistics (n, mean, standard deviation)
2. Runs a one-way ANOVA to test whether regeneration differs across the
   three groups overall, and prints the full ANOVA summary table (SS, df,
   MS, F, p) in the same format as the paper's Table 2
3. Runs a Fisher's LSD post-hoc test** on every pair of groups
   (control vs. low, control vs. high, low vs. high), reporting the mean
   difference, 95% confidence interval, p-value, and whether it's
   significant. This matches the paper's Table 3.
4. Calculates effect size (Cohen's d) for each pair
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
