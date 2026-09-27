"""
analysis.py

Reproduces the statistical analysis from:
    Casto, E. G., Shah, M. C., & Lorenz, C. M.
    "The effects of perfluorononanoic acid on pluripotent stem cell
    regeneration in planaria." Rock Canyon High School, 2025.

Three groups of planaria (control / low-concentration PFNA / high-concentration
PFNA) were compared on regeneration length (mm) after a 5-day head-decapitation
regenerative assay. The original paper used a one-way ANOVA followed by a
Fisher's Least Significant Difference (LSD) post-hoc test. This script
rebuilds that exact pipeline in Python.

Written to be readable for someone new to Python -- each section says WHAT
it does and WHY.

How to run:
    python analysis.py

Expected input:
    data/sample_data.csv  (see README.md for an important note on this data)
"""

import itertools

import pandas as pd               # for loading and organizing the data
import numpy as np                # for the math behind the ANOVA/LSD by hand
from scipy import stats           # for the ANOVA F-test and t-distribution
import matplotlib.pyplot as plt   # for the charts


# ---------------------------------------------------------------------------
# 1. LOAD THE DATA
# ---------------------------------------------------------------------------
df = pd.read_csv("data/sample_data.csv")

GROUP_ORDER = ["control", "low", "high"]
LABELS = {"control": "Control (0 mg/L)", "low": "Low (29.5 mg/L)", "high": "High (59 mg/L)"}

groups = {g: df.loc[df["group"] == g, "regeneration_mm"] for g in GROUP_ORDER}

print("Summary statistics by group:")
for g in GROUP_ORDER:
    v = groups[g]
    print(f"  {LABELS[g]:<18} n={len(v):<3} mean={v.mean():.3f} mm  sd={v.std(ddof=1):.3f}")
print()


# ---------------------------------------------------------------------------
# 2. ONE-WAY ANOVA
# ---------------------------------------------------------------------------
# The question: across all three groups, is at least one group's mean
# different from the others by more than we'd expect from random noise?
# ANOVA doesn't say WHICH group differs -- just whether there's a
# significant difference SOMEWHERE among the three. That's what the
# post-hoc LSD test (step 3) is for.
f_stat, p_value = stats.f_oneway(*(groups[g] for g in GROUP_ORDER))

# Rebuild the classic ANOVA summary table (sum of squares, df, mean square)
# the same way the original paper reports it, so the numbers are checkable
# line by line against Table 2 in the paper.
all_vals = df["regeneration_mm"]
grand_mean = all_vals.mean()

ss_between = sum(len(groups[g]) * (groups[g].mean() - grand_mean) ** 2 for g in GROUP_ORDER)
ss_within = sum((len(groups[g]) - 1) * groups[g].var(ddof=1) for g in GROUP_ORDER)
ss_total = ss_between + ss_within

df_between = len(GROUP_ORDER) - 1
df_within = len(all_vals) - len(GROUP_ORDER)

ms_between = ss_between / df_between
ms_within = ss_within / df_within

print("One-way ANOVA summary table:")
print(f"  {'Source':<10}{'df':<6}{'SS':<10}{'MS':<10}{'F':<10}{'p'}")
print(f"  {'Group':<10}{df_between:<6}{ss_between:<10.3f}{ms_between:<10.3f}{f_stat:<10.3f}{p_value:.4g}")
print(f"  {'Error':<10}{df_within:<6}{ss_within:<10.3f}{ms_within:<10.3f}")
print(f"  {'Total':<10}{df_between+df_within:<6}{ss_total:<10.3f}")
print()

if p_value < 0.05:
    print(f"  -> At least one group differs significantly (p = {p_value:.4g})")
else:
    print(f"  -> No significant difference detected across groups (p = {p_value:.4g})")
print()


# ---------------------------------------------------------------------------
# 3. POST-HOC: FISHER'S LEAST SIGNIFICANT DIFFERENCE (LSD) TEST
# ---------------------------------------------------------------------------
# ANOVA told us SOMETHING differs. Fisher's LSD checks EVERY pair of groups
# (control vs. low, control vs. high, low vs. high) individually, using the
# pooled error variance (MS_within) from the ANOVA above -- this is more
# powerful than running three separate t-tests because it borrows strength
# from all three groups' data to estimate the noise level.
t_crit = stats.t.ppf(0.975, df_within)  # two-tailed, 95% confidence

print("Fisher's LSD pairwise comparisons (95% confidence):")
print(f"  {'Comparison':<20}{'Diff (mm)':<12}{'95% CI':<22}{'p-value':<10}{'Significant?'}")

for g1, g2 in itertools.combinations(GROUP_ORDER, 2):
    n1, n2 = len(groups[g1]), len(groups[g2])
    mean1, mean2 = groups[g1].mean(), groups[g2].mean()
    diff = mean1 - mean2

    se = np.sqrt(ms_within * (1 / n1 + 1 / n2))
    margin = t_crit * se
    ci_low, ci_high = diff - margin, diff + margin

    t_stat = diff / se
    p_pair = 2 * stats.t.sf(abs(t_stat), df_within)

    # "Significant" here means the 95% CI doesn't cross zero
    sig = "Yes" if (ci_low > 0 or ci_high < 0) else "No"

    label = f"{g1} vs {g2}"
    ci_str = f"({ci_low:.3f}, {ci_high:.3f})"
    print(f"  {label:<20}{diff:<12.3f}{ci_str:<22}{p_pair:<10.4f}{sig}")
print()


# ---------------------------------------------------------------------------
# 4. EFFECT SIZE (Cohen's d) FOR EACH PAIR
# ---------------------------------------------------------------------------
# Significance alone doesn't say how BIG the effect is. Cohen's d expresses
# the gap between two group means in standardized units.
# Rough guide: 0.2 = small, 0.5 = medium, 0.8+ = large.
print("Effect sizes (Cohen's d):")
for g1, g2 in itertools.combinations(GROUP_ORDER, 2):
    n1, n2 = len(groups[g1]), len(groups[g2])
    pooled_sd = np.sqrt(
        ((n1 - 1) * groups[g1].var(ddof=1) + (n2 - 1) * groups[g2].var(ddof=1)) / (n1 + n2 - 2)
    )
    d = (groups[g1].mean() - groups[g2].mean()) / pooled_sd
    print(f"  {g1} vs {g2}: d = {d:.2f}")
print()


# ---------------------------------------------------------------------------
# 5. PERCENT DIFFERENCE VS. CONTROL
# ---------------------------------------------------------------------------
# The paper reports how much greater the control's regeneration was than
# each PFNA group, in percentage terms. Reproduced here from this dataset.
control_mean = groups["control"].mean()
for g in ["low", "high"]:
    pct = (control_mean - groups[g].mean()) / groups[g].mean() * 100
    print(f"Control mean is {pct:.1f}% greater than the {g}-concentration group")
print()


# ---------------------------------------------------------------------------
# 6. CHARTS
# ---------------------------------------------------------------------------
# (a) Boxplot -- shows spread/median per group, like Figure 14 in the paper.
fig1, ax1 = plt.subplots(figsize=(6, 5))
ax1.boxplot(
    [groups[g] for g in GROUP_ORDER],
    tick_labels=[LABELS[g] for g in GROUP_ORDER],
    patch_artist=True,
)
ax1.axhline(0, color="gray", linewidth=0.8, linestyle="--")
ax1.set_ylabel("Regeneration, Day 5 (mm)")
ax1.set_title("Distribution of Planarian Regeneration by Group")
fig1.tight_layout()
fig1.savefig("output/regeneration_boxplot.png", dpi=150)

# (b) Bar chart with error bars (SD) -- like Figure 15 in the paper.
fig2, ax2 = plt.subplots(figsize=(6, 5))
means = [groups[g].mean() for g in GROUP_ORDER]
sds = [groups[g].std(ddof=1) for g in GROUP_ORDER]
ax2.bar([LABELS[g] for g in GROUP_ORDER], means, yerr=sds, capsize=6)
ax2.axhline(0, color="gray", linewidth=0.8)
ax2.set_ylabel("Mean regeneration, Day 5 (mm)")
ax2.set_title("Mean Regeneration by Group (error bars = SD)")
fig2.tight_layout()
fig2.savefig("output/regeneration_bar_chart.png", dpi=150)

print("Charts saved to output/regeneration_boxplot.png and output/regeneration_bar_chart.png")
