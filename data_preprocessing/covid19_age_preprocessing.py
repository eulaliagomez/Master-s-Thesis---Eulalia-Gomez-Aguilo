"""
covid19_age_preprocessing.py

Preprocesses the RIVM COVID-19 case line-list into weekly age-stratified
counts for two age groups:
    young: < 60 years  (RIVM groups: 0-9, 10-19, 20-29, 30-39, 40-49, 50-59)
    old:   >= 60 years (RIVM groups: 60-69, 70-79, 80-89, 90+)

Records with Agegroup in {'<50', 'Unknown'} are excluded (<0.02% of records).

Output files (same working directory):
    covid19_weekly_young_NL.csv   — weekly cases for <60 group
    covid19_weekly_old_NL.csv    — weekly cases for >=60 group
    covid19_weekly_allages_NL.csv — total (for cross-check with original)
    covid19_age_summary.txt       — sanity checks and wave-level totals

Usage:
    python covid19_age_preprocessing.py

Requires:
    COVID-19_casus_landelijk_tm_03102021.csv  (RIVM case line-list)
    pandas, numpy

Sources:
    RIVM COVID-19 case line-list:
        https://data.rivm.nl/covid-19/
    Age-group split justification:
        CBS mid-year 2020 population: N_young=12,736,000; N_old=3,359,000
        Age 60 chosen as cut-off to align with Dutch vaccination priority groups
        and with the Prem et al. 2021 contact matrix used in the SEIR model.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_FILE = "COVID-19_casus_landelijk_tm_03102021.csv"
YOUNG_GROUPS = {"0-9", "10-19", "20-29", "30-39", "40-49", "50-59"}
OLD_GROUPS   = {"60-69", "70-79", "80-89", "90+"}
EXCL_GROUPS  = {"<50", "Unknown"}  # ambiguous or missing — excluded

# Wave boundaries (inclusive, week-start dates)
WAVES = {
    "wave1": ("2020-02-24", "2020-07-05"),
    "wave2": ("2020-09-07", "2021-01-17"),
    "wave3": ("2021-02-01", "2021-06-20"),
    "wave4": ("2021-06-21", "2021-10-03"),
}

# Population sizes (CBS mid-year 2020)
N_YOUNG = 12_736_000
N_OLD   =  3_359_000

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading RIVM data...")
df = pd.read_csv(INPUT_FILE, sep=";", low_memory=False,
                 parse_dates=["Date_statistics"])
print(f"  Total records: {len(df):,}")

# ── Assign age group ──────────────────────────────────────────────────────────
df["age_cat"] = df["Agegroup"].map(
    {g: "young" for g in YOUNG_GROUPS} |
    {g: "old"   for g in OLD_GROUPS}   |
    {g: "excl"  for g in EXCL_GROUPS}
)

n_excl = (df["age_cat"] == "excl").sum()
print(f"  Excluded (ambiguous/unknown age): {n_excl:,} ({n_excl/len(df)*100:.3f}%)")
df = df[df["age_cat"] != "excl"].copy()
print(f"  Usable records: {len(df):,}")

# ── Filter: Date_statistics_type (use DOO preferred, as in the original analysis)
# Keep all Date_statistics_type values as the original script did —
# the date field already represents the best-estimate case date.

# ── Weekly aggregation (ISO weeks starting Monday) ────────────────────────────
df["year_week"] = df["Date_statistics"].dt.to_period("W-MON")

weekly_all = (
    df.groupby("year_week")
      .size()
      .reset_index(name="cases_total")
)
weekly_y = (
    df[df["age_cat"] == "young"]
      .groupby("year_week")
      .size()
      .reset_index(name="cases_young")
)
weekly_o = (
    df[df["age_cat"] == "old"]
      .groupby("year_week")
      .size()
      .reset_index(name="cases_old")
)

# Merge all three
weekly = (
    weekly_all
    .merge(weekly_y, on="year_week", how="left")
    .merge(weekly_o, on="year_week", how="left")
    .fillna(0)
    .assign(
        cases_young = lambda x: x["cases_young"].astype(int),
        cases_old   = lambda x: x["cases_old"].astype(int),
        week_start  = lambda x: x["year_week"].apply(
            lambda p: p.start_time.date()
        )
    )
    .sort_values("week_start")
    .reset_index(drop=True)
)

weekly["check"] = weekly["cases_young"] + weekly["cases_old"]
weekly["diff"]  = weekly["cases_total"] - weekly["check"]

print(f"\nWeekly aggregation complete: {len(weekly)} weeks")
print(f"  Max discrepancy (should be 0): {weekly['diff'].abs().max()}")

# ── Extract waves ──────────────────────────────────────────────────────────────
def extract_wave(df_w, start, end):
    mask = (df_w["week_start"] >= pd.Timestamp(start).date()) & \
           (df_w["week_start"] <= pd.Timestamp(end).date())
    out = df_w[mask].copy().reset_index(drop=True)
    out["week"] = range(1, len(out) + 1)
    return out[["week", "week_start", "cases_young", "cases_old", "cases_total"]]

waves_data = {name: extract_wave(weekly, s, e) for name, (s, e) in WAVES.items()}

# ── Save individual wave files ────────────────────────────────────────────────
for wname, wdf in waves_data.items():
    # Young
    out_y = wdf[["week", "week_start", "cases_young"]].rename(
        columns={"cases_young": "reports"})
    out_y.to_csv(f"covid19_{wname}_young_NL.csv", index=False)
    # Old
    out_o = wdf[["week", "week_start", "cases_old"]].rename(
        columns={"cases_old": "reports"})
    out_o.to_csv(f"covid19_{wname}_old_NL.csv", index=False)
    print(f"  Saved {wname}: T={len(wdf)} weeks")

# Full weekly series (for reference)
weekly[["week_start","cases_young","cases_old","cases_total"]].to_csv(
    "covid19_weekly_allages_NL.csv", index=False)
print("\nSaved covid19_weekly_allages_NL.csv")

# ── Summary report ────────────────────────────────────────────────────────────
lines = []
lines.append("=" * 65)
lines.append("COVID-19 Netherlands — Age-stratified preprocessing summary")
lines.append("=" * 65)
lines.append(f"\nInput:  {INPUT_FILE}")
lines.append(f"Records: {len(df):,} usable ({n_excl} excluded)")
lines.append(f"\nPopulation split (CBS 2020 mid-year):")
lines.append(f"  N_young (<60):  {N_YOUNG:>12,}  ({N_YOUNG/(N_YOUNG+N_OLD)*100:.1f}%)")
lines.append(f"  N_old   (>=60): {N_OLD:>12,}  ({N_OLD/(N_YOUNG+N_OLD)*100:.1f}%)")

lines.append("\n" + "-" * 65)
lines.append(f"{'Wave':<8} {'T':>4} {'Young total':>12} {'Old total':>12} "
             f"{'Young peak':>12} {'Old peak':>12}")
lines.append("-" * 65)

for wname, wdf in waves_data.items():
    lines.append(
        f"{wname:<8} {len(wdf):>4} "
        f"{wdf['cases_young'].sum():>12,} "
        f"{wdf['cases_old'].sum():>12,} "
        f"{wdf['cases_young'].max():>12,} "
        f"{wdf['cases_old'].max():>12,}"
    )

lines.append("\n" + "-" * 65)
lines.append("Week-by-week detail per wave:")
for wname, wdf in waves_data.items():
    lines.append(f"\n{wname.upper()}")
    lines.append(f"  {'wk':>3}  {'date':>12}  {'young':>8}  {'old':>8}  {'total':>8}  {'%old':>6}")
    for _, row in wdf.iterrows():
        pct = row["cases_old"] / row["cases_total"] * 100 if row["cases_total"] > 0 else 0
        lines.append(
            f"  {int(row['week']):>3}  {str(row['week_start']):>12}  "
            f"{int(row['cases_young']):>8,}  {int(row['cases_old']):>8,}  "
            f"{int(row['cases_total']):>8,}  {pct:>5.1f}%"
        )

report = "\n".join(lines)
print("\n" + report)
with open("covid19_age_summary.txt", "w") as f:
    f.write(report)
print("\nSaved covid19_age_summary.txt")
