"""
covid19_aggregate.py
--------------------
Aggregates the RIVM individual-case COVID-19 file
(COVID-19_casus_landelijk_tm_03102021.csv) to a clean daily time series
of reported cases in the Netherlands.

Input  : COVID-19_casus_landelijk_tm_03102021.csv  (semicolon-delimited)
Output : covid19_daily_cases_NL.csv                (comma-delimited, two columns)

Columns in output:
  date   – calendar date, ISO 8601 format (YYYY-MM-DD)
  cases  – number of lab-confirmed COVID-19 cases assigned to that date

Usage:
  python covid19_aggregate.py

Requirements: Python >= 3.8, pandas >= 1.3
"""

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
INPUT_PATH  = "COVID-19_casus_landelijk_tm_03102021.csv"
OUTPUT_PATH = "covid19_daily_cases_NL.csv"

# ── Step 1: load only the columns needed ─────────────────────────────────
# The file is semicolon-delimited (Dutch convention).
# We load only Date_statistics; all other columns are discarded.
df = pd.read_csv(
    INPUT_PATH,
    sep=";",
    usecols=["Date_statistics"],
    dtype=str,
)

# ── Step 2: parse dates ───────────────────────────────────────────────────
# Date_statistics is stored as YYYY-MM-DD strings; parse to datetime.
# errors="coerce" converts any unparseable value to NaT instead of raising.
df["Date_statistics"] = pd.to_datetime(df["Date_statistics"], errors="coerce")

# Safety check: report any rows with unparseable dates.
n_invalid = df["Date_statistics"].isna().sum()
if n_invalid > 0:
    print(f"Warning: {n_invalid} rows have unparseable Date_statistics values "
          f"and will be excluded from the aggregation.")
    df = df.dropna(subset=["Date_statistics"])

# ── Step 3: count records per day ─────────────────────────────────────────
# Each row in the source file represents one confirmed case.
# Grouping by date and counting rows therefore gives daily case counts.
daily = (
    df.groupby("Date_statistics")
    .size()
    .reset_index(name="cases")
)

# ── Step 4: fill gaps in the date sequence with zero ─────────────────────
# Some dates in the early pandemic (January–February 2020) have no reported
# cases. Reindexing to the complete date range makes the time series
# contiguous and explicit about zero-case days, which is required for
# fitting continuous-time epidemic models.
full_range = pd.date_range(
    start=daily["Date_statistics"].min(),
    end=daily["Date_statistics"].max(),
    freq="D",
)
daily = (
    daily.set_index("Date_statistics")
    .reindex(full_range, fill_value=0)
    .reset_index()
    .rename(columns={"index": "date"})
)

# ── Step 5: format date as ISO 8601 string ────────────────────────────────
daily["date"] = daily["date"].dt.strftime("%Y-%m-%d")

# ── Step 6: write output ──────────────────────────────────────────────────
daily.to_csv(OUTPUT_PATH, index=False)

# ── Summary ───────────────────────────────────────────────────────────────
print(f"Input records  : {len(df):>10,}")
print(f"Output rows    : {len(daily):>10,}  (one per calendar day)")
print(f"Date range     : {daily['date'].iloc[0]}  to  {daily['date'].iloc[-1]}")
print(f"Total cases    : {daily['cases'].sum():>10,}")
print(f"Zero-case days : {(daily['cases'] == 0).sum():>10,}  "
      f"(all in Jan–Feb 2020, pre-outbreak)")
print(f"Output saved to: {OUTPUT_PATH}")
