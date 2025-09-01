#!/usr/bin/env python3
"""
rank_normalization.py — Rank-based inverse normal transformation for antibody titres.

What it does
------------
- Loads an Excel file with multiplex serology data.
- Applies rank-based inverse normal transformation (INT) to a specified
  set of antibody columns.
- Writes the augmented dataframe (with new *_ranknorm columns) to TSV.

By default, applies to columns 6–15 (0-based indexing), matching your original script.
You can override this with --cols.

Usage
-----
# Default: normalize columns 6–15 from Excel
python rank_normalization.py \
  --input data/Multiplex_Serology_GoalPost_Data.xlsx \
  --output results/Multiplex_Normalized.tsv

# Explicit antibody columns (by name)
python rank_normalization.py \
  --input data/Multiplex_Serology_GoalPost_Data.xlsx \
  --output results/Multiplex_Normalized.tsv \
  --cols Antibody1 Antibody2 Antibody3
"""

from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import rankdata, norm


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Rank-based inverse normal transform for antibody titres.")
    p.add_argument("--input", required=True, help="Input Excel file (multiplex serology).")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--cols", nargs="*", default=None,
                   help="Antibody columns to normalize. "
                        "If omitted, defaults to columns 6–15 (0-based index).")
    return p.parse_args()


def rank_normalize(series: pd.Series) -> pd.Series:
    """Apply rank-based inverse normal transformation."""
    ranks = rankdata(series, method="average")
    n = len(ranks)
    return pd.Series(norm.ppf((ranks - 0.5) / n), index=series.index)


def main():
    args = parse_args()

    # Load input Excel
    df = pd.read_excel(args.input)

    # Determine antibody columns
    if args.cols:
        antibody_cols = args.cols
    else:
        antibody_cols = df.columns[6:16]  # default: cols 6–15 inclusive

    print(f"Applying rank normalization to: {list(antibody_cols)}")

    # Apply rank normalization
    for col in antibody_cols:
        df[f"{col}_ranknorm"] = rank_normalize(df[col])

    # Save
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"✅ Rank normalization complete. Output saved to: {out}")


if __name__ == "__main__":
    main()
