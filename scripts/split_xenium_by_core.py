#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path

def load_bounds(bounds_file):
    return pd.read_csv(bounds_file)

def split_and_save(df, bounds_df, x_col, y_col, prefix, outdir):
    for _, row in bounds_df.iterrows():
        sample = row["Selection"]
        x_min, x_max = row["x_min"], row["x_max"]
        y_min, y_max = row["y_min"], row["y_max"]

        df_sub = df[
            (df[x_col] >= x_min) & (df[x_col] <= x_max) &
            (df[y_col] >= y_min) & (df[y_col] <= y_max)
        ]

        out_path = outdir / f"{prefix}_{sample}.csv.gz"
        df_sub.to_csv(out_path, index=False, compression="gzip")
        print(f"Saved {out_path} ({len(df_sub)} rows)")

def main():
    parser = argparse.ArgumentParser(
        description="Split Xenium transcripts and cells files using bounding boxes."
    )
    parser.add_argument("--transcripts", required=False, help="Path to transcripts.parquet")
    parser.add_argument("--cells", required=False, help="Path to cells.parquet")
    parser.add_argument("--bounds", required=True, help="CSV with Selection,x_min,x_max,y_min,y_max")
    parser.add_argument("--outdir", default="split_output", help="Output directory (default: split_output)")

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading bounding boxes...")
    bounds_df = load_bounds(args.bounds)

    if args.transcripts:
        print("Loading transcripts...")
        transcripts = pd.read_parquet(args.transcripts)

        print("Splitting transcripts by x_location and y_location...")
        split_and_save(transcripts, bounds_df, "x_location", "y_location", "transcripts", outdir)

    if args.cells:
        print("Loading cells...")
        cells = pd.read_parquet(args.cells)

        print("Splitting cells by x_centroid and y_centroid...")
        split_and_save(cells, bounds_df, "x_centroid", "y_centroid", "cells", outdir)

    print("Done.")

if __name__ == "__main__":
    main()

