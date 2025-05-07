import pandas as pd
import argparse

"""
This script reads a CSV file containing coordinate selections and computes the
bounding box (x_min, x_max, y_min, y_max) for each unique selection.

Expected input CSV format:
Selection,X,Y
MDS-RCC1C1,6530.38,6212.10
MDS-RCC1C1,8479.34,6212.10
... (5 rows per selection)

Each selection must consist of exactly 5 rows of (X, Y) coordinates describing a polygon.

Usage:
    python bounding_box_calculator.py --input input_file.csv --output output_file.csv
"""

def compute_bounding_boxes(input_path: str, output_path: str):
    # Load data
    df = pd.read_csv(input_path)

    # Validate required columns
    if not {'Selection', 'X', 'Y'}.issubset(df.columns):
        raise ValueError("Input CSV must contain 'Selection', 'X', and 'Y' columns")

    # Group and compute bounding boxes
    bbox_data = df.groupby('Selection').agg(
        x_min=('X', 'min'),
        x_max=('X', 'max'),
        y_min=('Y', 'min'),
        y_max=('Y', 'max')
    ).reset_index()

    # Save result
    bbox_data.to_csv(output_path, index=False)
    print(f"Bounding boxes saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Compute bounding boxes from coordinate selections")
    parser.add_argument('--input', required=True, help='Path to input CSV file')
    parser.add_argument('--output', required=True, help='Path to output CSV file')
    args = parser.parse_args()

    compute_bounding_boxes(args.input, args.output)


if __name__ == '__main__':
    main()

