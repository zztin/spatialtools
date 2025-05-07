import pandas as pd
import argparse
import sys
from pathlib import Path

def convert_parquet_to_csv_gz(parquet_path, output_path):
    print(f"Reading: {parquet_path}")
    df = pd.read_parquet(parquet_path)
    
    print(f"Writing to: {output_path}")
    df.to_csv(output_path, index=False, compression="gzip")
    print("Done.")

def main():
    parser = argparse.ArgumentParser(
        description="Convert transcripts.parquet to transcripts.csv.gz"
    )
    parser.add_argument(
        "input", help="Path to transcripts.parquet"
    )
    parser.add_argument(
        "-o", "--output", default="transcripts.csv.gz", 
        help="Output file path (default: transcripts.csv.gz)"
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: Input file {input_path} not found.", file=sys.stderr)
        sys.exit(1)

    convert_parquet_to_csv_gz(input_path, output_path)

if __name__ == "__main__":
    main()
