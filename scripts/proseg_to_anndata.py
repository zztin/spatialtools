#!/usr/bin/env python3

## Code has not been tested yet. Drafted together with chatGPT.
import os
import glob
import pandas as pd
import numpy as np
from anndata import AnnData
import argparse

def load_and_aggregate_transcripts(transcript_dir):
    print("Loading transcript metadata...")
    files = glob.glob(os.path.join(transcript_dir, "*_transcript-metadata.csv.gz"))
    all_dfs = []
    
    for f in files:
        sample_name = os.path.basename(f).replace("_transcript-metadata.csv.gz", "")
        df = pd.read_csv(f)
        df["sample_name"] = sample_name
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, ignore_index=True)

    # Keep only assigned transcripts (ProSeg uses 4294967295 for unassigned)
    df_all = df_all[df_all["assignment"] != 4294967295]

    # Group by sample + cell + gene to get counts
    counts = (
        df_all.groupby(["sample_name", "assignment", "gene"])
        .size()
        .reset_index(name="count")
    )

    # Pivot to wide format
    expr_matrix = counts.pivot_table(
        index=["sample_name", "assignment"],
        columns="gene",
        values="count",
        fill_value=0
    ).reset_index()

    return expr_matrix

def load_cell_metadata(cell_dir):
    print("Loading cell metadata...")
    files = glob.glob(os.path.join(cell_dir, "*_cell-metadata.csv.gz"))
    all_dfs = []

    for f in files:
        sample_name = os.path.basename(f).replace("_cell-metadata.csv.gz", "")
        df = pd.read_csv(f)
        df["sample_name"] = sample_name
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, ignore_index=True)

    return df_all

def build_anndata(expr_matrix, cell_metadata):
    print("Merging metadata and expression...")
    # Create unique cell ID index
    expr_matrix["unique_id"] = expr_matrix["sample_name"] + "_" + expr_matrix["assignment"].astype(str)
    cell_metadata["unique_id"] = cell_metadata["sample_name"] + "_" + cell_metadata["original_cell_id"].astype(str)

    # Set index for merging
    expr_matrix = expr_matrix.set_index("unique_id")
    cell_metadata = cell_metadata.set_index("unique_id")

    # Match ordering
    cell_metadata = cell_metadata.loc[expr_matrix.index]

    # Extract gene expression matrix
    gene_cols = [col for col in expr_matrix.columns if col not in ["sample_name", "assignment"]]
    X = expr_matrix[gene_cols]

    # Build AnnData
    adata = AnnData(X=X.values, obs=cell_metadata.copy(), var=pd.DataFrame(index=gene_cols))
    adata.obs_names = expr_matrix.index
    adata.var_names.name = "gene"
    adata.obs_names.name = "cell_id"

    # Add spatial coordinates if available
    if "centroid_x" in adata.obs.columns and "centroid_y" in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[["centroid_x", "centroid_y"]].to_numpy()

    return adata

def main():
    parser = argparse.ArgumentParser(description="Build AnnData from ProSeg output.")
    parser.add_argument("--input-dir", required=True, help="Directory with ProSeg *_transcript-metadata.csv.gz and *_cell-metadata.csv.gz files")
    parser.add_argument("--output", default="adata.h5ad", help="Output .h5ad file (default: adata.h5ad)")
    args = parser.parse_args()

    expr_matrix = load_and_aggregate_transcripts(args.input_dir)
    cell_metadata = load_cell_metadata(args.input_dir)
    adata = build_anndata(expr_matrix, cell_metadata)

    print(f"Saving AnnData to: {args.output}")
    adata.write(args.output)
    print("Done.")

if __name__ == "__main__":
    main()

