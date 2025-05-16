#!/usr/bin/env python3

import pandas as pd
import anndata
import numpy as np
import argparse
import os
from scipy.sparse import csr_matrix

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("output", help="Path to output h5ad file")
    ap.add_argument("--input-dir", required=True, help="Directory containing donor count and metadata files")
    args = vars(ap.parse_args())

    input_dir = args["input_dir"]
    all_adatas = []

    # Find all donor files (assuming _counts and _cells naming pattern)
    donors = set()
    for fname in os.listdir(input_dir):
        if "-counts" in fname:
            donor = fname.replace("_expected-counts.csv.gz", "").replace("-counts.parquet", "")
            donors.add(donor)
    print(donors)

    for donor in sorted(donors):
        print(f"Processing donor: {donor}")
        count_path = os.path.join(input_dir, f"{donor}_expected-counts.csv.gz")
        meta_path = os.path.join(input_dir, f"{donor}_cell-metadata.csv.gz")

        # Handle parquet optionally
        if not os.path.exists(count_path):
            count_path = os.path.join(input_dir, f"{donor}_counts.parquet")
        if not os.path.exists(meta_path):
            meta_path = os.path.join(input_dir, f"{donor}_cells.parquet")

        count_ext = os.path.splitext(count_path)[1]
        meta_ext = os.path.splitext(meta_path)[1]

        count_matrix = pd.read_parquet(count_path) if count_ext == ".parquet" else pd.read_csv(count_path)
        cell_metadata = pd.read_parquet(meta_path) if meta_ext == ".parquet" else pd.read_csv(meta_path)

        obsm = {
            "spatial": np.stack(
                [cell_metadata.centroid_x.to_numpy(), cell_metadata.centroid_y.to_numpy()],
                axis=1
            )
        }
        var = pd.DataFrame(index=count_matrix.columns)
        obs = cell_metadata.set_index(pd.Index(cell_metadata.cell, str, True, "cell_d"))
        X = csr_matrix(count_matrix.to_numpy())

        adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm)
        adata.obs["donor"] = donor
        per_donor_output = os.path.join(input_dir, f"{donor}.h5ad")
        adata.write(per_donor_output)
        print(f"Saved: {per_donor_output}")
        # add to all as a list
        all_adatas.append(adata)
    # Concatenate all donors
    print("Merging all donors into one AnnData...")
    combined = anndata.concat(all_adatas, join="outer", label="donor", keys=None)

    print(f"Writing combined data to {args['output']}")
    combined.write_h5ad(args["output"])
    print("Done.")

