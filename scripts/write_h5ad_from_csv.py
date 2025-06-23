
import pandas as pd
import anndata as ad
import numpy as np


# Read the compressed CSV files
counts = pd.read_csv("/hpc/pmc_holstege/rstudio/liting/spatial/20250502_cite_healthy_BM_ref_Annotated_BM_map_counts_matrix.csv.gz", index_col=0)
normalized = pd.read_csv("/hpc/pmc_holstege/rstudio/liting/spatial/20250502_cite_healthy_BM_ref_Annotated_BM_map_normalized_matrix.csv.gz", index_col=0)
metadata = pd.read_csv("/hpc/pmc_holstege/rstudio/liting/spatial/20250502_cite_healthy_BM_ref_Annotated_BM_map_metadata.csv.gz", index_col=0)

print("all files read-in.")

# If you want to extract gene names (columns) from counts:
gene_names = list(counts.columns)


# Ensure metadata and counts match
print("Counts shape", counts.shape)
print("Gene names", len(gene_names))
print("metadata shape", metadata.shape[0])
#assert counts.shape[0] == metadata.shape[0], "Mismatch: metadata and expression matrix cell count"
#assert counts.shape[1] == len(gene_names), "Mismatch: gene names and expression matrix gene count"

# Add raw counts as a separate layer
adata.layers["counts"] = counts.values
# Normalized: (newly added)
# adata.layers["normalized"] = normalized.values
# Save AnnData
adata.write("/hpc/pmc_holstege/rstudio/liting/spatial/20250502_cite_healthy_BM_ref_Annotated_BM_map.h5ad")


