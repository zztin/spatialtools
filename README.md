# spatialtools

A set of scripts to aid spatial trnascriptomics (Xenium)  preprocessing.

## Output files description

1. coordinates_all_cores.csv:  Export from Xenium Explorer 3.2.0: Selections -> Retangular selection tools -> Rename Selection -> Download all selections as csv.
2. coordinates_min_max.csv: Converted coordinates from convert_XE_selection_to_min_max_coor.py
3. transcripts.csv.gz: Converted from Parquet to csv.gz with pandas using convert_parquet_csv.py 
4. transcripts_XXXX.csv.gz: Created by split_xenium_by_core.py
5. cells_XXXX.csv.gz: Created by split_xenium_by_core.py


