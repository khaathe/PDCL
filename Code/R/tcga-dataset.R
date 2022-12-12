# This file contain all the variables useful for working with the TCGA dataset

TCGA_DATASET_PARAM <- list(
  "hbm_count_path" = "../Data/tcga-dataset/control/human-brain_map/rna_counts_human_brain_map.csv",
  "hbm_metadata_path" = "../Data/tcga-dataset/control/human-brain_map/metadata_human_brain_map.csv",
  "tcga_count_path" = "../Data/tcga-dataset/TCGA-GBM-v32.0/count_matrix.txt",
  "tcga_metadata_path" = "../Data/tcga-dataset/TCGA-GBM-v32.0/gdc_sample_sheet.2022-07-18.tsv",
  "sample_to_remove" = c(
    "S010492_L8.LB25",
    "S010498_L9.LB3",
    "S010504_L6.LB11",
    "S010508_L5.LB5",
    "S020656_L7.LB18",
    "S020671_L7.LB16",
    "S020671_L7.LB16b",
    "S020697_L1.LB3",
    "S020722_L4.LB25",
    "477a4ae1-84c0-49f2-b3b3-4015b8c26f18"
    ),
  "penda_result" = "penda_result_tcga_gbm.rds"
) 
