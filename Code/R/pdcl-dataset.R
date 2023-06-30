# This file contain all the variables useful for working with the PDCL dataset

PDCL_DATASET_PARAM <- list(
  ##### Control Parameters
  "control_count_path" = "~/PhD/Project/PDCL/Data/pdcl-dataset/control/astrocyte_GSE109001_filter_non_unique_commercial_cell_included.csv",
  "control_sample_names" = c(
    "AF22_NES_Astro_Br1_d29_37_S46",
    "AF22_NES_Astro_Br2_d29_38_S56",
    "AF22_NES_Astro_Br3_d29_39_S66",
    "CCF.STTG1_p24_Br1_S16",
    "CCF.STTG1_p24_Br2_S17",
    "CCF.STTG1_p24_Br3_S18",
    "CDIAstrocytes_p2_Br1_S19",
    "CDIAstrocytes_p2_Br2_S20",
    "CDIAstrocytes_p2_Br3_S21",
    "phaAstrocyte_p2_Br1_S1",
    "phaAstrocyte_p2_Br2_S2",
    "phaAstrocyte_p2_Br3_S3"
  ),
  ##### PDCL Parameters
  "pdcl_count_path" = "~/PhD/Project/PDCL/Data/pdcl-dataset/PDCL/lignees_count_genes_PDCL.csv",
  "pdcl_sample_names" = c(
    "4339-p21",
    "4371-p37",
    "5706-p14",
    "6190-p43",
    "6240-p12",
    "7015-p17",
    "7060-p18", 
    "7142-p14",
    "N13-1300",
    "N13-1520-p9",
    "N14-0072",
    "N14-0870",
    "N14-1208",
    "N14-1525",
    "N15_0460",
    "N15_0516",
    "N15-0385",
    "N15-0661",
    "N16_0535",
    "N16-0240"
  ),
  ##### Penda Parameters
  "penda_result" = "~/PhD/Project/PDCL/Data/pdcl-dataset/results_penda_astrocytes_2021.rds"
)