# This file contain all the variables useful for working with the PDCL dataset

PDCL_DATASET_PARAM <- list(
  ##### Control Parameters
  "control_count_path" = "~/PhD/Project/PDCL/Data/pdcl-dataset/control/controls_pdcl_datasets.csv",
  "control_sample_names" = c(
    "AF22_NES_Astro_Br1_d0_01_S1",
    "AF22_NES_Astro_Br1_d29_37_S46",
    "AF22_NES_Astro_Br2_d0_02_S12",
    "AF22_NES_Astro_Br2_d29_38_S56",
    "AF22_NES_Astro_Br3_d0_03_S22",
    "AF22_NES_Astro_Br3_d29_39_S66",
    "C1_NES_Astro_Br1_d0_04_S32",
    "C1_NES_Astro_Br2_d0_05_S42",
    "C1_NES_Astro_Br3_d0_06_S52",
    "C9_NES_Astro_Br1_d0_07_S62",
    "C9_NES_Astro_Br2_d0_08_S72",
    "C9_NES_Astro_Br3_d0_09_S2"
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
  "penda_result" = "penda_results_pdcl.rds"
)