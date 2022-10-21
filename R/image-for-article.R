library(ggplot2)
library(plotly)
library(tidyverse)
library(reticulate)
reticulate::use_python("/usr/bin/python3")
reticulate::py_run_string("import sys")

setwd("~/PhD/PDCL/")
source("R/process-result.R")

all_enrich <- read.csv("all_results.csv", as.is = T)

glioma_pathways <- all_enrich %>% 
  dplyr::filter(grepl("\\bglioma\\b", description, ignore.case = T, perl = T)) %>% 
  dplyr::select(!genes)

data <- glioma_pathways %>% dplyr::filter( !(sample %in% c("tcga_gbm", "gbm")) )

glioma_fdr_boxplot <- ggplot(data, aes(x = method, y = fdr)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0.05, color = "red") + 
  geom_hline(yintercept = 0.1, color = "blue") + 
  labs(
    title = "FDR of Glioma pathways in enrichment", 
    x = "Enrichment Tools", 
    y = "False Discovery Rate (FDR)"
  )

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/glioma-fdr-boxplot.png", plot = glioma_fdr_boxplot, width = 7, height = 7)

global_and_tcga <- all_enrich %>%
  dplyr::filter( sample %in% c("PDCL", "TCGA-GBM") )

common_global_tcga <- global_and_tcga %>%
  dplyr::inner_join(global_and_tcga, by = c("pathway_id", "description", "sample", "db")) %>%
  dplyr::filter(method.x == "G:Profiler", method.y == "GSEA") %>%
  dplyr::rename(fdr_gprofiler = fdr.x, fdr_gsea = fdr.y) %>%
  dplyr::select(pathway_id, description, sample, db, fdr_gprofiler, fdr_gsea)

common_global_tcga <- common_global_tcga %>%
  dplyr::inner_join(common_global_tcga, by = c("pathway_id", "description", "db")) %>%
  dplyr::filter(sample.x == "PDCL", sample.y == "TCGA-GBM") %>%
  dplyr::rename(
    fdr_gprofiler_pdcl = fdr_gprofiler.x, 
    fdr_gsea_pdcl = fdr_gsea.x, 
    fdr_gprofiler_tcga = fdr_gprofiler.y,
    fdr_gsea_tcga = fdr_gsea.y
  ) %>%
  dplyr::select(pathway_id, description, db, fdr_gprofiler_pdcl, fdr_gprofiler_tcga, fdr_gsea_pdcl, fdr_gsea_tcga)

threshold <- 0.1
common_global_tcga <- common_global_tcga %>% 
  dplyr::filter(fdr_gprofiler_pdcl < threshold, fdr_gprofiler_tcga < threshold, fdr_gsea_pdcl < threshold, fdr_gsea_tcga < threshold)

kegg_gs_metadata <- read.csv("Data/gene-set/kegg_geneset_metadata.csv", as.is = T)
reactome_gs_metadata <- read.csv("Data/gene-set/reactome_geneset_metadata.csv", as.is = T)

gs_metadata <- rbind(kegg_gs_metadata, reactome_gs_metadata)
row.names(gs_metadata) <- gs_metadata$pathway_id

common_global_tcga$category <- gs_metadata[common_global_tcga$pathway_id, "category"]
common_global_tcga <- common_global_tcga %>% dplyr::arrange(category, description)

data <- global_and_tcga[global_and_tcga$pathway_id %in% common_global_tcga$pathway_id, ]
data$description <- factor(data$description, levels = common_global_tcga$description, ordered = T)

heatmap_global_tcga <- ggplot(data, aes(x = method, y = description, fill = fdr)) +
  geom_tile() +
  facet_wrap(~sample) +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  labs(
    title = "FDR for Pathway Common between PDCL and TCGA datasets",
    x = "Pathway Enrichment Tool",
    y = "Pathways",
    fill = "FDR"
  )

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/heatmap-fdr-global-tcga.png", plot = heatmap_global_tcga, width = 10, height = 7)

common_pdcl <- all_enrich %>%  
  dplyr::filter(
    !(sample %in% c("TCGA-GBM", "PDCL")),
    fdr < 0.1
  ) %>%
  dplyr::count(pathway_id, description, sample) %>%
  dplyr::filter(n==2)

frequently_dereg <- common_pdcl %>%
  dplyr::count(pathway_id, description) %>%
  dplyr::inner_join(gs_metadata, by = c("pathway_id", "description")) %>%
  dplyr::filter( !grepl(pattern = "disease", x = category, perl = T, ignore.case = T) ) %>%
  dplyr::arrange(desc(n))

freq_dereg_enrich <- all_enrich %>%
  dplyr::filter(pathway_id %in% frequently_dereg$pathway_id[frequently_dereg$n>2], !(sample %in% c("TCGA-GBM", "PDCL")) )

data_commons <- freq_dereg_enrich %>%
  dplyr::filter(fdr<0.1) %>%
  dplyr::count(pathway_id, description, sample) %>%
  dplyr::filter(n>1)

merged_categ <- c(
  "Autophagy" = "Autophagy",
  "Cell Cycle" = "Cell Cycle",
  "Cell-Cell communication" = "Cell-Cell communication",
  "Cellular responses to stimuli" = "Cellular responses to stimuli",
  "Chromatin organization" = "DNA Processess",
  "Circadian Clock" = "Circadian Clock",
  "DNA Repair" = "DNA Processess",
  "DNA Replication" = "DNA Processess",
  "Developmental Biology" = "Developmental Biology",
  "Digestion and absorption" = "Organismal Systems",
  "Disease" = "Disease",
  "Drug ADME" = "Drug",
  "Extracellular matrix organization" = "Extracellular matrix organization",
  "Gene expression (Transcription)" = "DNA Processess",
  "Hemostasis" = "Hemostasis",
  "Immune System" = "Organismal Systems",
  "Metabolism" = "Metabolism",
  "Metabolism of RNA" = "Metabolism of RNA and proteins",
  "Metabolism of proteins" = "Metabolism of RNA and proteins",
  "Muscle contraction" = "Muscle contraction",
  "Neuronal System" = "Organismal Systems",
  "Organelle biogenesis and maintenance" = "",
  "Programmed Cell Death" = "Cellular Processes",
  "Protein localization" = "Protein localization",
  "Reproduction" = "Reproduction",
  "Sensory Perception" = "Organismal Systems",
  "Signal Transduction" = "Signal Transduction and Transport",
  "Transport of small molecules" = "Signal Transduction and Transport",
  "Vesicle-mediated transport" = "Signal Transduction and Transport",
  "Cellular Processes" = "Cellular Processes",
  "Drug Development" = "Drug",
  "Environmental Information Processing" = "Signal Transduction and Transport",
  "Genetic Information Processing" = "DNA Processess",
  "Human Diseases" = "Disease",
  "Metabolism" = "Metabolism",
  "Organismal Systems" = "Organismal Systems",
  "Unknown" = "Unknown"
)

pdcl_names <- c(
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
)

data <- common_pdcl %>%
  dplyr::inner_join(gs_metadata, by = c("pathway_id", "description")) %>%
  dplyr::mutate(category_merged = merged_categ[category], sample = factor(sample, levels = pdcl_names) )
# We remove the only unknown pâthway as it complikcates the viewing and does not
# really affect the result
# its  KEGG entry is path:hsa01240 - Biosynthesis of cofactors
data <- data[data$category_merged != "Unknown",]

barplot_categ_pdcl <- ggplot(data, aes(sample, fill = category_merged)) +
  geom_bar() +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90) ) +
  labs(
    title = "Count of deregulated pathways by samples and category",
    x = "Sample Name",
    y = "Count Pathways",
    fill = "Biological Category"
  ) +
  scale_x_discrete(drop = F)

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/barplot-categ-pdcl.png", plot = barplot_categ_pdcl, width = 12, height = 7)

data <- common_pdcl %>%
  dplyr::inner_join(gs_metadata, by = c("pathway_id", "description")) %>%
  dplyr::mutate(category_merged = merged_categ[category], sample = factor(sample, levels = pdcl_names) ) %>%
  dplyr::count(category_merged)

piechart_categ_pdcl <- plot_ly(data, labels =~category_merged, values =~n, type = 'pie') %>%
  layout(
    title = "Count of the number of pathways per category"
  )

# save_image(piechart_categ_pdcl, "/home/spinicck/PhD/paper/analysis-pdcl/manuscript/img/piechart-categ-pdcl.png", scale = 2.0)

selected_pathways <- c(
  "path:hsa04115",
  "R-HSA-69563",
  "R-HSA-109581",
  "R-HSA-1650814",
  "R-HSA-8948216",
  "R-HSA-1474290",
  "path:hsa04512",
  "R-HSA-1566948",
  "path:hsa04510",
  "path:hsa00190",
  "R-HSA-5389840",
  "R-HSA-611105",
  "R-HSA-69242",
  "R-HSA-69206",
  "R-HSA-69481",
  "R-HSA-73884",
  "R-HSA-191273",
  "path:hsa00100",
  "R-HSA-5358508",
  "path:hsa03430",
  "R-HSA-5693532",
  "path:hsa03030",
  "R-HSA-69002",
  "R-HSA-69190"
)

data <- all_enrich %>%
  dplyr::filter(pathway_id %in% selected_pathways, !(sample %in% c("TCGA-GBM", "PDCL"))) %>%
  dplyr::mutate(category = gs_metadata[pathway_id, "category"]) %>%
  dplyr::arrange(category, description)

data$description <- sub(" - Homo sapiens \\(human\\)", "", data$description, ignore.case = T, perl = T)
data$description <- str_c( data$description, data$pathway_id, sep = " - ")

data_commons <- common_pdcl[common_pdcl$pathway_id %in% selected_pathways, ]
data_commons$description <- sub(" - Homo sapiens \\(human\\)", "", data_commons$description, ignore.case = T, perl = T)
data_commons$description <- str_c( data_commons$description, data_commons$pathway_id, sep = " - ")
data_commons <- data_commons %>%
  dplyr::mutate(category = gs_metadata[pathway_id, "category"]) %>%
  dplyr::arrange(category, description)

data$description <- factor(data$description, levels = unique(data$description), order = T)
data_commons$description <- factor(data_commons$description, levels = unique(data_commons$description), order = T)

heatmap_fdr_pathway <- heatmap.fdr(data, data_commons) +
  scale_y_discrete(limits = rev) +
  facet_wrap(vars(method)) +
  labs(
    title = "FDR value for each samples in the PDCL datasets in G:Profiler and GSEA",
    x = "Sample"
  )

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/heatmap-fdr-pathway.png", plot = heatmap_fdr_pathway, width = 12, height = 7)

all_commons <- all_enrich %>%  
  dplyr::filter(fdr < 0.1) %>%
  dplyr::count(pathway_id, description, sample) %>%
  dplyr::filter(n>=2) %>%
  dplyr::mutate(category = gs_metadata[pathway_id, "category"]) %>%
  dplyr::mutate(category_merged = merged_categ[category])

count_tcga_commons <- all_commons %>%
  dplyr::filter(sample == "TCGA-GBM") %>%
  dplyr::count(category_merged)
count_pdcl_commons <- all_commons %>%
  dplyr::filter(sample == "PDCL") %>%
  dplyr::count(category_merged)
count_pdcl_personnalized_commons <- all_commons %>%
  dplyr::filter(!(sample %in% c("TCGA-GBM", "PDCL"))) %>%
  dplyr::distinct(pathway_id, description, category_merged) %>%
  dplyr::count(category_merged)

piechart_count_category <- plot_ly() %>% 
  add_pie(data = count_tcga_commons, labels =~category_merged, values =~n, domain = list(row = 0, column = 0), textinfo ="value", 
          textfont = list(size = 24)) %>%
  add_pie(data = count_pdcl_commons, labels =~category_merged, values =~n, domain = list(row = 0, column = 1), textinfo ="value", 
          textfont = list(size = 24)) %>% 
  add_pie(data = count_pdcl_personnalized_commons, labels =~category_merged, values =~n, domain = list(row = 0, column = 2), textinfo ="value", 
          textfont = list(size = 24)) %>%
  layout(
    title = list(
      text = "Count pathways per category for each analysis",
      font = list(size = 32),
      y = 1.0,
      yref = "paper",
      yanchor = "top"
    ),
    legend = list(
      title = list(text="Category"),
      font = list(size = 24),
      y = 0.5,
      yref = "paper",
      yanchor = "middle"
    ),
    grid = list(rows = 1, columns = 3),
    annotations = list(
      list(
        text = "TCGA",
        x = 0.15,
        y = 0.9,
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "top",  
        showarrow = FALSE,
        font = list(size = 24)
      ),
      list(
        text="PDCL",
        x = 0.5,
        y = 0.9,
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "top",  
        showarrow = FALSE,
        font = list(size = 24)
      ),
      list(
        text="PDCL Personnalized",
        x = 0.85,
        y = 0.9,
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "top",  
        showarrow = FALSE,
        font = list(size = 24)
      )
    )
  )
# save_image(piechart_count_category, "/home/spinicck/PhD/paper/analysis-pdcl/manuscript/img/complete-piechart-categ-pdcl.png", width = 2000, height = 800, scale = 2)

piechart_count_category <- plot_ly() %>% 
  add_pie(data = count_tcga_commons, labels =~category_merged, values =~n, domain = list(row = 0, column = 0), textinfo ="value" ) %>%
  add_pie(data = count_pdcl_commons, labels =~category_merged, values =~n, domain = list(row = 0, column = 1), textinfo ="value" ) %>% 
  layout(
    title = list(
      text = "Count pathways per category for each datasets",
      font = list(size = 32),
      y = 1.0,
      yref = "paper",
      yanchor = "top"
    ),
    legend = list(
      title = list(text="Category"),
      font = list(size = 24),
      y = 0.5,
      yref = "paper",
      yanchor = "middle"
    ),
    grid = list(rows = 1, columns = 2),
    annotations = list(
      list(
        text = "TCGA-GBM",
        x = 0.25,
        y = 0.95,
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "top",  
        showarrow = FALSE,
        font = list(size = 24)
      ),
      list(
        text="PDCL",
        x = 0.75,
        y = 0.95,
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "top",  
        showarrow = FALSE,
        font = list(size = 24)
      )
    )
  )
# save_image(piechart_count_category, "/home/spinicck/PhD/paper/analysis-pdcl/manuscript/img/complete-piechart-categ-pdcl.png", width = 2000, height = 1000, scale = 2)