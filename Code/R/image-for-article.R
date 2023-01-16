library(ggplot2)
library(plotly)
library(tidyverse)
library(reticulate)
reticulate::use_python("/usr/bin/python3")
reticulate::py_run_string("import sys")

pdcl_global_rds <- "~/PhD/tmp/all_enrichment_pdcl_global.rds"
tcga_global_rds <- "~/PhD/tmp/all_enrichment_tcga_global.rds"
pdcl_personnalized <- "~/PhD/tmp/all_enrichment_pdcl_personnalized.rds"
kegg_gs_md_file <- "~/PhD/Project/PDCL/Data/gene-set/kegg_geneset_metadata.csv"
reactome_gs_md_file <- "~/PhD/Project/PDCL/Data/gene-set/reactome_geneset_metadata.csv"

####### Load Data

all_enrichment_pdcl_global <- readRDS(pdcl_global_rds)
all_enrichment_tcga_global <- readRDS(tcga_global_rds)
all_enrichment_pdcl_personnalized <- readRDS(pdcl_personnalized)

all_enrich <- rbind(all_enrichment_pdcl_global, all_enrichment_tcga_global, all_enrichment_pdcl_personnalized)

kegg_gs_metadata <- read.csv(kegg_gs_md_file, as.is = T)
reactome_gs_metadata <- read.csv(reactome_gs_md_file, as.is = T)

gs_metadata <- rbind(kegg_gs_metadata, reactome_gs_metadata)
row.names(gs_metadata) <- gs_metadata$pathway_id

threshold <- 0.05

####### Global Analysis

all_enrichment_global <- rbind(all_enrichment_pdcl_global, all_enrichment_tcga_global)

get.common.pathways <- function(x){
  x %>%
    inner_join(x, by=c("pathway_id", "description", "sample", "db")) %>%
    dplyr::filter(fdr.x<threshold, fdr.y<threshold, method.x == "G:Profiler", method.y == "GSEA") %>%
    dplyr::rename(
      p_value_gprofiler = p_value.x,
      p_value_gsea = p_value.y,
      fdr_gprofiler = fdr.x,
      fdr_gsea = fdr.y, 
      phenotype_gsea = phenotype.y,
      gene_gprofiler = genes.x,
      gene_gsea = genes.y
    ) %>%
    dplyr::select(sample, pathway_id, description, db, p_value_gprofiler, fdr_gprofiler, p_value_gsea, fdr_gsea, phenotype_gsea, gene_gprofiler, gene_gsea)
}

get.signficant.pathways <- function(x, gs, threshold){
  x %>% 
    dplyr::inner_join(gs, by = c("pathway_id", "description")) %>% 
    dplyr::filter(fdr<threshold) %>%
    dplyr::mutate(category = factor(category, levels = sort(unique(gs$category))))
}

common_pathway_pdcl_global <- get.common.pathways(all_enrichment_pdcl_global)
common_pathway_tcga_global <- get.common.pathways(all_enrichment_tcga_global)
common_pathway_global <- rbind(common_pathway_pdcl_global, common_pathway_tcga_global)
common_pdcl_and_tcga_global_pathway_id <- common_pathway_pdcl_global$pathway_id[common_pathway_pdcl_global$pathway_id %in% common_pathway_tcga_global$pathway_id]

####### Barplot

barplot.category <- function(x){
  ggplot(x, aes(category, fill = common)) +
    geom_bar(position = "dodge2") + 
    # scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90) ) +
    labs(
      x = "Biological Categories",
      y = "Count Pathways",
      fill = "Method"
    ) +
    scale_x_discrete(drop = F)
}


kegg_significant_global_pdcl <- get.signficant.pathways(all_enrichment_pdcl_global, kegg_gs_metadata, threshold)
kegg_significant_global_pdcl$sample <- factor(kegg_significant_global_pdcl$sample, levels = c("gbm", "tcga_gbm"), labels = c("PDCL", "TCGA-GBM"))
kegg_significant_global_pdcl$common <- factor(kegg_significant_global_pdcl$method, levels = c("Common", "G:Profiler", "GSEA"))
kegg_significant_global_pdcl$common[kegg_significant_global_pdcl$pathway_id %in% common_pathway_pdcl_global$pathway_id] <- "Common"
barplot_kegg_categ_pdcl_global <- barplot.category(kegg_significant_global_pdcl) + labs(title = "Count of deregulated pathways by samples and Kegg category - PDCL")

reactome_significant_global_pdcl <- get.signficant.pathways(all_enrichment_pdcl_global, reactome_gs_metadata, threshold)
reactome_significant_global_pdcl$sample <- factor(reactome_significant_global_pdcl$sample, levels = c("gbm", "tcga_gbm"), labels = c("PDCL", "TCGA-GBM"))
reactome_significant_global_pdcl$common <- factor(reactome_significant_global_pdcl$method, levels = c("Common", "G:Profiler", "GSEA"))
reactome_significant_global_pdcl$common[reactome_significant_global_pdcl$pathway_id %in% common_pathway_pdcl_global$pathway_id] <- "Common"
barplot_reactome_categ_pdcl_global <- barplot.category(reactome_significant_global_pdcl) + labs(title = "Count of deregulated pathways by samples and Reactome category - PDCL")

####### Barplot Reactome category

kegg_significant_global_tcga <- get.signficant.pathways(all_enrichment_tcga_global, kegg_gs_metadata, threshold)
kegg_significant_global_tcga$sample <- factor(kegg_significant_global_tcga$sample, levels = c("gbm", "tcga_gbm"), labels = c("PDCL", "TCGA-GBM"))
kegg_significant_global_tcga$common <- factor(kegg_significant_global_tcga$method, levels = c("Common", "G:Profiler", "GSEA"))
kegg_significant_global_tcga$common[kegg_significant_global_tcga$pathway_id %in% common_pathway_tcga_global$pathway_id] <- "Common"
barplot_kegg_categ_tcga_global <- barplot.category(kegg_significant_global_tcga) + labs(title = "Count of deregulated pathways by samples and Kegg category - TCGA")

reactome_significant_global_tcga <- get.signficant.pathways(all_enrichment_tcga_global, reactome_gs_metadata, threshold)
reactome_significant_global_tcga$sample <- factor(reactome_significant_global_tcga$sample, levels = c("gbm", "tcga_gbm"), labels = c("PDCL", "TCGA-GBM"))
reactome_significant_global_tcga$common <- factor(reactome_significant_global_tcga$method, levels = c("Common", "G:Profiler", "GSEA"))
reactome_significant_global_tcga$common[reactome_significant_global_tcga$pathway_id %in% common_pathway_tcga_global$pathway_id] <- "Common"
barplot_reactome_categ_tcga_global <- barplot.category(reactome_significant_global_tcga) + labs(title = "Count of deregulated pathways by samples and Reactome category - TCGA")

barplot_categ_global <- cowplot::plot_grid(
  barplot_kegg_categ_pdcl_global,
  barplot_kegg_categ_tcga_global,
  barplot_reactome_categ_pdcl_global,
  barplot_reactome_categ_tcga_global,
  labels = c("A", "B", "C", "D"),
  nrow = 2
)

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/barplot-categ-global.png", plot = barplot_categ_global, width = 20, height = 10)

########## Heatmap FDR

fdr.breaks <- c(0.0, 0.01, 0.025, 0.04, 0.05, 0.1, 0.25, 0.50, 0.75, 1.00)
fdr.labels <- c("0.01", "0.025", "0.04", "0.05", "0.1", "0.25", "0.50", "0.75", "1.0")
heatmap.fdr <- function(x){
  ggplot(x, aes(x = method, y = description, fill = fdr)) +
    geom_tile() +
    facet_wrap(~sample) +
    scale_y_discrete(limits = rev) +
    labs(
      x = "Pathway Enrichment Tool",
      y = "Pathways",
      fill = "FDR"
    ) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_fill_brewer(
      palette = "RdBu",
      na.value = "gray30"
    )
}

########## Heatmap FDR Kegg Selected

kegg_selected_pathways <- c(
  "path:hsa05214",
  "path:hsa03030",
  "path:hsa04510",
  "path:hsa04360",
  "path:hsa04512",
  "path:hsa04724",
  "path:hsa04725",
  "path:hsa00100",
  "path:hsa00190",
  "path:hsa04612",
  "path:hsa04723",
  "path:hsa03008",
  "path:hsa04012",
  "path:hsa04072",
  "path:hsa04392",
  "path:hsa04514",
  "path:hsa04010",
  "path:hsa04012",
  "path:hsa04014",
  "path:hsa04115",
  "path:hsa04150",
  "path:hsa04150",
  "path:hsa04151",
  "path:hsa04330",
  "path:hsa04370"
)

kegg_selected_pathways_global <- all_enrichment_global %>% 
  dplyr::filter(pathway_id %in% kegg_selected_pathways) %>%
  dplyr::mutate(
    sample = factor(sample, labels = c("PDCL", "TCGA-GBM")),
    fdr = cut(fdr, breaks = fdr.breaks, labels = fdr.labels),
    description = sub(" - Homo sapiens \\(human\\)", "", description, perl = T)
  )

heatmap_fdr_kegg_selected_global <- heatmap.fdr(kegg_selected_pathways_global) +
  labs(title = "Heatmap FDR of KEGG pathways for analysis at Population scale")

########## Heatmap FDR Reactome Selected

reactome_selected_pathways <- c(
  "R-HSA-112315",
  "R-HSA-1650814",
  "R-HSA-5576891",
  "R-HSA-8948216",
  "R-HSA-112314",
  "R-HSA-191273",
  "R-HSA-6790901",
  "R-HSA-6794362",
  "R-HSA-109581",
  "R-HSA-1296071",
  "R-HSA-1566948",
  "R-HSA-157579",
  "R-HSA-174417",
  "R-HSA-187037",
  "R-HSA-204998"
)

reactome_selected_pathways_global <- all_enrichment_global %>% 
  dplyr::filter(pathway_id %in% reactome_selected_pathways) %>%
  dplyr::mutate(
    sample = factor(sample, labels = c("PDCL", "TCGA-GBM")),
    fdr = cut(fdr, breaks = fdr.breaks, labels = fdr.labels)
  )

heatmap_fdr_reactome_selected_global <- heatmap.fdr(reactome_selected_pathways_global) +
  labs(title = "Heatmap FDR of Reactome pathways for analysis at Population scale")

####### Combine FDR Heatmap

heatmap_fdr_global <- cowplot::plot_grid(heatmap_fdr_kegg_selected_global, heatmap_fdr_reactome_selected_global, labels = c("A", "B"), nrow = 2)

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/heatmap-fdr-global.png", plot = heatmap_fdr_global, width = 20, height = 10)

######## Personnalized analysis PDCL

################ Barplot Kegg Category

kegg_significant_pdcl_personnalized <- get.signficant.pathways(all_enrichment_pdcl_personnalized, kegg_gs_metadata, threshold)

order_kegg_significant_pdcl <- kegg_significant_pdcl_personnalized %>%
  dplyr::count(sample) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(sample)

kegg_significant_pdcl_personnalized$sample <- factor(kegg_significant_pdcl_personnalized$sample, levels = order_kegg_significant_pdcl, ordered = T)

barplot_kegg_categ_pdcl_personnalized <- ggplot(kegg_significant_pdcl_personnalized, aes(sample, fill = category)) +
  geom_bar() +
  # scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90) ) +
  labs(
    title = "Count of deregulated pathways by samples and Kegg category",
    x = "Sample Name",
    y = "Count Pathways",
    fill = "Biological Category"
  ) +
  scale_x_discrete(drop = F)

################ Barplot Reactome Category

reactome_significant_pdcl_personnalized <- get.signficant.pathways(all_enrichment_pdcl_personnalized, reactome_gs_metadata, threshold)

order_reactome_significant_pdcl <- reactome_significant_pdcl_personnalized %>%
  dplyr::count(sample) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(sample)

reactome_significant_pdcl_personnalized$sample <- factor(reactome_significant_pdcl_personnalized$sample, levels = order_reactome_significant_pdcl, ordered = T)

barplot_reactome_categ_pdcl_personnalized <- ggplot(reactome_significant_pdcl_personnalized, aes(sample, fill = category)) +
  geom_bar() +
  # scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90) ) +
  labs(
    title = "Count of deregulated pathways by samples and Reactome category",
    x = "Sample Name",
    y = "Count Pathways",
    fill = "Biological Category"
  ) +
  scale_x_discrete(drop = F)

################ Combine barplot

barplot_categ_pdcl_personnalized <- cowplot::plot_grid(barplot_kegg_categ_pdcl_personnalized, barplot_reactome_categ_pdcl_personnalized, nrow = 2)

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/barplot-categ-pdcl.png", plot = barplot_categ_pdcl_personnalized, width = 12, height = 7)

################ Frequently Kegg deregulated

freq_dereg_kegg_pdcl_personnalized <- kegg_significant_pdcl_personnalized %>% dplyr::filter(category != "Human Diseases") %>% dplyr::count(pathway_id, description, category) %>% dplyr::arrange(desc(n))

################ Frequently Reactome deregulated

freq_dereg_reactome_pdcl_personnalized <- reactome_significant_pdcl_personnalized %>% dplyr::filter(category != "Disease") %>% dplyr::count(pathway_id, description, category) %>% dplyr::arrange(desc(n))

########### TODO: 

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

data <- all_enrichment_pdcl_personnalized %>%
  dplyr::filter(pathway_id %in% selected_pathways) 

data$description <- sub(" - Homo sapiens \\(human\\)", "", data$description, ignore.case = T, perl = T)
data$description <- str_c( data$description, data$pathway_id, sep = " - ")
data$fdr <- cut(data$fdr, breaks = c(0.0, 0.01, 0.05, 0.1, 0.25, 0.50, 1.00), labels = c("0.01", "0.05", "0.1", "0.25", "0.50", "1.0"))

heatmap_fdr_pathway <- ggplot(
  data = data, 
  aes(x = sample, y = description, fill =  fdr)
) +
  geom_tile() +
  labs(
    title = "Pathways FDR value by PDCL",
    x = "PDCL",
    y = "Pathway",
    fill = "FDR"
  ) +
  scale_fill_brewer(
    palette = "RdBu",
    na.value = "gray30"
  ) +
  #We need to specify this. If not, unpopulated categories do not appear
  scale_x_discrete(drop = F) +
  # By default it is sorted in desc order, we reverse it so it is in asc order
  # scale_y_discrete(limits=rev) +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.background = element_rect(fill = "white")
  )

# ggsave(filename = "~/PhD/paper/analysis-pdcl/manuscript/img/heatmap-fdr-pathway.png", plot = heatmap_fdr_pathway, width = 12, height = 7)
