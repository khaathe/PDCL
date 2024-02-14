library(ggplot2)
library(plotly)
library(tidyverse)
library(viridis)
library(reticulate)
library(pheatmap)

source("~/PhD/Project/PDCL/Code/R/process-result.R")

####### Theming

base_theme <- theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    panel.background = element_rect(fill = "white"),
    panel.ontop = F
)

scatter_theme <- base_theme + theme(
    panel.grid.major = element_line(colour = "grey50", linetype = "solid", linewidth = 0.1),
)

bar_theme <- base_theme + theme(
    panel.grid.major.y = element_line(colour = "grey50", linetype = "solid", linewidth = 0.1)
)

####### Path to Data

pdcl_global_rds <- "~/PhD/Project/PDCL/Result/all_enrichment_pdcl_global.rds"
tcga_global_rds <- "~/PhD/Project/PDCL/Result/all_enrichment_tcga_global.rds"
pdcl_personnalized <- "~/PhD/Project/PDCL/Result/all_enrichment_pdcl_personnalized.rds"
tcga_personnalized <- "~/PhD/Project/PDCL/Result/TCGA-GBM_Personnalized/Result_2_tcga_enrichment_personnalized/all_enrichment_1c40087720a184b034f39ff6b6f5b872.rds"
kegg_gs_md_file <- "~/PhD/Project/PDCL/Data/gene-set/kegg_geneset_metadata.csv"
reactome_gs_md_file <- "~/PhD/Project/PDCL/Data/gene-set/reactome_geneset_metadata.csv"

####### Load Data

all_enrichment_pdcl_global <- readRDS(pdcl_global_rds)
all_enrichment_tcga_global <- readRDS(tcga_global_rds)
all_enrichment_pdcl_personnalized <- readRDS(pdcl_personnalized)
all_enrichment_tcga_personnalized <- readRDS(tcga_personnalized)

kegg_gs_metadata <- read.csv(kegg_gs_md_file, as.is = T)
reactome_gs_metadata <- read.csv(reactome_gs_md_file, as.is = T)
gs_metadata <- rbind(kegg_gs_metadata, reactome_gs_metadata)
row.names(gs_metadata) <- gs_metadata$pathway_id

all_enrichment_global <- rbind(all_enrichment_pdcl_global, all_enrichment_tcga_global)
summary_dereg <- readRDS("~/PhD/Project/PDCL/Data/tcga-dataset/penda_dereg_freq_tcga_gbm_project_only.rds")
pdcl_sample_names <- c( "4339-p21", "4371-p37", "5706-p14", "6190-p43", "6240-p12", "7015-p17", "7060-p18",  "7142-p14", "N13-1300", "N13-1520-p9", "N14-0072", "N14-0870", "N14-1208", "N14-1525", "N15_0460", "N15_0516", "N15-0385", "N15-0661", "N16_0535", "N16-0240")
threshold <- 0.05

####### Plotting functions

str_width <- 40 # Length for pathway desc in heatmap plots
fdr.breaks <- c(0.0, 0.01, 0.025, 0.04, 0.05, 0.1, 0.25, 0.50, 0.75, 1.00) # cut-off for FDR bins
fdr.labels <- c("0.01", "0.025", "0.04", "0.05", "0.1", "0.25", "0.50", "0.75", "1.0") # Lables for FDR bins
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

heatmap_pathways <- function(data){
    df <- data %>%
        dplyr::mutate(
            fdr = cut(fdr, breaks = fdr.breaks, labels = fdr.labels)
        )
    m <- df %>%
        dplyr::mutate(
            fdr = fdr.breaks[as.numeric(fdr)+1] # Add +1 to avoid 0 in fdr.breaks
        ) %>%
        dplyr::select(description, sample, fdr) %>%
        tidyr::pivot_wider(names_from = "sample", values_from = "fdr") %>%
        tibble::column_to_rownames("description") %>%
        as.matrix()
    
    row_clust <- hclust(dist(m))
    roworder <- row_clust$labels[row_clust$order]
    col_clust <- hclust(dist(t(m)))
    colorder <- col_clust$labels[col_clust$order]
    df$description <- factor(df$description, levels = roworder)
    df$sample <- factor(df$sample, levels = colorder)
    ggplot(
        data = df,
        aes(x = sample, y = description, fill =  fdr)
    ) +
        geom_tile() +
        labs(
            x = "Samples",
            y = "Pathway",
            fill = "FDR"
        ) +
        scale_fill_brewer(
            palette = "RdBu",
            na.value = "gray30"
        ) +
        base_theme
}

get.signficant.pathways <- function(x, gs, threshold){
    x %>% 
        dplyr::inner_join(gs, by = c("pathway_id", "description")) %>% 
        dplyr::filter(fdr<threshold) %>%
        dplyr::mutate(category = factor(category, levels = sort(unique(gs$category))))
}

####### Global Analysis

############## Barplot Global
common_pathway_global <- all_enrichment_global %>%
    tidyr::pivot_wider(names_from = "method", values_from = c("p_value", "fdr", "phenotype", "genes")) %>%
    dplyr::mutate(
        common = dplyr::case_when(
            (fdr_GSEA < threshold) & (`fdr_G:Profiler` < threshold) ~ "Common",
            `fdr_G:Profiler` < threshold ~ "G:Profiler",
            fdr_GSEA < threshold ~ "GSEA"
        )
    )

all_enrichment_global <- all_enrichment_global %>%
    dplyr::inner_join(common_pathway_global) %>%
    dplyr::select(pathway_id, description, p_value, fdr, phenotype, genes, sample, db, method, common) %>%
    dplyr::inner_join(gs_metadata) %>%
    dplyr::mutate(
        sample = factor(sample, levels = c("gbm", "tcga_gbm"), labels = c("PDCL", "TCGA-GBM"))
    )

data <- all_enrichment_global %>%
    dplyr::filter(db == "kegg", fdr<threshold) %>%
    dplyr::mutate(category = factor(category, levels = sort(unique(kegg_gs_metadata$category))))
barplot_categ_global_kegg <- ggplot(data, aes(category, fill = common)) +
    geom_bar(position = "dodge2") + 
    # scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90) ) +
    labs(
        title = "Enriched KEGG Categories",
        x = "Biological Categories",
        y = "Count Pathways",
        fill = "Method"
    ) +
    facet_grid(rows = vars(sample), scales = "free") +
    scale_x_discrete(drop = F) +
    bar_theme +
    theme(
        panel.grid.major = element_line(colour = "grey50", linetype = "solid", linewidth = 0.1)
    )
# ggsave(filename = "~/PhD/Article/63122c5e068151b939801c68/img/barplot-categ-global-kegg.png", plot = barplot_categ_global_kegg, width = 19, height = 14, units = "cm")

data <- all_enrichment_global %>%
    dplyr::filter(db == "reactome", fdr<threshold) %>%
    dplyr::mutate(category = factor(category, levels = sort(unique(reactome_gs_metadata$category))))
barplot_categ_global_reactome <- ggplot(data, aes(category, fill = common)) +
    geom_bar(position = "dodge2") + 
    # scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90) ) +
    labs(
        title = "Enriched Reactome Categories",
        x = "Biological Categories",
        y = "Count Pathways",
        fill = "Method"
    ) +
    facet_grid(rows = vars(sample), scales = "free") +
    scale_x_discrete(drop = F) +
    bar_theme +
    theme(
        panel.grid.major = element_line(colour = "grey50", linetype = "solid", linewidth = 0.1)
    )
# ggsave(filename = "~/PhD/Article/63122c5e068151b939801c68/img/barplot-categ-global-reactome.png", plot = barplot_categ_global_reactome, width = 19, height = 14, units = "cm")

########## Heatmap FDR Global

selected_pathways <- c(
  # KEGG
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
  "path:hsa04370",
  # Reactome
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

selected_pathways_global <- all_enrichment_global %>% 
  dplyr::filter(pathway_id %in% selected_pathways) %>%
  dplyr::mutate(
    sample = factor(sample, labels = c("PDCL", "TCGA-GBM")),
    fdr = cut(fdr, breaks = fdr.breaks, labels = fdr.labels),
    db = factor(db, levels = c("kegg", "reactome"), labels = c("KEGG", "Reactome")),
    description = sub(" - Homo sapiens \\(human\\)", "", description, perl = T),
    description = stringr::str_wrap(description, width=40)
  )

heatmap_fdr_selected_global <- ggplot(selected_pathways_global, aes(x = method, y = description, fill = fdr)) +
    geom_tile() +
    facet_grid(rows = vars(db), cols = vars(sample), scales = "free") +
    scale_y_discrete(limits = rev) +
    labs(
        title = stringr::str_wrap("Heatmap FDR of selected pathways for analysis at Population scale", width = 40),
        x = "Pathway Enrichment Tool",
        y = "Pathways",
        fill = "FDR"
    ) +
    scatter_theme +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_fill_brewer(
        palette = "RdBu",
        na.value = "gray30"
    )

# ggsave(filename = "~/PhD/Article/63122c5e068151b939801c68/img/heatmap-fdr-global.png", plot = heatmap_fdr_selected_global, width = 19, height = 32, units = "cm")

######## PDCL

################ KEGG
kegg_significant_pdcl_personnalized <- get.signficant.pathways(all_enrichment_pdcl_personnalized, kegg_gs_metadata, threshold) %>%
    dplyr::mutate(sample = factor(sample, levels = pdcl_sample_names))

order_kegg_significant_pdcl <- kegg_significant_pdcl_personnalized %>%
  dplyr::count(sample, sort = T, .drop = F) %>%
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
  scale_x_discrete(drop = F) +
  bar_theme

################ Reactome
reactome_significant_pdcl_personnalized <- get.signficant.pathways(all_enrichment_pdcl_personnalized, reactome_gs_metadata, threshold) %>%
    dplyr::mutate(sample = factor(sample, levels = pdcl_sample_names))

order_reactome_significant_pdcl <- reactome_significant_pdcl_personnalized %>%
  dplyr::count(sample, sort = T, .drop = F) %>%
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
  scale_x_discrete(drop = F) +
  bar_theme

################ Combine barplot
barplot_categ_pdcl_personnalized <- cowplot::plot_grid(barplot_kegg_categ_pdcl_personnalized, barplot_reactome_categ_pdcl_personnalized, nrow = 2, labels = "AUTO")
# ggsave(filename = "~/PhD/Article/63122c5e068151b939801c68/img/barplot-categ-pdcl.png", plot = barplot_categ_pdcl_personnalized, width = 32, height = 26, units = "cm")

######## TCGA-GBM

################ KEGG
kegg_significant_tcga_personnalized <- get.signficant.pathways(all_enrichment_tcga_personnalized, kegg_gs_metadata, threshold)

order_kegg_significant_pdcl <- kegg_significant_tcga_personnalized %>%
  dplyr::count(sample) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(sample)

kegg_significant_tcga_personnalized$sample <- factor(kegg_significant_tcga_personnalized$sample, levels = order_kegg_significant_pdcl, ordered = T)

barplot_kegg_categ_tcga_personnalized <- ggplot(kegg_significant_tcga_personnalized, aes(sample, fill = category)) +
  geom_bar() +
  # scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_blank()) +
  labs(
    title = "Count of deregulated pathways by samples and Kegg category",
    x = "Sample",
    y = "Count Pathways",
    fill = "Biological Category"
  ) +
  scale_x_discrete(drop = F) +
  bar_theme  +
  theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
  )

################ Reactome
reactome_significant_tcga_personnalized <- get.signficant.pathways(all_enrichment_tcga_personnalized, reactome_gs_metadata, threshold)

order_reactome_significant_tcga <- reactome_significant_tcga_personnalized %>%
  dplyr::count(sample) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(sample)

reactome_significant_tcga_personnalized$sample <- factor(reactome_significant_tcga_personnalized$sample, levels = order_reactome_significant_tcga, ordered = T)

barplot_reactome_categ_tcga_personnalized <- ggplot(reactome_significant_tcga_personnalized, aes(sample, fill = category)) +
  geom_bar() +
  # scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_blank() ) +
  labs(
    title = "Count of deregulated pathways by samples and Reactome category",
    x = "Sample",
    y = "Count Pathways",
    fill = "Biological Category"
  ) +
  scale_x_discrete(drop = F) +
  bar_theme  +
  theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
  )

################ Combine barplot
barplot_categ_tcga_personnalized <- cowplot::plot_grid(barplot_kegg_categ_tcga_personnalized, barplot_reactome_categ_tcga_personnalized, nrow = 2, labels = "AUTO")
# ggsave(filename = "~/PhD/Article/63122c5e068151b939801c68/img/barplot-categ-tcga.png", plot = barplot_categ_tcga_personnalized, width = 32, height = 26, units = "cm")

####### Frequently deregulated Pathway

############## PDCL
id_freq_dereg_pdcl <-
    all_enrichment_pdcl_personnalized %>%
    dplyr::inner_join(gs_metadata) %>%
    dplyr::filter(stringr::str_detect(category, stringr::regex("disease", ignore_case = T), negate = T), fdr<threshold) %>%
    dplyr::count(pathway_id, description, category, db) %>%
    dplyr::arrange(db, desc(n)) %>%
    dplyr::group_by(db) %>%
    dplyr::slice(1:15) %>%
    dplyr::pull(pathway_id)
data <- all_enrichment_pdcl_personnalized %>%
    dplyr::filter(pathway_id %in% id_freq_dereg_pdcl, db == "kegg") %>%
    dplyr::mutate(
        description = sub(" - Homo sapiens \\(human\\)", "", description, ignore.case = T, perl = T),
        description = stringr::str_wrap(description, width = str_width)
    )
heatmap_kegg_pathway_pdcl <- heatmap_pathways(data) +
    labs(
       title = "Heatmap frequently deregulated pathways in PDCL - KEGG"
    ) +
    theme(
        axis.text.x = element_text(angle = 90)
    )

data <- all_enrichment_pdcl_personnalized %>%
    dplyr::filter(pathway_id %in% id_freq_dereg_pdcl, db == "reactome") %>%
    dplyr::mutate(
        description = stringr::str_wrap(description, width = str_width)
    )
heatmap_reactome_pathway_pdcl <- heatmap_pathways(data) +
    labs(
       title = "Heatmap frequently deregulated pathways in PDCL - Reactome" 
    ) +
    theme(
        axis.text.x = element_text(angle = 90)
    )

heatmap_pathway_pdcl <- cowplot::plot_grid(heatmap_kegg_pathway_pdcl, heatmap_reactome_pathway_pdcl, labels = "AUTO", nrow = 2)
# Utiliser sinon la fonction d'export de RStudio qui donne de bons résultats
# ggsave("~/PhD/Article/63122c5e068151b939801c68/img/heatmap-pathways-pdcl.png", heatmap_pathway_pdcl, width = 20, height = 30, units = "cm")
# png("~/PhD/Article/63122c5e068151b939801c68/img/heatmap-pathways-pdcl.png", width = 1000, height = 1000, dpi = 125); heatmap_pathway_pdcl; dev.off()

############## TCGA-GBM
id_freq_dereg_tcga <-
    all_enrichment_tcga_personnalized %>%
    dplyr::inner_join(gs_metadata) %>%
    dplyr::filter(stringr::str_detect(category, stringr::regex("disease", ignore_case = T), negate = T), fdr<threshold) %>%
    dplyr::count(pathway_id, description, category, db) %>%
    dplyr::arrange(db, desc(n)) %>%
    dplyr::group_by(db) %>%
    dplyr::slice(1:15) %>%
    dplyr::pull(pathway_id)
data <- all_enrichment_tcga_personnalized %>%
    dplyr::mutate(
        description = sub(" - Homo sapiens \\(human\\)", "", description, ignore.case = T, perl = T),
        description = stringr::str_wrap(description, width = str_width)
    ) %>%
    dplyr::filter(pathway_id %in% id_freq_dereg_tcga, db == "kegg")
heatmap_kegg_pathway_tcga <- heatmap_pathways(data) +
    labs(
        title = "Heatmap frequently deregulated pathways in TCGA - KEGG" 
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

data <- all_enrichment_tcga_personnalized %>%
    dplyr::filter(pathway_id %in% id_freq_dereg_tcga, db == "reactome") %>%
    dplyr::mutate(
        description = stringr::str_wrap(description, width = str_width)
    )
heatmap_reactome_pathway_tcga <- heatmap_pathways(data) +
    labs(
        title = "Heatmap frequently deregulated pathways in TCGA - Reactome" 
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )
heatmap_pathway_tcga <- cowplot::plot_grid(heatmap_kegg_pathway_tcga, heatmap_reactome_pathway_tcga, labels = "AUTO", nrow = 2)
# Utiliser sinon la fonction d'export de RStudio qui donne de bons résultats
# ggsave("~/PhD/Article/63122c5e068151b939801c68/img/heatmap-pathways-tcga.png", heatmap_pathway_tcga, width = 20, height = 30, units = "cm")
# png("~/PhD/Article/63122c5e068151b939801c68/img/heatmap-pathways-tcga.png", width = 1000, height = 1000, dpi = 125); heatmap_pathway_tcga; dev.off()

####### Heatmap of the contribution

########### For Kegg pathways
kegg_contrib_matrix <- get.genes.in.pathway(kegg_significant_tcga_personnalized) %>%
  dplyr::count(pathway_id, gene) %>%
  dplyr::rename(rank = n) %>%
  tidyr::complete(pathway_id, gene, fill = list(rank = 0)) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rank) %>%
  tibble::column_to_rownames("pathway_id") %>%
  as.matrix()

mu_contrib_pathways <- rowSums(kegg_contrib_matrix)
mu_contrib_genes <- colSums(kegg_contrib_matrix)
kegg_contrib_matrix <- kegg_contrib_matrix[
  mu_contrib_pathways>quantile(mu_contrib_pathways, probs = 0.75), 
  mu_contrib_genes>quantile(mu_contrib_genes, probs = 0.75)
]

kegg_annot_row <- data.frame(pathway_id = row.names(kegg_contrib_matrix)) %>%
  dplyr::inner_join(kegg_gs_metadata) %>%
  tibble::column_to_rownames(var = "pathway_id") %>%
  dplyr::select(category) %>%
  dplyr::rename(Category = category)

kegg_annot_col <- data.frame(gene = colnames(kegg_contrib_matrix)) %>%
  dplyr::inner_join(summary_dereg) %>%
  dplyr::mutate(dereg = case_when(
    up>0.9 ~ "Frequently Up",
    down>0.9 ~ "Frequently Down",
    up+down<0.1 ~ "Frequently Conserved",
    .default = "None"
  )) %>%
  tibble::column_to_rownames(var = "gene") %>%
  dplyr::select(dereg) %>%
  dplyr::rename(Deregulation = dereg)

kegg_annot_colors <- list(
  Deregulation = setNames(c("#E41A1C", "#377EB8", "#4DAC26", "#D9D9D9"), c("Frequently Up", "Frequently Down", "Frequently Conserved", "None")),
  Category  = setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"), unique(kegg_annot_row$Category))
)

pheatmap::pheatmap(
  kegg_contrib_matrix,
  show_rownames = F,
  show_colnames = F,
  annotation_row = kegg_annot_row,
  annotation_col = kegg_annot_col,
  annotation_colors = kegg_annot_colors,
  legend = F,
  main = "TCGA-GBM - Heatmap Gene Contribution to KEGG Pathways",
  filename = "PhD/Article/analysis-pdcl/manuscript/img/gene_contrib_kegg_tcga.png",
  width = 14,
  height = 8
)

############# For Reactome Pathways
reactome_contrib_matrix <- get.genes.in.pathway(reactome_significant_tcga_personnalized) %>%
  dplyr::count(pathway_id, gene) %>%
  dplyr::rename(rank = n) %>%
  tidyr::complete(pathway_id, gene, fill = list(rank = 0)) %>%
  tidyr::pivot_wider(names_from = gene, values_from = rank) %>%
  tibble::column_to_rownames("pathway_id") %>%
  as.matrix()

mu_contrib_pathways <- rowSums(reactome_contrib_matrix)
mu_contrib_genes <- colSums(reactome_contrib_matrix)
reactome_contrib_matrix <- reactome_contrib_matrix[
  mu_contrib_pathways>quantile(mu_contrib_pathways, probs = 0.9), 
  mu_contrib_genes>quantile(mu_contrib_genes, probs = 0.75)
]

reactome_annot_row <- data.frame(pathway_id = row.names(reactome_contrib_matrix)) %>%
  dplyr::inner_join(reactome_gs_metadata) %>%
  tibble::column_to_rownames(var = "pathway_id") %>%
  dplyr::select(category) %>%
  dplyr::rename(Category = category)

reactome_annot_col <- data.frame(gene = colnames(reactome_contrib_matrix)) %>%
  dplyr::inner_join(summary_dereg) %>%
  dplyr::mutate(dereg = case_when(
    up>0.9 ~ "Frequently Up",
    down>0.9 ~ "Frequently Down",
    up+down<0.1 ~ "Frequently Conserved",
    .default = "None"
  )) %>%
  tibble::column_to_rownames(var = "gene") %>%
  dplyr::select(dereg) %>%
  dplyr::rename(Deregulation = dereg)

reactome_annot_colors <- list(
  Deregulation = setNames(c("#E41A1C", "#377EB8", "#4DAC26", "#D9D9D9"), c("Frequently Up", "Frequently Down", "Frequently Conserved", "None")),
  Category  = setNames(rainbow(n_distinct(reactome_annot_row$Category)), unique(reactome_annot_row$Category))
)

pheatmap::pheatmap(
  reactome_contrib_matrix,
  show_rownames = F,
  show_colnames = F,
  annotation_row = reactome_annot_row,
  annotation_col = reactome_annot_col,
  annotation_colors = reactome_annot_colors,
  legend = F,
  main = "TCGA-GBM - Heatmap Gene Contribution to Reactome Pathways",
  filename = "PhD/Article/analysis-pdcl/manuscript/img/gene_contrib_reactome_tcga.png",
  width = 14,
  height = 8
)

