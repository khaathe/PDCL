######## Load Library needed
library(ggplot2)
library(ggh4x)
library(tidyr)
library(purrr)

# Generic function turn a list into a data.frame.
# Used by the collapse.gem.list function.
map.list.to.df <- function(x, var_name){
  map2_dfr(x, names(x), function(el, el_name){
    el[,var_name] = el_name
    el
  })
}

# Collapse a list of GEM results into one data.frame with :
# - all the column from a GEM file (https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files)
# - db : Source of pathway informations (ex: Reactome)
# - pdcl : ID of the pdcl used for the enrichment analysis
# - method : the enrichment method that has generated this result (ex: GSEA, G:Profiler) 
collapse.enrichment <- function(all_enrichment){
  all_enrichment_collapsed <- all_enrichment
  
  all_enrichment_collapsed <-  map_depth(.x = all_enrichment_collapsed, .depth = 2, .f = function(x){
    map.list.to.df(x, "sample")
  })
  
  all_enrichment_collapsed <-  map_depth(.x = all_enrichment_collapsed, .depth = 1, .f = function(x){
    map.list.to.df(x, "db")
  })
  
  all_enrichment_collapsed <-  map_depth(.x = all_enrichment_collapsed, .depth = 0, .f = function(x){
    imap(all_enrichment_collapsed, function(x,y){
      if( all(isEmpty(x))){
        stop(y, " result is empty.")
      }
    })
    map.list.to.df(x, "method")
  })
  
  all_enrichment_collapsed$sample <- as.factor(all_enrichment_collapsed$sample)
  all_enrichment_collapsed$db <- as.factor(all_enrichment_collapsed$db)
  all_enrichment_collapsed$method <- as.factor(all_enrichment_collapsed$method)
  
  all_enrichment_collapsed
}

# Plot the number of enriched pathways (pathways with an associated FDR below
# a threshold) of all PDCL for each method present in the data.frame collapsed.by.method.
plot.count.pathways <- function(collapsed.by.method){
  ggplot(
    data = data, 
    aes(x = method, fill = db)
  ) +
    geom_bar(position = "dodge2") +
    labs(
      title = "Number of enriched pathways for each enrichment method",
      x = "Method",
      y = "Counts",
      fill = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}

# Turn a list of rank for GSEA analysis into a data.frame object.
collapse.rnk <- function(rnk){
  rnk <- imap_dfr(rnk, .f = function(rank,pdcl){
    data.frame(gene = names(rank), rank = rank, pdcl = pdcl, row.names = NULL)
  })
  rnk$pdcl <- factor(rnk$pdcl)
  rnk
}

# Get the genes in the leadingEdge column from GEM result.
# Return a data.frame where one line correspond to one gene
# in a pathway for one PDCL. All informations from the original
# data.frame are kept.
get.genes.in.pathway <- function(collapsed.by.method){
  genes.in.common.pathways <- collapsed.by.method %>%
    mutate(genes = strsplit(genes, ",")) %>%
    unnest(genes) %>%
    dplyr::rename(gene = genes)
}

# Plot the genes present in the leadingEdge of enriched pathways coloured
# by the rank submitted.
plot.leading.edge <- function(collapsed.by.method, rnk, nb.max.genes=100){
  data <- get.genes.in.pathway(collapsed.by.method) %>%
    dplyr::mutate(rank = rnk[gene])
  
  nb_total_genes <- n_distinct(data$gene)
  nb_pathways <- n_distinct(data$description)
  keep <- data %>% 
    distinct(gene, rank) %>% 
    arrange(desc(abs(rank))) %>% 
    slice(1:nb.max.genes) %>% 
    pull(gene)
  data <- filter(data, gene %in% keep) %>%
    mutate(gene = factor(gene, levels = keep))
  
  plot_title <- paste0("Top ", nb.max.genes , " genes in the ", nb_pathways ," most enriched pathways")
  plot_subtitle <- paste0("Total number of different genes in pathways leading edge : ", nb_total_genes)
  plot_leading_edge <- ggplot(data, aes(gene, description, fill = rank) ) +
    geom_tile() + 
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Gene Symbol",
      y = "Pathways"
    ) +
    scale_fill_viridis_c(option = "plasma") +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_y_discrete(limits = rev)
  plot_leading_edge
}

plot.number.pathways.per.genes <- function(collapsed.by.method, nb.max.genes=100){
  data <- get.genes.in.pathway(collapsed.by.method) %>%
    dplyr::count(gene) %>%
    arrange(desc(n)) %>%
    mutate(gene=factor(gene, levels=gene))
  nb_genes <- nrow(data)
  data <- data %>% slice(1:nb.max.genes)
  plot_title <- paste0("Top ", nb.max.genes , " genes with the biggest count of enriched pathways")
  plot_subtitle <- paste0("Total number of different genes in pathways leading edge : ", nb_genes)
  plot_count_per_gene <- ggplot(data, aes(gene, n) ) +
    geom_col() + 
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Gene Symbol",
      y = "Pathways"
    ) + 
    theme(axis.text.x = element_text(angle = 90) )
  plot_count_per_gene
}

plot.nes.gsea <- function(gseares){
  data <- gseares %>%
    dplyr::mutate(description = reorder(description, NES))
  
  nb_pathways <- nrow(data)
  
  plot_title <- paste0("NES for the top ", nb_pathways , " most significant pathways")
  
  nes_plot <- ggplot(data, aes(NES, description, fill = padj) ) +
    geom_col() +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    labs(
      title = plot_title,
      x = "NES",
      y = "Pathways"
    )
  nes_plot
}
