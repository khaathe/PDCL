######## Load Library needed
library(ggplot2)
library(S4Vectors)
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
  
  nb_pathways <- dplyr::n_distinct(data$pathway)
  
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

plot.enriched.categories <- function(collapsed.by.method, category, common.pathways = NULL){
  data <- collapsed.by.method %>%
    mutate(category = category)
  
  if (!is.null(common.pathways)){
    common <- which(data$pathway_id %in% common.pathways)
    data$method <- factor(data$method, level = c("Common Pathway", levels(data$method) ))
    data$method[common] <- "Common Pathway"
    data <- data %>% dplyr::distinct(pathway_id, category, method)
  }
  
  plot_title <- paste0("Count of enriched pathways per categories")
  
  ggplot(data, aes(category, fill = method)) +
    geom_bar(position = "dodge2") + 
    labs(
      title = plot_title,
      x = "Category",
      y = "Count"
    ) + 
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_x_discrete(drop = F)
}

deep_merge <- function(el1, el2, level = 1, max.level = NULL) {
  is_el1_list <- any(class(el1) == "list")
  is_el2_list <- any(class(el2) == "list")
  is_max_depth <- all(!is.null(max.level), level>=max.level)
  if ( !is_el1_list | !is_el2_list | is_max_depth){
    new_el <- c(el1, el2)
    return(new_el)
  }
  keys <- unique(c(names(el1), names(el2)))
  new_el <- lapply(keys, function(x){
    deep_merge(el1[[x]], el2[[x]], level+1, max.level)
  })
  names(new_el) <- keys
  new_el
}

# Generate a vector of breaks for ggplot graph with n values
# inferior and superior a midpoint.
generate.breaks <- function(min, max, midpoint, n){
  left <- seq(from = min, to = midpoint, length.out = n ) 
  right <- seq(from = midpoint, to = max, length.out = n )[-1]
  c(left, right)
}

# Turn the numeric FDR col into a factor with correct labels. 
fdr.to.factor <- function(x, n){
  break.vect <- generate.breaks(0, 1, threshold, n)
  cut.vect <- cut(x$fdr, breaks = break.vect, labels = break.vect[-1])
  cut.vect <- factor(cut.vect, levels = rev(levels(cut.vect)))
  x$fdr <- cut.vect
  x
}

heatmap.fdr <- function(collapsed.by.method, common.pathways){
  data <- fdr.to.factor(collapsed.by.method, 5)
  heatmap_fdr <- ggplot(
    data = data, 
    aes(x = sample, y = description, fill =  fdr)
  ) +
    geom_tile() +
    geom_tile(
      data = common.pathways,
      mapping = aes(x = sample, y = description),
      inherit.aes = F,
      colour = 'black',
      fill = NA,
      size = 0.5
    ) +
    labs(
      title = "FDR value each pathways by PDCL and method",
      x = "PDCL",
      y = "Pathway",
      fill = "FDR"
    ) +
    facet_wrap(vars(method)) +
    scale_fill_brewer(
      palette = "RdBu",
      na.value = "gray30",
      direction = -1
    ) +
    #We need to specify this. If not, unpopulated categories do not appear
    scale_x_discrete(drop = F) +
    # By default it is sorted in desc order, we reverse it so it is in asc order
    # scale_y_discrete(limits=rev) +
    theme(
      axis.text.x = element_text(angle = 90),
      panel.background = element_rect(fill = "white")
    )
  heatmap_fdr
}