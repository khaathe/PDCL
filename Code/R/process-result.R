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
plot.leading.edge <- function(collapsed.by.method, rnk, nb.max.genes=100, keep_genes=NULL){
  data <- get.genes.in.pathway(collapsed.by.method) %>%
    dplyr::mutate(rank = rnk[gene])
  
  nb_total_genes <- n_distinct(data$gene)
  nb_pathways <- n_distinct(data$description)
  if(is.null(keep_genes)){
      keep <- data %>% 
          distinct(gene, rank) %>% 
          arrange(desc(abs(rank))) %>% 
          slice(1:nb.max.genes) %>% 
          pull(gene)
      data <- filter(data, gene %in% keep) %>%
          mutate(gene = factor(gene, levels = keep))
  } else {
      data <- filter(data, gene %in% keep_genes) %>%
          mutate(gene = factor(gene, levels = keep_genes))
  }
  
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
    dplyr::distinct(pathway_id, description, gene) %>%
    dplyr::count(gene) %>%
    arrange(desc(n)) %>%
    mutate(gene=factor(gene, levels=gene))
  nb_genes <- nrow(data)
  data <- data %>% dplyr::slice(1:nb.max.genes)
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

heatmap.fdr <- function(collapsed.by.method, highlight.pathways=NULL){
  data <- fdr.to.factor(collapsed.by.method, 5)
  heatmap_fdr <- ggplot(
    data = data, 
    aes(x = sample, y = description, fill =  fdr)
  ) +
    geom_tile() +
    labs(
      title = "FDR value each pathways by PDCL and method",
      x = "PDCL",
      y = "Pathway",
      fill = "FDR"
    ) +
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
  if( !is.null(highlight.pathways) ){
    heatmap_fdr <- heatmap_fdr +
      geom_tile(
        data = highlight.pathways,
        mapping = aes(x = sample, y = description),
        inherit.aes = F,
        colour = 'black',
        fill = NA,
        size = 0.5
      )
  }
  heatmap_fdr
}

penda.res.sort <- function(penda.result.list, gene.domain, keep.sample = TRUE){
  lapply(penda.result.list, function(p){
    u <- p$up_genes[,keep.sample] # UP penda-matrix with samples filtered
    d <- p$down_genes[,keep.sample] # DOWN penda-matrix with samples filtered
    i <- which( !(gene.domain %in% row.names(u)) )
    g <- gene.domain[i] 
    # Create a matrix of the missing genes
    mg <- matrix(FALSE, nrow = length(g), ncol = ncol(u), dimnames = list(g, colnames(u)) )
    # Add missing genes to u and d
    u <- rbind(u, mg)
    d <- rbind(d, mg)
    
    u <- u[order(match(row.names(u), gene.domain)), ]
    d <- d[order(match(row.names(d), gene.domain)), ]
    
    list(up_genes = u, down_genes = d)
  })
}

get.penda.final.matrix <- function(penda.result.list, gene.domain = NULL, sort = T, ...){
  if (is.null(gene.domain)){
    gene.domain <- c()
    # for(p in penda.result.list) { gene.domain <- c(gene.domain, row.names(p$up_genes)) }
    gene.domain <- sapply(penda.result.list, function(p){ gene.domain <- c(gene.domain, row.names(p$up_genes)) })
    gene.domain <- sort(unique(gene.domain))
  }
  
  po.list <- if (sort) penda.res.sort(penda.result.list, gene.domain = gene.domain, ...) else penda.result.list
  
  up_genes <- po.list[[1]]$up_genes
  down_genes <- po.list[[1]]$down_genes
  for (po in po.list[2:length(po.list)]){
    up_genes <- (up_genes & po$up_genes)
    down_genes <- (down_genes & po$down_gene)
  }
  
  dereg_matrix <- matrix(0, nrow = nrow(up_genes), ncol = ncol(up_genes), dimnames = list(row.names(up_genes), colnames(up_genes)))
  dereg_matrix[up_genes] <- 1
  dereg_matrix[down_genes] <- -1
  
  dereg_matrix
}

plot.contrib.genes <- function(
        data,
        colvar = "pathway_id",
        q.pathways = NULL,
        q.genes = NULL,
        g.pattern = "\\w*",
        ...
){
    # Create the contribution matrix
    # We count the number of time a gene is associated to a pathway.
    # Is is better practice to filter the data to keep only the significant 
    # results. We can also apply a filter on the gene name to keep only the genes
    # we wants.
    contrib_matrix <- get.genes.in.pathway(data) %>%
        dplyr::filter(stringr::str_detect(gene, g.pattern)) %>%
        dplyr::count(.data[[colvar]], gene) %>%
        dplyr::rename(rank = n) %>%
        tidyr::complete(.data[[colvar]], gene, fill = list(rank = 0)) %>%
        tidyr::pivot_wider(names_from = gene, values_from = rank) %>%
        tibble::column_to_rownames(colvar) %>%
        as.matrix()
    # Filter the pathways on quantile to reduce the number of pathways displayed
    # in the heatmap (number of rows)
    if( !is.null(q.pathways) ){
        mu_contrib_pathways <- rowSums(contrib_matrix)
        contrib_matrix <- contrib_matrix[
            mu_contrib_pathways>quantile(mu_contrib_pathways, probs = q.pathways),
        ]
    }
    # Filter the genes on quantile to reduce the number of genes displayed
    # in the heatmap (number of columns)
    if ( !is.null(q.genes) ){
        mu_contrib_genes <- colSums(contrib_matrix)
        contrib_matrix <- contrib_matrix[, 
                                         mu_contrib_genes>quantile(mu_contrib_genes, probs = q.genes)
        ]
    }
    
    pheatmap::pheatmap(
        contrib_matrix,
        ...
    )
    
}
