source("~/PhD/PDCL/R/gsea-analysis.R")
source("~/PhD/PDCL/R/gprofiler-analysis.R")
source("~/PhD/PDCL/R/plot-result.R")

identifier <- "gene_symbol_human"
####################################### GSEA Jorge

ranked_lists <- readRDS(file = "~/PhD/Project/Code_Kevin/GSEA_fgsea/rna_seq_all_results_gsea.rds")
gsea_results <- readRDS("~/PhD/Project/Code_Kevin/GSEA_fgsea/gsea_results.rds")

###### INDICATE IDENTIFIER TO USE
ranking_metric <- "-log10_padj_log2FoldChange"

rank_jorge <- lapply(ranked_lists, function(x){
  x <- x %>% 
    select(all_of(identifier),log2FoldChange,stat,pvalue,padj) %>% 
    rename(Gene = identifier) %>% 
    drop_na() %>%
    mutate(
      `-log10_pvalue_sign_log2FoldChange` = -log10(pvalue)*sign(log2FoldChange),
      `-log10_padj_log2FoldChange` = -log10(padj)*log2FoldChange
    )
  gene_list <- x %>% select(all_of(ranking_metric)) %>% pull()
  names(gene_list) <- x$Gene
  #Sort the list in increasing order (for fgsea)
  gene_list <- sort(gene_list, decreasing = F)
  gene_list
})

####################################### G:Profiler Jorge

significant_results <- readRDS("~/PhD/Project/Code_Kevin/ORA_Gprofiler2/rna_seq_results_005.rds")
gprofiler_results <- readRDS("~/PhD/Project/Code_Kevin/ORA_Gprofiler2/enrichment_rna_reactome_bioplanet_gprofiler2_no_unwanted_categ.rds")
gprofiler_gene_list <- lapply(significant_results, function(x){
  x %>%
    dplyr::select(all_of(identifier)) %>%
    tidyr::drop_na()
})

####################################### Process Result with my functions
gsea_gem <- lapply(gsea_results, convert.gsea.to.gem)
gprofiler_gem <- lapply(gprofiler_results, convert.gost.to.gem)

gem_list <- list("G:Profiler" = gprofiler_gem, "GSEA" = gsea_gem)
# For my method
pdcl_names <- names(rank_jorge)

gem_collapsed <- collapse.gem.list(gem_list)
gem_collapsed_method_vs_method <- collapse.method.vs.method(gem_collapsed)
all_common_pathways <- find.method.common.pathways(gem_collapsed_method_vs_method, c(0.05, 0.1))
gem_collapsed_method_vs_method <- gem_collapsed_method_vs_method %>%
  inner_join(all_common_pathways, by=c("pathway_id", "pdcl", "method_1", "method_2"))
threshold <- 0.1
gprofiler_vs_gsea <- plot_method_vs_method(gem_collapsed_method_vs_method, "G:Profiler", "GSEA", threshold)
count_categories_gprofiler_vs_gsea_0.1 <- plot.count.categories.m.vs.m(gem_collapsed, all_common_pathways, threshold, "G:Profiler", "GSEA")
count_pathways_gprofiler_vs_gsea_0.1 <- plot.count.pathways.m.vs.m(gem_collapsed, all_common_pathways, threshold, "G:Profiler", "GSEA")

### Plot Heatmap FDR 
heatmap_pathways_gprofiler_vs_gsea_0.1 <- heatmap.pathways(gem_collapsed, all_common_pathways, threshold, "G:Profiler", "GSEA")

### Plot heatmap Categories
heatmap_categories_0.1 <- heatmap.categories(gem_collapsed, threshold)


############################################# PouPiDoup
run.gsea <- function(ranked.list, gmt.list){
  lapply(gmt.list, function(gmt){
    fgsea(gmt, ranked.list, minSize=15, maxSize=500, nperm = 1000000, nproc = 8)
  })
}

path_reactome_gmt <- "~/PhD/Project/Code_Kevin/Gmt_files/reactome_gmt_symbol_no_unwanted_categ_no_less_10.gmt"
reactome_gmt_file <- gmtPathways(path_reactome_gmt)
gsea_gmt_list <- list(reactome = reactome_gmt_file)
set.seed(179)
gsea_res_list <- lapply(rank_jorge, run.gsea, gsea_gmt_list)
gsea_gem_2 <- lapply(gsea_res_list, convert.gsea.to.gem)

reactome_gmt_token <- upload_GMT_file(gmtfile = path_reactome_gmt)
gprofiler_gmt_list <- list(reactome = reactome_gmt_token)

all_pdcl_genes <- c()
for (x in gprofiler_gene_list) { all_pdcl_genes <- c(all_pdcl_genes, x$gene_symbol_human) }
all_pdcl_genes <- unique(all_pdcl_genes)

gost_res_list <- lapply(gprofiler_gene_list, function(gene_list){
  res <- run.gost(
    gene_list$gene_symbol_human,
    all_pdcl_genes,
    gprofiler_gmt_list
  )
  return(res)
})

gost_gem_2 <- lapply(gost_res_list, convert.gost.to.gem)
gem_list_2 <- list("G:Profiler" = gost_gem_2, "GSEA" = gsea_gem_2)
gem_collapsed_2 <- collapse.gem.list(gem_list_2)
gem_collapsed_method_vs_method_2 <- collapse.method.vs.method(gem_collapsed_2)
all_common_pathways_2 <- find.method.common.pathways(gem_collapsed_method_vs_method_2, c(0.05, 0.1))
gem_collapsed_method_vs_method_2 <- gem_collapsed_method_vs_method_2 %>%
  inner_join(all_common_pathways_2, by=c("pathway_id", "pdcl", "method_1", "method_2"))
threshold <- 0.1
gprofiler_vs_gsea_2 <- plot_method_vs_method(gem_collapsed_method_vs_method_2, "G:Profiler", "GSEA", threshold)
count_categories_gprofiler_vs_gsea_0.1_2 <- plot.count.categories.m.vs.m(gem_collapsed_2, all_common_pathways_2, threshold, "G:Profiler", "GSEA")
count_pathways_gprofiler_vs_gsea_0.1_2 <- plot.count.pathways.m.vs.m(gem_collapsed_2, all_common_pathways_2, threshold, "G:Profiler", "GSEA")

### Plot Heatmap FDR 
heatmap_pathways_gprofiler_vs_gsea_0.1_2 <- heatmap.pathways(gem_collapsed_2, all_common_pathways_2, threshold, "G:Profiler", "GSEA")

### Plot heatmap Categories
heatmap_categories_0.1_2 <- heatmap.categories(gem_collapsed_2, threshold)

