######## Load Library needed
library(gprofiler2)
library(dplyr)

run.gost <- function(genes.symbol, domain, gmt.list){
  lapply(gmt.list, function(gmt){
    gost(
      query = genes.symbol,
      organism = gmt,
      user_threshold = 0.05,
      significant = F,
      correction_method = "fdr",
      evcodes = T,
      custom_bg = domain
    )
  })
}

# Convert a list of Gost Analysis result to a list of Data Frame containing the column in the Generic Enrichment Map (GEM)
# as described in G:Profiler documentation (https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis)
convert.gost.to.gem <- function(res){
  lapply(res, function(x){
    x$result %>% dplyr::rename(pathway_id = term_id, description = term_name, genes = intersection) %>%
      mutate(fdr = p_value, phenotype = 1) %>%
      select(pathway_id,description, p_value, fdr, phenotype, genes)
  })
}

# Convert a list of Gost Analysis result to a list of Data Frame with the same column defined in Jorge Scripts
convert.gost.to.enrichment <- function(res){
  lapply(res, function(x){
    x$result %>% dplyr::rename(pathway_name = term_id,pathway_id = term_name, pathway_size = term_size, p_value_adjusted = p_value) %>%
      select(pathway_id,pathway_name, p_value_adjusted,source, pathway_size, query_size,intersection_size)
  })
}

# Save a list of Generic Enrichment Map to tab delimited txt files
# Files name conttain the tool used for differential expression analysis, which patient
# has been compared to the control (analysis), the source of the gmt file (db.source)
save.gem.list.to.txt <- function(gem.list, save.dir, de.analysis){
  for (analysis in names(gem.list) ) {
    for (db.source in names(gem.list[[analysis]]) ){
      f <- paste0(save.dir, "gprofiler_", de.analysis, "_",analysis, "_", db.source, ".txt")
      write.table(gem.list[[analysis]][[db.source]], file = f, row.names = F, col.names = T, sep = "\t", quote = F)
    }
  }
}
