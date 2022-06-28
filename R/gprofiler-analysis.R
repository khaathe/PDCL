######## Load Library needed
library(gprofiler2)
library(dplyr)

# Run Multiple G:Profiler analysis and return a list of result per sample for each GMTs
# Return : list[[GMTs]][[sample]]
run.gost <- function(all.gene.list, gmt.list, domain, alpha = 0.05, only_significant = F, adjustment.method = "fdr", return.genes = T){
  lapply(gmt.list, function(gmt){
    lapply(all.gene.list, function(gene.list){
      gost(
        query = gene.list,
        organism = gmt,
        user_threshold = alpha,
        significant = only_significant,
        correction_method = adjustment.method,
        evcodes = return.genes,
        custom_bg = domain
      )  
    })
  })
}

# Convert a list of Gost Analysis result to a list of Data Frame containing the column in the Generic Enrichment Map (GEM)
# as described in G:Profiler documentation (https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis)
convert.gost.to.enrichment <- function(res){
  lapply(res, function(x){
    x$result %>% dplyr::rename(pathway_id = term_id, description = term_name, genes = intersection) %>%
      mutate(fdr = p_value, phenotype = 1) %>%
      select(pathway_id,description, p_value, fdr, phenotype, genes)
  })
}

# Save a list of Generic Enrichment Map to tab delimited txt files
# Files name conttain the tool used for differential expression analysis, which patient
# has been compared to the control (analysis), the source of the gmt file (db.source)
save.gost.gem.list.to.txt <- function(all.enrich, save.dir){
  for (db.source in names(all.enrich) ) {
    for (sample.name in names(all.enrich[[db.source]]) ){
      f <- paste0(save.dir, "gprofiler_", sample.name, "_", db.source, ".txt")
      write.table(all.enrich[[db.source]][[sample.name]], file = f, row.names = F, col.names = T, sep = "\t", quote = F)
    }
  }
}
