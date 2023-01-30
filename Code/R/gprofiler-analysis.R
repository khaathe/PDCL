######## Load Library needed
library(gprofiler2)
library(dplyr)

# Run one multi-query G:Profiler analysis and return a list of result per sample for each GMTs
# Return : list[[GMTs]][[sample]]
run.multi.gost <- function(all.gene.list, gmt.list, domain, alpha = 0.05, only_significant = F, adjustment.method = "fdr", return.genes = T){
  lapply(gmt.list, function(gmt){
    gost(
      query = all.gene.list,
      organism = gmt,
      user_threshold = alpha,
      significant = only_significant,
      correction_method = adjustment.method,
      evcodes = return.genes,
      custom_bg = domain
    )
  })
}

# Run Multiple G:Profiler analysis and return a list of result per sample for each GMTs
# Return : list[[GMTs]][[sample]]
run.parralelized.gost <- function(all.gene.list, gmt.list, domain, alpha = 0.05, only_significant = F, adjustment.method = "fdr", return.genes = T){
  # Create a list to loop over for BiocParallel to avoid creating a nested 
  # for-loop  to parallelize
  iter <- list()
  imap(gmt.list, function(gmt, db){
    imap(all.gene.list, function(gene.list, sample){
      el <- list(gene.list = gene.list, gmt = gmt, db = db, sample = sample)
      # append insert all element of a list at the end
      # We store el in a list to add it as one element of the list.
      iter <<- BiocGenerics::append(iter, list(el))
    })
  })
  
  # Loop over the list iter to run parallelized gost analysis.
  bioc_results <- BiocParallel::bplapply(iter, function(x, domain, alpha, only_significant, adjustment.method, return.genes){
    with(x,{
      res <- list(gostres = NULL, db = db, sample = sample)
      res$gostres <- gost(
        query = gene.list,
        organism = gmt,
        user_threshold = alpha,
        significant = only_significant,
        correction_method = adjustment.method,
        evcodes = return.genes,
        custom_bg = domain
      )
      res
    })
  }, domain, alpha, only_significant, adjustment.method, return.genes)
  
  # Loop over the gost results and recreate the hierarchy used by other
  # functions to process results.
  results <- list()
  lapply(bioc_results, function(b){
    with(b, {
      results[[db]][[sample]] <<- gostres
    })
  })
  
  return(results)
}

# Run Parallelized Multiple G:Profiler analysis and return a list of result per sample for each GMTs
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

# Convert a list of Gost Analysis result to a list of Data Frame containing the column in the Generic Enrichment Map (GEM)
# as described in G:Profiler documentation (https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis)
convert.multi.gost.to.enrichment <- function(res,samples){
  lapply(res, function(x){
    l <- lapply(samples, function(s){
      x$result %>% 
        dplyr::filter(grepl(s, query, perl = TRUE)) %>%
        dplyr::rename(pathway_id = term_id, description = term_name, genes = intersection) %>%
        mutate(fdr = p_value, phenotype = 1) %>%
        select(pathway_id,description, p_value, fdr, phenotype, genes)
    })
    names(l) <- samples
    l
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
