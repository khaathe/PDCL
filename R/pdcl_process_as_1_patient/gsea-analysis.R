######## Load Library needed
library(fgsea)
library(dplyr)
library(tidyr)

# Run Multiple GSEA analysis and return a list of result per sample for each GMTs
# Set the seed to allow for reproductible results as GSEA do multiple permutation
# to compute its null-distribution.
# Return : list[[GMTs]][[sample]]
run.gsea <- function(all.ranked.list, gmt.list, gmt.desc, seed.num = 1, min.size = 15, max.size = 500, n.perm = 10000, n.proc=8){
  set.seed(seed.num)
  imap(gmt.list, .f = function(gmt, db){
    lapply(all.ranked.list, function(ranked.list){
      res <- fgsea(gmt, ranked.list, minSize=min.size, maxSize=max.size, nperm = n.perm, nproc = n.proc)
      res$description <- sapply(gmt.desc[[db]][res$pathway], function(x){ x })
      res
    })
  })
}

# Convert a list of GSEA Analysis result to a list of Data Frame containing the column in the Generic Enrichment Map (GEM)
# as described in GEM documentation (https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files)
# path.id.to.desc is list with pathway_id as key and the description as value
convert.gsea.to.enrichment <- function(res){
  lapply(res, function(x){
    x %>% dplyr::rename(pathway_id = pathway, p_value= pval, fdr=padj) %>%
      mutate(phenotype = sign(ES), genes = sapply(leadingEdge, function(e){ paste0(e, collapse = ",")} ) ) %>%
      select(pathway_id, description, p_value, fdr, phenotype, genes) %>%
      tidyr::drop_na()
  })
}

# Read a gmt file
# Return a list object indexing using the id of the pathway (1st column e). 
# Each element in the list contains the description of the pathway (2nd
# column) and the list of genes.
read.gmt <- function(gmt.file, delim="\t"){
  gmt.list <- list()
  for (l in strsplit(readLines(con = gmt.file), delim, perl=T) ) {
    gmt.list[[l[1]]] <- list("description"=l[2], "genes"= as.vector(l[-2:-1]) )
  }
  return(gmt.list)
}

# Save a list of Generic Enrichment Map to tab delimited txt files
# Files name conttain the tool used for differential expression analysis, which patient
# has been compared to the control (analysis), the source of the gmt file (db.source)
save.gsea.gem.list.to.txt <- function(all.enrich, save.dir){
  for (db.source in names(all.enrich) ) {
    for (sample.name in names(all.enrich[[db.source]]) ){
      f <- paste0(save.dir, "gsea_", sample.name, "_", db.source, ".txt")
      write.table(all.enrich[[db.source]][[sample.name]], file = f, row.names = F, col.names = T, sep = "\t", quote = F)
    }
  }
}
