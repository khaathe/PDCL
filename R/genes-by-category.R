source("~/PhD/PDCL/R/gsea-analysis.R")
source("~/PhD/PDCL/R/gprofiler-analysis.R")
source("~/PhD/PDCL/R/plot-result.R")

all_pathways <- gmtPathways("Data/gene-set/jorge-gmt/reactome_gmt_symbol.gmt")
all_pathway_names <- names(all_pathways)

categories_re <- c(
  "^\\bAutophagy\\b$",
  "^\\bCell Cycle\\b$",
  "^\\bCell-Cell communication\\b$",
  "^\\bCellular responses to external stimuli\\b$",
  "^\\bChromatin organization\\b$",
  "^\\bCircadian Clock\\b$",
  "^\\bDNA Repair\\b$",
  "^\\bDNA Replication\\b$",
  "^\\bExtracellular matrix organization\\b$",
  "^\\bGene expression",
  "^\\bMetabolism\\b$",
  "^\\bMetabolism of proteins\\b$",
  "^\\bMetabolism of RNA\\b$",
  "^\\bOrganelle biogenesis and maintenance\\b$",
  "^\\bProgrammed Cell Death\\b$",
  "^\\bProtein localization\\b$",
  "^\\bVesicle-mediated transport\\b$",
  "^\\bSignal Transduction\\b$"
)

pattern <- paste0(categories_re, collapse = "|")
match_index <- grepl(pattern, all_pathway_names, ignore.case = T, perl = T)
pathways_of_interest <- all_pathways[match_index]

lapply(pathways_of_interest, function(x){ length(x) })

pdcl_dereg_genes <- imap_dfr(penda_res, function(x, y){
  dereg_index <- which(x != 0)
  data.frame(gene = penda_genes[dereg_index], dereg = x[dereg_index], pdcl = y)
})

category_genes <- imap_dfr(pathways_of_interest, function(x, y){
  data.frame(category = y, gene = x)
})

pdcl_dereg_categories <- pdcl_dereg_genes %>% inner_join(category_genes, by = "gene")
table_pdcl_dereg_categories <- table( select(pdcl_dereg_categories, category, pdcl) )

list_table_dereg_per_pdcl <- lapply(pdcl_names, function(x){
  addmargins(
    table( pdcl_dereg_categories %>% filter(pdcl == x) %>% select(category, dereg))
  )
})
names(list_table_dereg_per_pdcl) <- pdcl_names

gene_of_interest <- c("HIF1A", "HIF3A", "SLC2A1", "SLC16A1", "PSMA1", "VHL", "LDHA", "LDHB", "LDHC", "LDHD", "PDHA1", "PDHB", "PDK1", "PDK2", "PDK3", "PDK4")
penda_res[gene_of_interest,]
