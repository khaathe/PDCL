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
categories <- all_pathways[match_index]

lapply(categories, function(x){ length(x) })

penda_genes <- row.names(penda_res)
pdcl_dereg_genes <- imap_dfr(penda_res, function(x, y){
  dereg_index <- which(x != 0)
  data.frame(gene = penda_genes[dereg_index], dereg = x[dereg_index], pdcl = y)
})

category_genes <- imap_dfr(categories, function(x, y){
  data.frame(category = y, gene = x)
})

pdcl_dereg_categories <- pdcl_dereg_genes %>% inner_join(category_genes, by = "gene")
table_pdcl_dereg_categories <- table( select(pdcl_dereg_categories, category, pdcl) )

super_deregulated_genes <- pdcl_dereg_genes %>%
  count(gene, dereg) %>%
  filter(n==length(pdcl_names)) %>%
  arrange(gene)
table(super_deregulated_genes$dereg)

list_table_dereg_per_pdcl <- lapply(pdcl_names, function(x){
  addmargins(
    table( pdcl_dereg_categories %>% filter(pdcl == x) %>% select(category, dereg))
  )
})
names(list_table_dereg_per_pdcl) <- pdcl_names

pdcl_rank_per_category <- lapply(comp_deseq2_penda, function(x){
  x <- x %>% inner_join(category_genes, by = "gene") %>% filter(penda != 0)
})

plot_list <- imap(pdcl_rank_per_category, function(x, y){
  ggplot(x, aes(category, log2FoldChange)) + 
    geom_boxplot()+
    labs(title = y, x = "category", y = "log2FoldChange") +
    theme(
      axis.text.x = element_text(angle = 90)
    )
})

summary_pdcl_rank_per_category <- lapply(pdcl_rank_per_category, function(x){
  x %>% 
    group_by(category) %>%
    summarise(
      num_up = sum(penda>0),
      num_down = sum(penda<0),
      mean_up = mean(log2FoldChange[log2FoldChange>0]),
      mean_down = mean(log2FoldChange[log2FoldChange<0]),
      sum_up = sum(log2FoldChange[log2FoldChange>0]),
      sum_down = sum(log2FoldChange[log2FoldChange<0]),
      q1 = quantile( log2FoldChange, probs = 0.25),
      median = quantile( log2FoldChange, probs = 0.5),
      q3 = quantile( log2FoldChange, probs = 0.75)
    )
})


gene_of_interest <- c("HIF1A", "HIF3A", "SLC2A1", "SLC16A1", "PSMA1", "VHL", "LDHA", "LDHB", "LDHC", "LDHD", "PDHA1", "PDHB", "PDK1", "PDK2", "PDK3", "PDK4")

pathways_of_interest <- c(
  "Cellular response to hypoxia",
  "Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha",
  "Pyruvate metabolism",
  "Glycolysis",
  "Citric acid cycle (TCA cycle)"
)

selected_pathway_list <- all_pathways[pathways_of_interest]
selected_pathway_genes <- imap_dfr(selected_pathway_list, function(x, y){
  data.frame(category = y, gene = x)
})
gene_of_interest <- unique(as.character(selected_pathway_genes$gene))
