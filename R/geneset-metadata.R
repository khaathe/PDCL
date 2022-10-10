kegg_gmt <- read.gmt("Data/gene-set/kegg.gmt")
reactome_gmt <- read.gmt("Data/gene-set/reactome_no_space.gmt")

kegg_geneset_metadata <- data.frame(
  pathway_id = names(kegg_gmt), 
  description = sapply(kegg_gmt, function(x){x$description}),
  size = sapply(kegg_gmt, function(x){length(x$genes)}),
  category = rep("Unknown", length(kegg_gmt))
)

reactome_geneset_metadata <- data.frame(
  pathway_id = names(reactome_gmt), 
  description = sapply(reactome_gmt, function(x){x$description}),
  size = sapply(reactome_gmt, function(x){length(x$genes)}),
  category = rep("Unknown", length(reactome_gmt))
)

kegg_entries <- get.kegg.entries(kegg_geneset_metadata$pathway_id)
kegg_categories <- get.kegg.categ(kegg_entries)
kegg_categories <- sapply(kegg_categories, function(x){ paste0(x, collapse = ",") })
kegg_geneset_metadata$category <- kegg_categories

react_ancestors_list <- get.list.reactome.ancestors( reactome_geneset_metadata$pathway_id )
react_categories <- sapply(react_ancestors_list, function(x){
  ancestors <- sapply(x, function(y){
    y$displayName
  })
  ancestors <-rev(ancestors)
  ancestors <- paste0(ancestors, collapse = ",")
  ancestors
})
reactome_geneset_metadata$category <- react_categories

kegg_gs_metadata <- read.csv("Data/gene-set/kegg_geneset_metadata.csv", as.is = T)
kegg_gs_metadata$category <- sapply( strsplit(kegg_gs_metadata$category, ","), function(x){x[1]})
ggplot(kegg_gs_metadata, aes(x = category)) + geom_bar() + theme(axis.text.x = element_text(angle = 90) )

react_gs_metadata <- read.csv("Data/gene-set/reactome_geneset_metadata.csv", as.is = T)
react_gs_metadata$category <- sapply( strsplit(react_gs_metadata$category, ","), function(x){x[1]})
ggplot(react_gs_metadata, aes(x = category)) + geom_bar() + theme(axis.text.x = element_text(angle = 90) )

