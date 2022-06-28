library(KEGGREST)
library(stringr)

KEGG_TOP_CATEGORIES <- c(
  "Cellular Processes",
  "Drug Development",
  "Environmental Information Processing",
  "Genetic Information Processing",
  "Human Diseases",
  "Metabolism",
  "Organismal Systems",
  "Unkown"
)

get.kegg.entries <- function(kegg.ids){
  entries <- list()
  ids <- kegg.ids
  while ( !isEmpty(ids) ) {
    n <- min(10, length(ids))
    keggres <- keggGet(ids[0:n])
    entries <- c(entries, keggres)
    ids <- ids[0:-n]
  }
  names(entries) <- kegg.ids
  entries
}

get.kegg.categ <- function(kegg.entries){
  categories <- lapply(kegg.entries, function(x){
    str <- x$CLASS
    if(is.null(str)){
      str <- "Unkown"
    }
    str
  })
  categories <- str_split(categories, ";\\s*")
  names(categories) <- names(kegg.entries)
  categories
}
