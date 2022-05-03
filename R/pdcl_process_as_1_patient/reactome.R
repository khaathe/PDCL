library(httr)

REACTOME_BASE_URL <- "https://reactome.org/ContentService/data/"

REACTOME_TOP_LEVEL_PATHWAYS = c(
  "Autophagy",
  "Cell Cycle",
  "Cell-Cell communication",
  "Cellular responses to stimuli",
  "Chromatin organization",
  "Circadian Clock",
  "DNA Repair",
  "DNA Replication",
  "Developmental Biology",
  "Digestion and absorption",
  "Disease",
  "Drug ADME",
  "Extracellular matrix organization",
  "Gene expression (Transcription)",
  "Hemostasis",
  "Immune System",
  "Metabolism",
  "Metabolism of RNA",
  "Metabolism of proteins",
  "Muscle contraction",
  "Neuronal System",
  "Organelle biogenesis and maintenance",
  "Programmed Cell Death",
  "Protein localization",
  "Reproduction",
  "Sensory Perception",
  "Signal Transduction",
  "Transport of small molecules",
  "Vesicle-mediated transport",
  "Unknown"
)

get.reactome.ancestors <- function(react.id){
  query_url <- paste0(REACTOME_BASE_URL, "event/", react.id, "/ancestors")
  response <- GET(query_url)
  stop_for_status(response)
  # Reactome return a list containing one element, which is
  # the list of ancestors. Thus we only take the first and
  # onloy item of the list.
  ancestors <- httr::content(response, type = "application/json")[[1]]
  ancestors
}

get.list.reactome.ancestors <- function(react.ids){
  # Unfortunately, Reactome does not provide a
  # method to send multiple IDs and get their ancestors.
  # Thus have to make a request for each IDs.
  ancestors_list <- lapply(react.ids, function(id){
    ancestors <- get.reactome.ancestors(id)
  })
  names(ancestors_list) <- react.ids
  ancestors_list
}

get.reactome.categories <- function(ancestors.list){
  categories <- sapply(ancestors.list, function(x){
    category <- "Unknown"
    lastAncestor <- x[[length(x)]]
    # The last ancestor is usually the TopLevelPathway
    # or the categorie we want to map our pathway to.
    if ( lastAncestor$schemaClass == "TopLevelPathway"){
      category <- lastAncestor$displayName
    }
    category
  })
  categories
}
