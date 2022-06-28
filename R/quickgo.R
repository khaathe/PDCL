library(tidyverse)
library(httr)
library(png)
QUICK_GO_BASE_URL <- "https://www.ebi.ac.uk/QuickGO/services/"

# Send a GET request to retrieve the chart
# associated with the GO terms (go.ids).
# The chart is save in a png image file
# defined by the user.
download.go.chart <- function(go.ids, file){
  requestURL <- paste0(
    QUICK_GO_BASE_URL,
    "ontology/go/terms/%7Bids%7D/chart?ids=",
    paste0(go.ids, collapse = ",")
  )
  response <- GET(requestURL, accept("image/png"))
  stop_for_status(response)
  response_content <- content(response, as = "parsed")
  writePNG(response_content, target = file)
}
