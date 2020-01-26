#' This function takes in a chemical formula and will return all possible isomers based on the
#' KEGG database. Note that this function will access the KEGGREST API instead of the compound library
#' `CompMZ` data frame supplied in with this package.
#' @author Daniel Braas
#' @param formula is the chemical formula entered with quotation marks, ie as string
#' @return a data frame of all potential isomers of that particular chemical formula
#' @export
#' @seealso \code{\link{keggFind}}

isomer_KEGG <- function(formula){
  compounds <- data.frame(KEGGREST::keggFind("compound", formula,'formula'))
  names(compounds)[1] <- 'Formula'
  compounds$KEGG.ID <- gsub('cpd:','', rownames(compounds))
  compounds$KEGG.Entry <- ''
  compounds$Name <- ''
  for (i in 1:nrow(compounds)){
    compounds$KEGG.Entry[i] <- KEGGREST::keggFind('compound', compounds$KEGG.ID[i])
    compounds$Name <- gsub(';(.)*','',compounds$KEGG.Entry)
  }
  select(compounds, Name, KEGG.ID, Formula, KEGG.Entry) %>%
    filter(Formula == formula) %>%
    return()
}

#' This function uses the mapvalues function from plyr without the need to load the package.
#' @param data The data (column) to be modified
#' @param from The original value for what is to be renamed
#' @param to The value something is mapped to
#' @return A vector (within a data table) with re-mapped values
#' @export

mapvalues <- function(data, from, to){
  data <- plyr::mapvalues(data, from = from, to = to)
  return(data)
}
