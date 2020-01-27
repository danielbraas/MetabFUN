#' This function takes in a chemical formula and will return all isomers found in the
#' compound library supplied in this package.
#' @author Daniel Braas
#' @param formula The value for the Type (e.g. H2O for Type Formula)
#' @return a data frame of all potential isomers of that particular chemical formula
#' @export

isomer <- function(formula = 'H2O'){

  if (value == '') stop('You have to add a formula')

  res <- filter(CompMZ, Formula == formula) %>%
    select(ID:MonoisotopicMass) %>%
    distinct()

  if (nrow(res) < 1){
    return('Entered formula does not exist in current compound library')
  } else {
    return(res)
  }
}

#' This function takes in a chemical formula and will return all possible isomers based on the
#' KEGG database. Note that the function will pull the data from the KEGG API, thus will require
#' an internet connection.
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

#' This function takes in an m/z value and will return all possible isobaric compounds
#' based on the supplied compound data frame supplied in this package.
#' @author Daniel Braas
#' @param MZ The m/z value to be searched
#' @param ppm The m/z deviation allowed to search the compound library
#' @param pol The polarity of the ion to be searched
#' @return a data frame of all potential isomers of that particular chemical formula
#' @export

isobar <- function(MZ = 18.0106, ppm = 5, pol = 'positive'){

  if (tolower(pol) %in% c('positive', 'pos')){
    pol <- 'positive'
  } else {
    pol <- 'negative'
  }

  if (MZ == '')  stop('You have to select an m/z value')

  MZ <- suppressWarnings(as.numeric(MZ))

  if (is.na(MZ)) stop('You have to specify a number for the m/z value')

  ppm_max <- MZ + ppm * MZ / 1e6
  ppm_min <- MZ - ppm * MZ / 1e6

  res <- filter(CompMZ, mz >= ppm_min, mz <= ppm_max, Mode == pol) %>%
    select(mz:Adduct)

  if (nrow(res) < 1) {
    print('m/z value not found in compound library')
  } else{
    return(res)
  }
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
