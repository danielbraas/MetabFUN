#' This function produces a data frame with fractional contribution data.
#' @author Daniel Braas
#' @param DF The input data. This should be a data frame that contains the MID envelope. The data frame should have the following
#' columns: Used_ID, Iso and each sample separated into individual columns.
#' @param max.iso The maximal number of isotopologues in the data (default is 50).
#' @return A data frame with fractional contribution data.
#' @export

make_unlabeled <- function(DF, max.iso=50){
  
  if (exists('Title')==F) stop('Title not specified')
  DF$Iso <- factor(DF$Iso, levels=paste0('M', 0:max.iso))
  DF <- DF %>% 
    gather(Sample, Value,-Used_ID, -Iso) %>% 
    arrange(Used_ID, Sample, Iso)
  
  flat <- function(dat){
    for (i in 1:nrow(dat)){
      if (is.na(dat$Value[i])) {
        dat$Value[i:nrow(dat)] <- NA
        return(dat)
      }
    }
  }
  
  DF <- DF %>% 
    split(.[c('Used_ID','Sample')]) %>% 
    map(~ flat(.)) %>% 
    do.call(rbind,.)
  
  rownames(DF) <- NULL
  
  DF <- DF %>% 
    spread(Sample, Value)
  
  write.csv(DF, paste(Title, '-Uncorrected raw data.csv'), row.names = F)
  return(DF)
}