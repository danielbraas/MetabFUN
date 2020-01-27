#' Create a tab-delimited file with MDVs and a csv file with FC or \%labeled for INCA
#' @author Daniel Braas
#' @param iso_data The input which contains the isotopologue data
#' @param QE A set of metabolites to be returned
#' @param label_data Either percent labeled or fractional contribution
#' @return Two csv files with a) isotopologue information for INCA and b) a data frame
#' to be used with Cytoscape
#' @export

make_INCA <- function(iso_data, QE, label_data){

  if (dir.exists('INCA_data') == F) dir.create('INCA_data')
  test <- iso_data %>%
    filter(Name %in% QE) %>%
    mutate(Name = factor(Name, levels = QE),
           Av = as.character(round(Norm_Av / 100, 4))) %>%
    select(Name, Condition, Iso, Av) %>%
    spread(Iso, Av) %>%
    arrange(Condition)

  for(i in 1:nrow(test)){
    for(j in 3:length(test)){
      if (is.na(test[i,j]) == T) test[i,j] <- 99
      else if (test[i,j] == 'NaN') test[i,j] <- 0
    }
  }

  test[,3:length(test)] <- apply(test[,3:length(test)], 2, as.numeric)

  for (i in 1:length(test$Condition)){
    filter(test, Condition == levels(test$Condition)[i]) %>%
      write.csv(., paste0('INCA_data/Matlab-13C corrected-INCA-',levels(test$Condition)[i],'.csv'), row.names=F)
  }

  label_data %>%
    select(KEGG.ID, Condition, Norm_Av) %>%
    spread(Condition, Norm_Av) %>%
    write.csv(., paste0('INCA_data/Labeling-flux map', Title, '.csv'), na='NaN', row.names=F)
}
