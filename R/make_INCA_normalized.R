#' Creates a csv file with flux estimates relative to a control flux vector
#' @details
#' This function takes in the output of INCA flux estimates and normalizes them relative
#' to a chosen condition.
#' @author Daniel Braas
#' @param flux_data A file name for a data table with reaction IDs, Equation and flux estimates
#' @param condition The condition to be normalized to
#' @param label_data Either percent labeled or fractional contribution
#' @return Two data frames with flux information.
#' @export

make_INCA_normalized <- function(flux_data, condition, label_data){
  flux_data %>%
    filter(Name %in% QE) %>%
    mutate(Name = factor(Name, levels = QE),
           Av=round(Norm_Av/100, 4)) %>%
    select(Name, Condition, Iso, Av) %>%
    spread(Iso, Av) %>%
    arrange(Condition) %>%
    write.table(., 'INCA_13C corrected.txt', na='NaN', sep='\t', row.names=F)

  label_data %>%
    select(KEGG.ID, Condition, Norm_Av) %>%
    spread(Condition, Norm_Av) %>%
    write.csv(., 'Label data-INCA.csv', na='NaN', row.names=F)
}

