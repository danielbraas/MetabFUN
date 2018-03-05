#' This function produces a data frame with mass isotopologue distribution (MID) data that is corrected for naturally occurring 13C.
#' @author Daniel Braas
#' @param data_frame The input data. This should be a data frame with Name, Iso, Condition, Time, Exp, Value, Used_ID, KEGG.ID,
#' Nr.C, Rt and Formula columns. The function assumes that there is a column (Time), which is a factor variable.
#' @param iso The isotopologues that are to be plotted.
#' @return A line graph for the specified metabolites and isotopologues.
#' @export

make_line_graph <- function(data_frame, iso){
  if (exists('Title')==F) stop('Title not specified')
  test <- filter(data3, Name %in% c(glycolysis, PPP))

  test <- test %>%
    mutate(Condition = str_replace(Condition, '-','.')) %>%
    separate(Condition, c('Condition','Time'), sep='-') %>%
    mutate(Condition = str_replace(Condition,'\\.','-'),
           Condition = factor(Condition, levels=unique(Condition)),
           Time = factor(Time, levels = unique(Time)))

  p <- c('M1','M2','M3')
  graph <- ggplot(test)

  for (i in 1:length(p)){
    graph <- graph +
      geom_errorbar(aes(Time, Norm_Av, group=Condition, ymin = Norm_Av - Norm_Std, ymax = Norm_Av + Norm_Std),
                width=.05, color='black', data=filter(test, Iso==p[i]))+
      geom_line(aes(Time, Norm_Av, group=Condition, color=Condition, linetype=Iso), size=1.3, data=filter(test, Iso==p[i]))+
      geom_point(aes(Time, Norm_Av), color='black', data=filter(test, Iso==p[i]))
  }

  graph + facet_wrap(~ Name, scales = 'free')+
    scale_color_manual(name='Condition', values=colors[c(3,6,9,12)])+
    geom_point(aes(Time, Norm_Av), color='black', data=filter(test, Iso=='M2'))+
    theme_bw()+
    theme(text = element_text(size=14, face='bold'))+
    labs(x='Time',y='Percent [C-13] labeled from 1,2-[C-13]Glc')
}
