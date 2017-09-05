#' This function takes in a data and produces a volcano plot relative to a specified condition level.
#' @author Daniel Braas
#' @param DF The data frame to be used.
#' @param Ctrl The condition that everything should be measured against.
#' @param sig_val The ratio cutoff (default is 1.5).
#' @param p_val The ANOVA p-value. Can also be a regular p-value from a t-test (default is 0.05).
#' @return A pdf file with a volcano plot.
#' @export

make_AOV_volcano <- function(DF, Ctrl, sig_val=1.5, p_val=0.05){
  if (exists('Title')==F) stop('Title not specified')


  DF <- DF %>%
    arrange(Name, Condition) %>%
    group_by(Name) %>%
    mutate(Ratio = Av / Av[Ctrl]) %>%
    ungroup() %>%
    select(Name, KEGG.ID, Condition, Av, Ratio, ANOVA)

  write.csv(DF, paste0(Title, '-Volcano Plot using ANOVA normalized to ', levels(DF$Condition)[Ctrl],'.csv'), row.names=F)
  pdf(paste0(Title, '-Volcano Plot using ANOVA normalized to ', levels(DF$Condition)[Ctrl],'.pdf'), width=14, height=10)
  plot <- ggplot(DF, aes(log2(Ratio), -log(ANOVA, 10), fill=Condition, label=Name))+
    theme_bw()+
    theme(text=element_text(face='bold'))+
    labs(x=paste0('log2(Ratio relative to ',levels(DF$Condition)[Ctrl],')'),
         y='log10(ANOVA p-value)')+
    geom_hline(yintercept=-log(p_val,10), linetype=2) +
    geom_vline(xintercept=log2(sig_val), linetype=2) +
    geom_vline(xintercept=-log2(sig_val), linetype=2) +
    geom_point(size=2, shape=21, alpha=0.3)+
    geom_point(data=filter(DF, abs(log2(Ratio)) >= log2(sig_val), -log(ANOVA, 10) >= -log(p_val,10)),
               shape=21, size=5)+
    scale_fill_manual('Condition', values=colors)+
    geom_text(data=filter(DF, abs(log2(Ratio)) >= log2(sig_val), -log(ANOVA, 10) >= -log(p_val,10)),
              vjust=-0.5, size=4, fontface='bold')
  print(plot)
  dev.off()
}
