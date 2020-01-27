#' This function produces a bar plot using data from function bar.
#' @author Daniel Braas
#' @param a is the ggplot variable generated with the bar function
#' @param met is the data frame that contains what metabolites should be plotted
#' @param Title is the title for the plot
#' @param x is the title for the x-axis
#' @param y is the title for the y-axis
#' @param axis.text.x tells ggplot whether or not to label x-axis ticks
#' @param scales fixed or free scales
#' @return a data frame of all potential isomers of that particular chemical formula
#' @export

bar_plot = function(a, met, Title, x, y, axis.text.x, scales){
  a + geom_bar(position="dodge", stat="identity", width=0.9) +
    geom_bar(position="dodge", stat="identity", colour="black", width=0.9) +
    facet_wrap( ~ Name, scales=scales) +
    theme_bw() +
    labs(list(x=x, y=y, title=Title, fill=element_blank())) +
    theme(
      plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
      axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
      axis.text=element_text(size=11, face="bold"),
      #axis.text.x=element_blank(),
      axis.text.x=axis.text.x,
      legend.title=element_text(face="bold", size=12),
      legend.text=element_text(face="bold",size=12),                  #sets legend text
      strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
      panel.grid.major=element_blank()) +
    scale_fill_manual("Conditions", values = colors)  +
    geom_errorbar(aes(ymin=Norm_Av, ymax=Norm_Av+Norm_Std), position=position_dodge(0.9), width=.2)+
    geom_text(vjust=-0.1065, color='darkblue', fontface='bold')
}


#' This function produces a bar plot using data from function bar.
#' @author Daniel Braas
#' @param metabolites Metabolic pathway and what type of data to be plotted, for instance: 'glycolysis' for glycolytic metabolites. Individual metabolites can also be entered as vectors.
#' @param bar.type The type of data to be plotted. Can be 'tot' (relative amounts), 'iso' (isotopologue data), 'lab' (percent labeled) or 'FC' (fractioanl contribution)
#' @param repeats (Only needed in previous versions of ggplot2) Number of conditions to be plotted
#' @return a plot initialized with ggplot
#' @import tidyverse
#' @export

bar <- function(metabolites, bar.type, repeats){
  if (length(metabolites) > 1) {
    ending <- 'select metabolites'
  }
  else if (metabolites == 'glycolysis') {
    metabolites = glycolysis
    ending='glycolytic metabolites'
  }
  else if (metabolites=='TCA') {
    metabolites = TCA
    ending <- 'TCA metabolites'
  }
  else if (metabolites=='PPP') {
    metabolites = PPP
    ending <- 'PPP metabolites'
  }
  else if (metabolites=='Curr') {
    metabolites <- Curr
    ending <- 'Currency metabolites'
  }
  else if (metabolites=='Cys'){
    metabolites <- Cys
    ending <- 'Cysteine metabolites'
  }
  else if (metabolites=='Adenine'){
    metabolites <- Adenine
    ending <- 'Adenosine derivatives'
  }
  else if (metabolites=='Cytosine'){
    metabolites <- Cytosine
    ending <- 'Cytidine derivatives'
  }
  else if (metabolites=='Guanine'){
    metabolites <- Guanine
    ending <- 'Guanine derivatives'
  }
  else if (metabolites=='Thymine'){
    metabolites <- Thymine
    ending <- 'Thymine derivatives'
  }
  else if (metabolites=='Uracil'){
    metabolites <- Uracil
    ending <- 'Uracil derivatives'
  }
  else if (metabolites=='AA'){
    metabolites <- AA
    ending <- 'Amino Acids'
  }
  else if (metabolites=='Hex'){
    metabolites <- Hex
    ending <- 'Hexosamine metabolites'
  }
  else if (metabolites=='FA'){
    metabolites <- FA
    ending <- 'Fatty Acids intermediates'
  }
  else if (metabolites=='Fru'){
    metabolites <- Fru
    ending <- 'Fructose Metabolism'
  }
  else if (metabolites=='CoAs'){
    metabolites <- CoAs
    ending <- 'CoA metabolism'
  }
  else if (metabolites=='Neurotrans'){
    metabolites <- Neurotrans
    ending <- 'Neurotransmitter levels'
  }
  else ending = ''

  if (sum(grepl('MID', names(bar.type))) >= 1){
    met = subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Iso=paste(Iso, Sig, sep='\n'),
                   Sig='') %>%
            mutate(Iso = factor(Iso, levels = paste(rep(paste('M', 0:50, sep=''), each=4),
                                              c('','*','**','***'), sep='\n')))

    Title = paste0("Isotopologue distribution of ",ending)
    x <- 'Isotopologue'
    y <- '% Labeled'
    a <-ggplot(met, aes(Iso, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_text(size=11, face="bold")
    bar_plot(a, met, Title, x, y, axis.text.x, scales='free')
  }
  else if (sum(grepl('Exp', names(bar.type))) >= 1){
    met = subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title = paste0("Relative amounts of ",ending)
    x=''
    y='Relative Amounts'
    a <-ggplot(met, aes(Condition, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot(a, met, Title, x, y, axis.text.x, scales='free')
  }

  else if (sum(grepl('Labeled', names(bar.type))) >= 1){
    met <- subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title <- paste0('Percent labeled in ', ending)
    x <- ''
    y <- '% Labeled'
    a <- ggplot(met, aes(Condition, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot(a, met, Title, x, y, axis.text.x, scales='fixed')
  }
  else if (sum(grepl('FC', names(bar.type))) >= 1){
    met <- subset(bar.type, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title <- paste0('Fractional Contribution to ', ending)
    x <- ''
    y <- '% Fractional Contribution'
    a <- ggplot(met, aes(Condition, Norm_Av, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot(a, met, Title, x, y, axis.text.x, scales='fixed')
  }
}
