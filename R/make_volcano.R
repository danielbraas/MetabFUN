#' This function produces a volcano plot using data with two different conditions to be compared.
#' @author Daniel Braas
#' @param data input data frame
#' @param a first condition to evaluate
#' @param b second condition to evaluate
#' @return a data frame of all potential isomers of that particular chemical formula
#' @export
make_volcano <- function(data, a, b){
  conditions <- levels(data$Condition)[c(a, b)]
  Name <- ''
  if (sum(grepl('Exp', names(data))) > 0) {
    Name <- 'Relative Amounts - '
    log.val <- 1.5
  }
  if (sum(grepl('FC', names(data))) > 0) {
    Name <- 'Fractional Contribution - '
    log.val <- 1.1
  }
  volc_Title <- paste0('Comparison - ',Name, conditions[1],' over ',conditions[2],'.csv')
  volc_data <- data %>%
    select(Name, Condition, grep('Exp|MID|FC', names(data))) %>%
    filter(Condition %in% conditions) %>%
    mutate(Condition=factor(Condition, levels=conditions)) %>% 
    gather(Exp, Value, -Name, -Condition) %>%
    group_by(Name) %>%
    mutate(Sum=sum(Value, na.rm=T)) %>%
    filter(Sum != 0) %>%
    mutate(Sum = NULL) %>%
    ungroup()
    
  volc_data$Value[volc_data$Value==0] <- NA
  volc_data <- volc_data %>%
    group_by(Name, Condition) %>%
    arrange(Name, Condition) %>% 
    mutate(Av = mean(Value, na.rm=T)) %>% 
    spread(Exp, Value) %>%
    group_by(Name) %>%
    mutate(Ratio = Av[1]/Av[2]) %>% 
    ungroup() %>%
    gather(Exp, Value, -Name, -Condition, -Av, -Ratio)

  volc_data$Ratio <- as.numeric(mapvalues(volc_data$Ratio, from=c('Inf'), to=c('NA')))
  test1 <- split(volc_data, volc_data[c('Name','Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Value)))
  volc_data <- volc_data %>%
    arrange(Condition, Name) %>%
    mutate(Value=new.Value) 
  data8=split(volc_data, volc_data[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Value~x$Condition))$Pr[1]))
  
  volc_data <- volc_data %>%
    select(Name, Ratio) %>%
    distinct() %>%
    mutate(p.val=ANOVA)

  volc_data$Sig <- 'not significant'
  volc_data$Label <- ''
  for (i in 1:nrow(volc_data)){
    if (log(volc_data$Ratio[i]) == '-Inf')  volc_data$Ratio[i]=1
    else if (log(volc_data$Ratio[i]) == 'NA')  volc_data$Ratio[i]=1
    else if (log(volc_data$Ratio[i]) == 'NaN')  volc_data$Ratio[i]=1
  }
  for (i in 1:nrow(volc_data)){
    if (abs(log(volc_data$Ratio[i], 2)) >= log(log.val, 2) & -log(volc_data$p.val[i], 10) > -log(0.05,10)){
      volc_data$Sig[i] = 'significant'
      volc_data$Label[i] = volc_data$Name[i]
    }
  }

  write.csv(volc_data, volc_Title, row.names=F)
  xlabs <- paste0('Ratio: ', Name,conditions[1],' / ',conditions[2], ' (log2)')
  
  ggplot(volc_data, aes(log2(Ratio), -log(p.val,10), color=Sig, label=Label)) +
    geom_point(size=2) +
    geom_hline(yintercept=-log(0.05,10), linetype=2) +
    geom_vline(xintercept=0, linetype=2) +
    theme_bw() +
    labs(list(x=xlabs, y = '- log10 (p-value)', title=gsub('.csv','',volc_Title))) +
    theme(
      plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
      axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
      axis.text=element_text(size=11, face="bold"),
      legend.title=element_text(face="bold", size=12),
      legend.text=element_text(face="bold",size=12),                  #sets legend text
      strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
      panel.grid.major=element_blank()) +
    annotate("rect", xmin = -log(log.val, 2), xmax = log(log.val, 2), alpha = .2,
             ymin = 0, ymax = max(-log(volc_data$p.val, 10))) +
    scale_color_manual('Sig', values = c('grey80',colors[a])) +
    geom_text(vjust=-0.5, color='black')
  ggsave(gsub('.csv','.pdf', volc_Title), height=10, width=14)
}
