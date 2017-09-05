#' This function produces a data frame with fractional contribution data.
#' @author Daniel Braas
#' @param FC The input data. This should be a data frame that contains MID data generated with the 'make_MID' function.
#' @return A data frame with fractional contribution data.
#' @export

make_FC <- function(DF){
  if (exists('Title')==F) stop('Title not specified')
  FC <- suppressWarnings(DF %>%
    inner_join(., Abbrev, by=c('Name'='Abb', 'KEGG.ID'='KEGG.ID')) %>%
    select(Name, KEGG.ID, Condition, Iso, Nr.C, starts_with('MID')) %>%
    gather(Exp, MID, -Name, -KEGG.ID, -Condition, -Iso, -Nr.C) %>%
    mutate(i=as.numeric(gsub('M','',.$Iso)),
           iSi=i*MID) %>%
    group_by(Name, Condition, Exp) %>%
    mutate(FC=sum(iSi, na.rm=T)/Nr.C) %>%
    filter(Iso=='M0') %>%
    ungroup() %>%
    mutate(FC=mapvalues(FC, 0, NA)) %>%                #this is important when a sample is missing or all MID values were 0
    select(Name, KEGG.ID, Condition, Exp, FC) %>%
    group_by(Name, Condition) %>%
    mutate(Norm_Av=mean(FC, na.rm=T),
           Norm_Std=sd(FC, na.rm=T),
           CV=Norm_Std/Norm_Av,
           Av = Norm_Av) %>%
    ungroup())

  FC$Exp <- gsub('MID','FC',FC$Exp)
  test1 <- split(FC, FC[c('Name', 'Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$FC)))

  FC <- FC %>%
    arrange(Condition, Name) %>%
    mutate(FC=new.Value)

  data8=split(FC, FC[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$FC~x$Condition))$Pr[1]))
  ANOVA=rep(ANOVA,1,each=length(unique(info$Condition)))

  FC <- spread(FC, Exp, FC) %>%
    arrange(Name) %>%
    mutate(Sig='NA')
  FC$ANOVA <- ANOVA

  for (i in 1:nrow(FC)){
    if (FC$ANOVA[i] == "NaN") FC$Sig[i]=""
    else if (FC$ANOVA[i] <= 0.001) FC$Sig[i]="***"
    else if (FC$ANOVA[i] <= 0.01) FC$Sig[i]="**"
    else if (FC$ANOVA[i] <= 0.05) FC$Sig[i]="*"
    else FC$Sig[i]=""
  }

  write.csv(FC, file=paste0(Title, '-fractional contribution.csv'), row.names=F)
  save(FC, file='FC.rdata')
  return(FC)
}
