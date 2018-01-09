#' This function takes in a data matrix and produces a PCA plot and a correlation circle plot.
#' @author Daniel Braas
#' @param data The data matrix to be used.
#' @param a The PC used for the x-axis.
#' @param b The PC used for the y-axis.
#' @param cutoff The cutoff for correlation.
#' @return A pdf file with a scree plot, a pair plot showing the first five PCs as well as a PCA plot with PCs a and b as specified in the function call and the corresponding top 30 loadings as bar plot.
#' @export

make_PCA2 <- function(matrix, a=1, b=2, cutoff = 0.5){

  if (exists('Title')==F) stop('Title not specified')

  if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Exp') {
    ext = 'Relative Amounts'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='MID') {
    ext = 'MIDs'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='FC') {
    ext = 'Fractional Contribution'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Labeled') {
    ext = 'Percent Labeled'
  }

  pca <- prcomp(t(matrix), center=T, scale=T)
  var_PCs=round(summary(pca)$imp[2,]*100,1)
  CCP <- cor(scale(t(matrix), center=T), pca$x, use='pairwise') %>%
    data.frame() %>%
    mutate(Metabolite = rownames(.))

  CCP$Corr <- sqrt((CCP[,a])^2 + (CCP[,b])^2)

  PC <- pca$x %>%
    as.data.frame() %>%
    mutate(ROWNAMES = rownames(.),
           Sample.Name = gsub('Exp|MID|FC|Labeled','', ROWNAMES)) %>%
    left_join(., samples, by='Sample.Name')
  loadings <- data.frame(pca$rotation)
  loadings$Name <- rownames(loadings)
  loadings <- select(loadings, Name, everything())
  write.csv(loadings, file=paste0(Title,'-Loadings-',ext,'.csv'), row.names=T)
  scores=pca$x
  write.csv(scores, file=paste0(Title,'-Scores-',ext,'.csv'), row.names=T)

  PC.title=paste(Title,'-PCA Plots2-', ext, '.pdf', sep='')
  pdf(file = PC.title, width=16, height=10)

  plot(var_PCs[1:min(10, length(pca$sdev))], type='b', pch=20, col='blue', ylab='Variance explained (%)', xlab='Principal Component', main='Screeplot')
  print(lattice::splom(data=PC, ~PC[,1:5],
                       groups = PC$Condition,
                       par.settings=list(superpose.symbol=list(col=colors, pch=19)),
                       auto.key=list(columns=4), pch=19))

  plot <- ggplot(PC, aes(PC[,a], PC[,b], fill=Condition, label=Sample.Name))+
    geom_text(color='black', fontface='bold', size=4, vjust=-0.2)+
    geom_point(size=5, shape=21, color='black')+
    labs(list(x=paste('PC',a, ' (',var_PCs[a],'%)',sep=''), y=paste('PC',b, ' (',var_PCs[b],'%)',sep=''), title=paste('PC',a,' vs. PC',b,': All samples', sep='')))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          text = element_text(face='bold'))+
    scale_fill_manual('Condition', values=colors)

  CCP_plot <- filter(CCP, Corr >= cutoff) %>%
    ggplot(., aes(.[,a], .[,b], label=Metabolite))+
    geom_point(size=2, color='grey90')+
    geom_text(vjust=-1, color='navy', size=3)+
    labs(list(x=paste('PC',a, ' (',var_PCs[a],'%)',sep=''),
              y=paste('PC',b, ' (',var_PCs[b],'%)',sep=''),
              title = 'Correlation circle plot'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          text = element_text(face='bold'))+
    xlim(-1,1)+
    ylim(-1,1)

  gridExtra::grid.arrange(plot, CCP_plot, nrow=1)
  dev.off()

  CCP1 <- suppressWarnings(CCP %>%
    rename(Norm_Av = PC1) %>%
    right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
    select(KEGG.ID, Norm_Av))
  CCP1$Norm_Av[is.na(CCP1$Norm_Av)] <- 0
  write.csv(CCP1, paste0('CCP-PC1-', ext,'.csv'), row.names=F)
  CCP2 <- suppressWarnings(CCP %>%
                             rename(Norm_Av = PC2) %>%
                             right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
                             select(KEGG.ID, Norm_Av))
  CCP2$Norm_Av[is.na(CCP2$Norm_Av)] <- 0
  write.csv(CCP2, paste0('CCP-PC2-', ext,'.csv'), row.names=F)
  CCP3 <- suppressWarnings(CCP %>%
                             rename(Norm_Av = PC3) %>%
                             right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
                             select(KEGG.ID, Norm_Av))
  CCP3$Norm_Av[is.na(CCP3$Norm_Av)] <- 0
  write.csv(CCP3, paste0('CCP-PC3-', ext,'.csv'), row.names=F)
}
