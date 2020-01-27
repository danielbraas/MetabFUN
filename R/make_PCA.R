#' This function takes in a data matrix and produces a PCA plot including the top30 loadings.
#' @author Daniel Braas
#' @param matrix The data matrix to be used.
#' @param a The PC used for the x-axis.
#' @param b The PC used for the y-axis.
#' @return A pdf file with a scree plot, a pair plot showing the first five PCs as well as a PCA plot with PCs a and b as specified in the function call and the corresponding top 30 loadings as bar plot.
#' @export

make_PCA <- function(matrix, a, b){

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

  PC.title=paste(Title,'-PCA Plots-', ext, '.pdf', sep='')
  pdf(file = PC.title, width=16, height=10)

  plot(var_PCs[1:min(10, length(pca$sdev))], type='b', pch=20, col='blue', ylab='Variance explained (%)', xlab='Principal Component', main='Screeplot')
  color=colors
  print(lattice::splom(data=PC, ~PC[,1:5],
                 groups = PC$Condition,
                 par.settings=list(superpose.symbol=list(col=colors, pch=19)),
                 auto.key=list(columns=4), pch=19))


  plot <- ggplot(PC, aes(PC[,a], PC[,b], fill=Condition, label=Sample.Name))+
    geom_text(color='black', fontface='bold', size=4, vjust=-0.2)+
    geom_point(size=5, shape=21, color='black')+
    labs(list(x=paste('PC',a, ' (',var_PCs[a],'%)',sep=''), y=paste('PC',b, ' (',var_PCs[b],'%)',sep=''), title=paste('PC',a,' vs. PC',b,': All samples', sep='')))+
    theme_bw()+
    theme(panel.grid.major=element_blank())+
    scale_fill_manual('Condition', values=colors)

  pan_A <- loadings %>%
    arrange(desc(abs(PC1))) %>%
    .[1:30,] %>%
    mutate(col = PC1 > 0) %>%
    ggplot(., aes(reorder(Name, PC1), PC1))+
    geom_bar(stat='identity', aes(fill=col), color='black', width=.9)+
    theme_bw()+
    theme(text=element_text(face='bold'),
          axis.text.x=element_text(angle=90, vjust=.3, hjust=1, size=9))+
    labs(x='', y='Loadings PC1')+
    scale_fill_manual(values=c("#CCEEFF", "#FFDDDD"), guide=F)

  pan_B <-  loadings %>%
    arrange(desc(abs(PC2))) %>%
    .[1:30,] %>%
    mutate(col = PC2 > 0) %>%
    ggplot(., aes(reorder(Name, PC2), PC2))+
    geom_bar(stat='identity', aes(fill=col), color='black', width=.9)+
    theme_bw()+
    theme(text=element_text(face='bold'),
          axis.text.y=element_text(size=8))+
    labs(x='', y='Loadings PC2')+
    coord_flip()+
    scale_fill_manual(values=c("#CCEEFF", "#FFDDDD"), guide=F)

  gridExtra::grid.arrange(plot, pan_B, pan_A, ncol=2, widths = c(5, 3), heights = c(8, 3))
  dev.off()

  LoadPC1 <- suppressWarnings(loadings %>%
    rename(Norm_Av = PC1) %>%
    right_join(., Abbrev, by=c('Name'='Abb')) %>%
    select(KEGG.ID, Norm_Av))
  LoadPC1$Norm_Av[is.na(LoadPC1$Norm_Av)] <- 0
  write.csv(LoadPC1, paste0('LoadPC1-', ext,'.csv'), row.names=F)
  LoadPC2 <- suppressWarnings(loadings %>%
    rename(Norm_Av = PC2) %>%
    right_join(., Abbrev, by=c('Name'='Abb')) %>%
    select(KEGG.ID, Norm_Av))
  LoadPC2$Norm_Av[is.na(LoadPC2$Norm_Av)] <- 0
  write.csv(LoadPC2, paste0('LoadPC2-', ext,'.csv'), row.names=F)
  LoadPC3 <- suppressWarnings(loadings %>%
    rename(Norm_Av = PC3) %>%
    right_join(., Abbrev, by=c('Name'='Abb')) %>%
    select(KEGG.ID, Norm_Av))
  LoadPC3$Norm_Av[is.na(LoadPC3$Norm_Av)] <- 0
  write.csv(LoadPC3, paste0('LoadPC3-', ext,'.csv'), row.names=F)
}
