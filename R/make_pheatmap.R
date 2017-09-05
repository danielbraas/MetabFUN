#' This function takes in a data matrix and produces an annotated heatmap using the pheatmap package.
#' @author Daniel Braas
#' @param matrix The data matrix to be used
#' @param samples A data frame with information about the Conditions and Cell.Number
#' @param Norv The internal standard
#' @param Title The title to be used. This title will also be part of the file name.
#' @param heat.color The color scheme to be used for the heatmap.
#' @return An annotated heatmap saved as a pdf file.
#' @export

make_pheatmap <- function(matrix, samples=samples, heat.color=normal){
  if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Exp') {
    ext = 'Relative Amounts'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='MID') {
    ext = 'MIDs'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='FC') {
    ext = 'Fractional Contribution'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Labeled') {
    ext = 'Percent Labeled'
  }

  if (exists('samples')==F) samples <- info
  ann <- select(samples, Condition, Cell.Number) %>%
    as.data.frame()
  rownames(ann) <- colnames(matrix)
  if (exists('Norv')==T) {
    ann$Norvaline <- Norv
  } else ann$Norvaline <- 1

  ann_colors = list(
    Condition=colors[1:length(unique(gsub('_(.)*','',colnames(matrix))))],
    Norvaline = c("white", "blue"),
    Cell.Number = c("white", "green")
    )
  names(ann_colors[['Condition']]) <- unique(gsub('_(.)*','',colnames(matrix)))

  matrix[is.na(matrix)] <- 0
  heatmap.title=paste(Title, '-Heatmap-',ext,'.pdf', sep='')
  pheatmap::pheatmap(matrix, cluster_row=T, cluster_col=T,
                     clustering_distance_rows='correlation',
                     clustering_distance_cols='correlation',
                     color = colorRampPalette(heat.color)(100),
                     border_color="black", scale="row",
                     cellwidth = 20, cellheight = 10,
                     annotation=ann, annotation_colors = ann_colors,
                     show_colnames = F, main=paste(Title,ext,sep='-'),
                     filename=heatmap.title)
  dev.off()
}
