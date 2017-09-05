#' This function converts .raw files into .mzXML files for analysis with Maven and XCMS.
#' @details
#' For this function to work, one needs to define where the msconvert.exe file is located eg
#' D:\\pwiz\\msconvert.exe. This is the location of the program
#' you will use to convert the MS files. Also, you will need to go into the folder that contains
#' the .raw files. The function will create directories based on the polarity paramter that is specified.
#' @author Daniel Braas
#' @param polarity the output polariy. This can be either 'pos', 'neg', 'both', or 'no' for no polarity separation.
#' @export
make_mzXML <- function(polarity){
  wd=getwd()
  if (polarity=='neg'){
    neg <- paste0(wd, '/neg')
    dir.create(neg)
    FILES <- list.files(pattern='.raw')
    for (i in 1:length(FILES)){
      shell(paste(msconvert, "--mzXML --filter \"peakPicking true 1-\" --filter \"polarity negative\" ", FILES[i]))}
    file.rename(paste(wd, list.files(pattern='mz'), sep='/'), paste(wd,'neg', list.files(pattern='mz'), sep='/'))
  }

  else if (polarity=='pos'){
    pos <- paste0(wd, '/pos')
    dir.create(pos)
    FILES <- list.files(pattern='.raw')
    for (i in 1:length(FILES)){
      shell(paste(msconvert, "--mzXML --filter \"peakPicking true 1-\" --filter \"polarity positive\" ", FILES[i]))}
    file.rename(paste(wd, list.files(pattern='mz'), sep='/'), paste(wd,'pos', list.files(pattern='mz'), sep='/'))
  }

  else if (polarity=='both'){
    neg <- paste0(wd, '/neg')
    dir.create(neg)
    pos <- paste0(wd, '/pos')
    dir.create(pos)
    FILES <- list.files(pattern='.raw')
    for (i in 1:length(FILES)){
      shell(paste(msconvert, "--mzXML --filter \"peakPicking true 1-\" --filter \"polarity negative\" ", FILES[i]))}
    file.rename(paste(wd, list.files(pattern='mz'), sep='/'), paste(wd,'neg', list.files(pattern='mz'), sep='/'))
    for (i in 1:length(FILES)){
      shell(paste(msconvert, "--mzXML --filter \"peakPicking true 1-\" --filter \"polarity positive\" ", FILES[i]))}
    file.rename(paste(wd, list.files(pattern='mz'), sep='/'), paste(wd,'pos', list.files(pattern='mz'), sep='/'))
  }

  else if (polarity=='no'){
    FILES <- list.files(pattern='.raw')
    for (i in 1:length(FILES)){
      shell(paste(msconvert, "--mzXML --filter \"peakPicking true 1-\" ", FILES[i]))}
  }
}
