#' This function is searching either the KEGG or HMDB databases for potential annotations of ions.
#' @author Daniel Braas
#' @param MZ MZ is the mass to charge (m/z) of the ion(s) of interest
#' @param ppm ppm defines the allowed m/z deviation from MZ
#' @param adduct Defines the scope of potential ion clusters that will be searched for. This can be set to 'neg' or 'pos' to include all potential ion clusters in a given polarity or can be more specified (e.g. 'M+H')
#' @param Charge This variable defines the charge state of the MZ
#' @param db This can be used to restrict the search to either KEGG or HMDB. The default is to look in both.
#' @return A data frame with information about the putatively identified feature
#' @export


data_annotate <- function(MZ, ppm=10, adduct, Charge=1, db=''){
  dat <- filter(CompMZ, mz > (MZ - ppm * 1e-6 * MZ), mz < (MZ + ppm * 1e-6 * MZ))
  if (length(adduct) > 1) {dat <- filter(dat, Adduct %in% adduct)
  } else if (grepl('neg', adduct)==T) {dat <- filter(dat, Mode == 'negative')
  } else if (grepl('pos', adduct)==T) {dat <- filter(dat, Mode == 'positive')
  } else if (grepl('neu', adduct)==T) {dat <- filter(dat, Mode == 'neutral')
  } else {dat <- filter(dat, Adduct %in% adduct)
  }
  if (db == 'KEGG') {dat <- filter(dat, grepl('^C', dat$ID))
  } else if (db == 'HMDB') {dat <- filter(dat, grepl('HMDB', dat$ID))
  }
  dat <- filter(dat, charge %in% Charge)
  dat$Dppm = round((MZ - dat$mz) / dat$mz *1e6 ,3)
  if (nrow(dat) == 0) {
    for (i in 1:length(dat)){dat[1,i] <- NA}
  }
  dat$medMz = MZ
  dat <- select(dat, medMz, mz, ID, Name, Formula, MonoisotopicMass, Adduct, adductMass.x,
                num_molecules, charge, Mode, Type, Dppm)
  return(dat)
}
