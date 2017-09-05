path <- keggList('pathway','hsa')
pathway <- as.data.frame(path)
pathway$ID <- gsub('path:','', names(path))
pathway$path <- gsub(' - H(.)*', '', pathway$path)
names(pathway)[1] <- 'Name'
rownames(pathway) <- NULL

paths <- list()
for (i in 1:nrow(pathway)){
  paths[[i]] <- keggGet(pathway$ID[i])
}

pathway$Class <- as.character(sapply(paths, function(x) x[[1]]$CLASS))

#find out which pathway has the most compounds
which.max(sapply(paths, function(x) length(x[[1]]$COMPOUND)))

path <- data.frame('Name'=pathway$Name[1], 
                   'Compound'=paths[[1]][[1]]$COMPOUND, 
                   'KEGG'=names(paths[[1]][[1]]$COMPOUND))
for (i in 2:nrow(pathway)){
  if (length(paths[[i]][[1]]$COMPOUND) > 0) {
    assign(pathway$Name[i], data.frame('Name'=pathway$Name[i], 
                                       'Compound'=paths[[i]][[1]]$COMPOUND, 
                                       'KEGG'=names(paths[[i]][[1]]$COMPOUND)))
  } else {
    assign(pathway$Name[i], data.frame('Name'=pathway$Name[i], 
                                       'Compound'= NA, 
                                       'KEGG'= NA))
  }
  path <- bind_rows(path, get(pathway$Name[i]))
}
compound.pathway <- left_join(pathway, path, by='Name')

# make data frame with pathways and genes ---------------------------------

path <- data.frame('Name'=pathway$Name[1], 
                   'Gene'=paths[[1]][[1]]$GENE[seq(2,length(paths[[1]][[1]]$GENE), by=2)], 
                   'Entry'=paths[[1]][[1]]$GENE[seq(1,length(paths[[1]][[1]]$GENE), by=2)])
for (i in 2:nrow(pathway)){
  if (length(paths[[i]][[1]]$GENE) > 0) {
    assign(pathway$Name[i], data.frame('Name'=pathway$Name[i], 
                                       'Gene'=paths[[i]][[1]]$GENE[seq(2,length(paths[[i]][[1]]$GENE), by=2)], 
                                       'Entry'=paths[[i]][[1]]$GENE[seq(1,length(paths[[i]][[1]]$GENE), by=2)]))
  } else {
    assign(pathway$Name[i], data.frame('Name'=pathway$Name[i], 
                                       'Gene'= NA, 
                                       'Entry'= NA))
  }
  path <- bind_rows(path, get(pathway$Name[i]))
}
gene.pathway <- left_join(pathway, path, by='Name')
gene.pathway$Symbol <- gsub(';(.)*','',gene.pathway$Gene)
gene.pathway$EC <- gsub('(.)*EC:|](.)*','', gene.pathway$Gene)
gene.pathway$KO <- gsub('(.)*KO:|](.)*','', gene.pathway$Gene)
gene.pathway$Description <- gsub('(.)*; | \\[(.)*','', gene.pathway$Gene)
gene.pathway$Gene <- NULL

metabolism.gene.pathway <- filter(gene.pathway, grepl('Metabolism;', gene.pathway$Class))
metabolism.gene.pathway <- group_by(metabolism.gene.pathway, ID) %>%
  mutate(Gene = paste0('Gene',1:n())) %>% 
  ungroup()
metabolism.gene.pathway$Gene <- factor(metabolism.gene.pathway$Gene, levels = unique(metabolism.gene.pathway$Gene))
metab.genes <- sort(unique(metabolism.gene.pathway$Symbol))

KEGG.gmt <- select(metabolism.gene.pathway, ID, Gene, Symbol) %>% 
  spread(Gene, Symbol) %>% 
  left_join(., metabolism.gene.pathway, by = 'ID') %>% 
  select(Name, ID, contains('Gene'), -Gene) %>% 
  distinct()

write.table(KEGG.gmt, 'KEGG_metabolism.gmt',sep='\t', quote=F, row.names=F) 

# getting all KEGG compounds ----------------------------------------------
comp <- keggList('compound')
COMP <- data.frame('ID'=rep(NA,length(comp)), 
                   'Name'=rep(NA,length(comp)), 
                   'Formula'=rep(NA,length(comp)), 
                   'Exact.Mass'=rep(NA,length(comp)))

for(i in 1001:length(comp)){
  print(i)
  COMP$ID[i] <- gsub('cpd:','',names(comp[i]))
  a <- keggGet(COMP$ID[i])
  COMP$Name[i] <- paste(a[[1]]$NAME, collapse = "")
  if (is.null(a[[1]]$FORMULA)==T) {COMP$Formula[i] <- NA
  } else {COMP$Formula[i] <- a[[1]]$FORMULA}
  if (is.null(a[[1]]$EXACT_MASS)==T) {COMP$Exact.Mass[i] <- NA
  } else {COMP$Exact.Mass[i] <- a[[1]]$EXACT_MASS}
}
compounds <- COMP
save(compounds, file='KEGG_compounds.rda')
