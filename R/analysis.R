mapGOtoData <- function(dataframe, genome, identifier, rowname,category) {
  library(goseq)
  GO <- getgo(dataframe$gene_id, genome, identifier, fetch.cats=c(paste0("GO:",category)))
  dataframe <- addgo(dataframe,GO,rowname)
  return(dataframe) 
}

mapGOtoElement <- function(element, genome, identifier,category) {
  library(goseq)
  return(getgo(element,genome,identifier,fetch.cats=c(paste0("GO:",category))))
}

mapAllGOtoData <- function(dataframe, genome, identifier, rowname.BP, rowname.CC, rowname.MF) {
  BP <- getgo(dataframe$gene_id, genome, identifier, fetch.cats=c("GO:BP"))
  CC <- getgo(dataframe$gene_id, genome, identifier, fetch.cats=c("GO:CC"))
  MF <- getgo(dataframe$gene_id, genome, identifier, fetch.cats=c("GO:MF"))
  dataframe <- addgo(dataframe,BP,rowname.BP)
  dataframe <- addgo(dataframe,CC,rowname.CC)
  dataframe <- addgo(dataframe,MF,rowname.MF)
}
