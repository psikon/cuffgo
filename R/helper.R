addgo <- function (dataframe, liste, spaltenname) {
  tmp <- lapply(liste, paste, collapse=",")
  dataframe <- cbind(dataframe, unlist(tmp), stringsAsFactors=FALSE)
  colnames(dataframe)[ncol(dataframe)] <- spaltenname
  dataframe
}