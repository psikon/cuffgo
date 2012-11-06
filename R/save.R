#' wrapper function to save dataframes in cufflinks gtf format
#' 
#' This function generates gtf files from dataframes in cufflinks format.  
#'
#' @param filename      path of the output file
#' @param program       output of which program have to be parsed
#' @param dataframe     dataframe to be written
#'
#' @export
#' 
#'
saveToGTF <- function(filename, program ="cufflinks", dataframe) {
  program <- match.arg(program,c("cufflinks","cuffcompare","cuffmerge"))
  switch(program, cufflinks=saveCufflinksGTF(filename, dataframe),
         cuffcompare=saveCuffCompareGTF(filename, dataframe),
         cuffmerge=saveCuffMergeGTF(filename, dataframe))
}

saveCufflinksGTF <- function(filename, dataframe) {
  c.return <- .Call(writeCufflinksGTF,as.character(filename),
                    as.character(dataframe$chr),
                    as.character(dataframe$program),
                    as.character(dataframe$feature),
                    as.integer(dataframe$start),
                    as.integer(dataframe$end),
                    as.integer(dataframe$score),
                    as.character(dataframe$strand),
                    as.character(dataframe$frame),
                    as.character(dataframe$gene_id),
                    as.character(dataframe$transcript_id),
                    as.integer(dataframe$exon_number),
                    as.numeric(dataframe$fpkm),
                    as.numeric(dataframe$frac),
                    as.numeric(dataframe$conf_lo),
                    as.numeric(dataframe$conf_hi),
                    as.numeric(dataframe$cov),
                    as.character(dataframe$full_read_support),
                    dim(dataframe)[1],PACKAGE="cuffGO")
}

saveCuffCompareGTF <- function(filename, dataframe) {
  c.return <- .Call(writeCuffCompareGTF,as.character(filename),
                    as.character(dataframe$chr),
                    as.character(dataframe$program),
                    as.character(dataframe$feature),
                    as.integer(dataframe$start),
                    as.integer(dataframe$end),
                    as.character(dataframe$score),
                    as.character(dataframe$strand),
                    as.character(dataframe$frame),
                    as.character(dataframe$gene_id),
                    as.character(dataframe$transcript_id),
                    as.integer(dataframe$exon_number),
                    as.character(dataframe$gene_name),
                    as.character(dataframe$oId),
                    as.character(dataframe$nearest_ref),
                    as.character(dataframe$class_code),
                    as.character(dataframe$tss_id),
                    dim(dataframe)[1],PACKAGE="cuffGO")
}

saveCuffMergeGTF <- function(filename, dataframe) {
  c.return <- .Call(writeCuffMergeGTF,as.character(filename),
                    as.character(dataframe$chr),
                    as.character(dataframe$program),
                    as.character(dataframe$feature),
                    as.integer(dataframe$start),
                    as.integer(dataframe$end),
                    as.character(dataframe$score),
                    as.character(dataframe$strand),
                    as.character(dataframe$frame),
                    as.character(dataframe$gene_id),
                    as.character(dataframe$transcript_id),
                    as.integer(dataframe$exon_number),
                    as.character(dataframe$gene_name),
                    as.character(dataframe$oId),
                    as.character(dataframe$contained_in),
                    as.character(dataframe$nearest_ref),
                    as.character(dataframe$class_code),
                    as.character(dataframe$tss_id),
                    dim(dataframe)[1],PACKAGE="cuffGO")
}

saveToCSV <- function(filename,dataframe) {
  write.csv(dataframe,filename,quote=FALSE, 
            na="NA",row.names=TRUE,)
}