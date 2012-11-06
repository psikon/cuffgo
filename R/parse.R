readCufflinksFiles <- function(filename, program = "cufflinks") {

  program <- match.arg(program, c("cufflinks", "cuffcompare", 
                                  "cuffmerge","cuffdiff"))
  
  switch(program,
         cufflinks=readCufflinks(filename),
         cuffcompare=readCuffCompare(filename),
         cuffmerge=readCuffMerge(filename),
         cuffdiff=readCuffDiff(filename)
         )
}

readCufflinks <- function(filename) {
      c.return <- .Call(parseCufflinks,as.character(filename),PACKAGE="cuffGO")
      dataframe <- data.frame(c.return$chr,c.return$program,
                                       c.return$feature,c.return$start,
                                       c.return$end,c.return$score,
                                       c.return$strand,c.return$frame,
                                       c.return$gene_id,c.return$transcript_id,
                                       c.return$exon_number,c.return$fpkm,
                                       c.return$frac,c.return$conf_lo,
                                       c.return$conf_hi,c.return$cov,
                                       c.return$frs,stringsAsFactors=FALSE)
      colnames(dataframe) <- c("chr","program","feature","start","end",
                                "score","strand","frame","gene_id",
                                "transcript_id","exon_number","fpkm","frac","conf_lo",
                                "conf_hi","cov","full_read_support") 
      return(dataframe)
}

readCufflinksTranscripts <- function(filename) {
  dataframe <- readCufflinks(filename)
  dataframe <- dataframe[grep("transcript",dataframe$feature),]
  return (dataframe)
}

readCufflinksExons <- function(filename) {
  dataframe <- readCufflinks(filename)
  dataframe <- dataframe[grep("exon",dataframe$feature),]
  return (dataframe)
}

readCuffMerge <- function(filename) {
  c.return <- .Call(parseCuffmerge,as.character(filename),PACKAGE="cuffGO")
  dataframe <- data.frame(c.return$chr,c.return$program,
                          c.return$feature,c.return$start,
                          c.return$end,c.return$score,
                          c.return$strand,c.return$frame,
                          c.return$gene_id,c.return$transcript_id,
                          c.return$exon_number,c.return$gene_name,
                          c.return$oId,c.return$contained_in,
                          c.return$nearest_ref,c.return$class_code,
                          c.return$tss_id,stringsAsFactors=FALSE)
  colnames(dataframe) <- c("chr","program","feature","start","end",
                           "score","strand","frame","gene_id",
                           "transcript_id","exon_number",
                           "gene_name","oId","contained_in",
                           "nearest_ref","class_code","tss_id") 
  return(dataframe)
}


readCuffCompare <- function(filename) {
  c.return <- .Call(parseCuffcompare,as.character(filename),PACKAGE="cuffGO")
  dataframe <- data.frame(c.return$chr,c.return$program,
                          c.return$feature,c.return$start,
                          c.return$end,c.return$score,
                          c.return$strand,c.return$frame,
                          c.return$gene_id,c.return$transcript_id,
                          c.return$exon_number,c.return$gene_name,
                          c.return$oId,c.return$nearest_ref,
                          c.return$class_code,c.return$tss_id,
                          stringsAsFactors=FALSE)
  colnames(dataframe) <- c("chr","program","feature","start","end",
                           "score","strand","frame","gene_id",
                           "transcript_id","exon_number",
                           "gene_name","oId","nearest_ref",
                           "class_code","tss_id") 
  return(dataframe)  
}


readCuffDiff <- function(filename)
{
  return(read.table(filename, header=TRUE, sep="\t",as.is=TRUE))
}

readTracking <- function(filename) {
  return(read.table(filename, header=TRUE, sep="\t",as.is=TRUE))
  
}

readRefmap <- function(filename) {
  dataframe <- read.table(filename, header=TRUE, sep="\t",as.is=TRUE)
  colnames(dataframe) <- c("ref_gene_id","ref_id","class_code","cuff_id_list")
  return(dataframe)
}

readTmap <- function(filename) {
  dataframe <- read.table(filename, header=TRUE, sep="\t",as.is=TRUE)
  colnames(dataframe) <- c("Reference_gene_name","Reference_transcript_id",
                           "Class_code","cufflinks_gene_id",
                           "cufflinks_transcript_id","FMI", "FPKM",
                           "FPKM_conf_lo","FPKM_conf_hi","coverage",
                           "length","major_isoform_ID", "ref_match_len")
  return(dataframe)
}

