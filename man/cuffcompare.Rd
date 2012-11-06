\name{cuffcompare}
          \docType{data}
          \alias{cuffcompare}
          \title{part of a cufflinks analysis}
          \description{
            This data set give the first 1000 lines of a RNA-Seq experiment in human B lymphoma cells. The cufflinks suite was used to assemble the transcripts and search for differenly expressed genes. For the identification of the transcripts and exons the Ensembl annotations were used.
          }
          \format{
  A data frame with 1000 observations on the following 16 variables from a cuffcompare run:
  \describe{
        \item{\code{chr}}{Chromosome or contig name}
    \item{\code{program}}{The name of the program that generated this file}
    \item{\code{feature}}{The type of record (always either "transcript" or "exon"}
    \item{\code{start}}{The leftmost coordinate of this record}
    \item{\code{end}}{The rightmost coordinate of this record}
    \item{\code{score}}{The most abundant isoform for each gene is assigned a score of 1000. Minor isoforms are scored by the ratio (minor FPKM/major FPKM)}
    \item{\code{strand}}{Cufflinks' guess for which strand the isoform came from. Always one of "+", "-", "."}
    \item{\code{frame}}{deprecated - field not used}
    \item{\code{gene_id}}{Cufflinks gene id}
    \item{\code{transcript_id}}{Cufflinks transcript id}
    \item{\code{exon_number}}{number of exon for this transcript}
    \item{\code{gene_name}}{The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.}
    \item{\code{oId}}{old identifier of this gene e.g. in cufflinks output}
    \item{\code{nearest_ref}}{The reference transcript to which the class code refers, if any}
    \item{\code{class_code}}{The type of relationship between the Cufflinks transcripts and the reference transcript}
    \item{\code{tss_id}}{The tss_id associated with the object, or "-" if not a transcript/primary transcript, or if tss_id isn't present}
}
}
\usage{data(cuffcompare)}
\references{Cuffcompare Output:\url{http://cufflinks.cbcb.umd.edu/manual.html#cuffcomp_output}}
          \keyword{datasets}