\name{cufflinks}
          \docType{data}
          \alias{cufflinks}
          \title{part of a cufflinks analysis}
          \description{
            This data set give the first 1000 lines of a RNA-Seq experiment in human B lymphoma cells. The cufflinks suite was used to assemble the transcripts and search for differenly expressed genes. For the identification of the transcripts and exons the Ensembl annotations were used.
          }
          \format{A data frame with 1000 observations on the following 16 variables from a cufflinks run:
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
    \item{\code{fpkm}}{Isoform-level relative abundance in Fragments Per Kilobase of exon model per Million mapped fragment}
    \item{\code{frac}}{deprecated value}
    \item{\code{conf_lo}}{Lower bound of the 95\% confidence interval of the abundance of this isoform, as a fraction of the isoform abundance. That is, lower bound = FPKM * (1.0 - conf_lo)}
    \item{\code{conf_hi}}{Upper bound of the 95\% confidence interval of the abundance of this isoform, as a fraction of the isoform abundance. That is, upper bound = FPKM * (1.0 + conf_lo)}
    \item{\code{cov}}{Estimate for the absolute depth of read coverage across the whole transcript}
    \item{\code{full_read_support}}{When RABT assembly is used, this attribute reports whether or not all introns and internal exons were fully covered by reads from the data.}
}
  }
  \usage{data(cufflinks)}
  \references{Cufflinks Output: \url{http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output}}
          \keyword{datasets}