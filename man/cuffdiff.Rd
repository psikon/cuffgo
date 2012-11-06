\name{cuffdiff}
          \docType{data}
          \alias{cuffdiff}
          \title{part of a cuffdiff analysis}
          \description{
            This data set give the first 1000 lines of a RNA-Seq experiment in human B lymphoma cells. The cufflinks suite was used to assemble the transcripts and search for differenly expressed genes. For the identification of the transcripts and exons the Ensembl annotations were used.
          }
          \format{A data frame with 999 observations on the following 16 variables from a cuffdiff (V2.0) run: 
          \describe{
    \item{\code{test_id}}{A unique identifier describing the gene being tested. }
    \item{\code{gene_id}}{The gene_name or gene_id}
    \item{\code{gene}}{The gene_name(s) or gene_id(s) being tested}
    \item{\code{locus}}{Genomic coordinates for easy browsing to the genes or transcripts being tested.}
    \item{\code{sample_1}}{Label (or number if no labels provided) of the first sample being tested}
    \item{\code{sample_2}}{Label (or number if no labels provided) of the second sample being tested}
    \item{\code{status}}{Can be one of OK (test successful), NOTEST (not enough alignments for testing), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents testing.}
    \item{\code{value_1}}{FPKM of the gene in sample x}
    \item{\code{value_2}}{FPKM of the gene in sample y}
    \item{\code{log2.fold_change.}}{The (base 2) log of the fold change y/x}
    \item{\code{test_stat}}{The value of the test statistic used to compute significance of the observed change in FPKM}
    \item{\code{p_value}}{The uncorrected p-value of the test statistic}
    \item{\code{q_value}}{The FDR-adjusted p-value of the test statistic}
    \item{\code{significant}}{Can be either "yes" or "no", depending on whether p is greater then the FDR after Benjamini-Hochberg correction for multiple-testing}
    }}
          \usage{data(cuffdiff)}
          \references{Cuffdiff Output:\url{http://cufflinks.cbcb.umd.edu/manual.html#gene_exp_diff}}
          \keyword{datasets}