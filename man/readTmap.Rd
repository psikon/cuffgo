\name{readTmap}
  \alias{readTmap}
  \title{reads output in .tmap format from cuffcompare}
  \description{This function reads output from the cuffcompare specific tmap format. It is an tab delimited file that lists the most closely matching reference transcripts for each Cufflinks transcript.}
  \usage{readTmap(filename)}
  \arguments{
  \item{filename}{location of the cuffcompare output file in .tmap format}
  }
  \value{
  \item{dataframe}{content of the parsed file as dataframe}
  }
  \format{
  Tab delimited file, where one row represents a Cufflinks transcripts with the following attributes:
  \describe{
  \item{\code{Reference_gene_name}}{The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.}
  \item{\code{Referenc_transcript_id}}{The transcript_id attribute of the reference GTF record for this transcript}
  \item{\code{Class_code}}{The type of relationship between the Cufflinks transcripts and the reference transcript}
  \item{\code{Cufflinks_gene_id}}{The Cufflinks internal gene id}
  \item{\code{Cufflinks_transcript_id}}{The Cufflinks internal transcript id}
  \item{\code{FMI}}{Fraction of major isoform: the expression of this transcript expressed as a fraction of the major isoform for the gene. Ranges from 1 to 100.}
  \item{\code{FPKM}}{The expression of this transcript expressed in FPKM}
  \item{\code{FPKM_conf_lo}}{The lower limit of the 95\% FPKM confidence interval}
  \item{\code{FPKM_conf_hi}}{The upper limit of the 95\% FPKM confidence interval}
  \item{\code{coverage}}{The estimated average depth of read coverage across the transcript.}
  \item{\code{length}}{The length of the transcript.}
  \item{\code{major_isoform_ID}}{The Cufflinks ID of the gene's major isoform}
  \item{\code{ref_match_len}}{The length of the matching positions between the cufflinks and the reference transcript}
  }
  }
  \references{Cufflinks Manual <cuff_in>.tmap: \url{http://cufflinks.cbcb.umd.edu/manual.html#tmap}}