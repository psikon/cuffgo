\name{readRefmap}
  \alias{readRefmap}
  \title{reads output in .refmap format from cuffcompare}
  \description{This function reads output from the cuffcompare specific refmap format. It is an tab delimited file that lists, for each reference transcript, which cufflinks transcript either fully or partially match.}
  \usage{readRefmap(filename)}
  \arguments{
  \item{filename}{location of the cuffcompare output file in .refmap format}
  }
  \value{
  \item{dataframe}{content of the parsed file as dataframe}
  }
  \format{
  Tab delimited file, where one row represents a reference transcript with the following attributes:
  \describe{
    \item{\code{ref_gene_id}}{The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.}
    \item{\code{ref_id}}{The transcript_id attribute of the reference GTF record for this transcript}
    \item{\code{class_code}}{The type of match between the Cufflinks transcripts and the reference transcript.}
    \item{\code{cuff_id_list}}{A comma separated list of Cufflinks transcript ids matching the reference transcript}
  }
  }
  \references{Cufflinks Manual <cuff_in>.refmap: \url{http://cufflinks.cbcb.umd.edu/manual.html#refmap}}
 