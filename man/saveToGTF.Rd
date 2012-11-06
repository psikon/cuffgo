\name{saveToGTF}
         \alias{saveToGTF}
         \title{saves dataframes to cufflinks specific GTF format}
         \description{This wrapper function creates from a dataframe processed in R an GTF File in cufflinks specific GTF format. There will be different fields in GTF Files for every program of the cufflinks suite, so you have to choose the right program to save your data.}
         \usage{
         saveToGTF(filename, program = "cufflinks", dataframe)
         }
         \arguments{
           \item{filename}{path to the save location}
           \item{program}{for which program the GTF File will be created                     (default="cufflinks")}
           \item{dataframe}{saves this dataframe}
         }
         \examples{
         ## save cufflinks GTF File
         saveToGTF("test.gtf", "cufflinks", cufflinks)
         
         ## save cuffmerge GTF File
         saveToGTF("test.gtf", "cuffmerge", cuffmerge)
      
         ##save cuffcompare GTF File
         saveToGTF("test.gtf", "cuffcompare",cuffcompare)
         }
         \format{
         All three programs have 8 standard tab delimited values that have to be provided by every programm and a variable attribute field that is separated by semicolons. this attribute field is different between all three programs.
         Standard fields:
         \describe{
         \item{\code{chr}}{Chromosome or contig name}
         \item{\code{program}}{The name of the program that generated this file}
         \item{\code{feature}}{The type of record (always either "transcript" or "exon"}
         \item{\code{start}}{The leftmost coordinate of this record}
         \item{\code{end}}{The rightmost coordinate of this record}
         \item{\code{score}}{The most abundant isoform for each gene is assigned a score of 1000. Minor isoforms are scored by the ratio (minor FPKM/major FPKM)}
         \item{\code{strand}}{Cufflinks' guess for which strand the isoform came from. Always one of "+", "-", "."}
         \item{\code{frame}}{deprecated - field not used}
         }
         cufflinks specific attributes:
         \describe{
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
         cuffcompare specific attributes:
         \describe{
         \item{\code{gene_id}}{Cufflinks gene id}
    \item{\code{transcript_id}}{Cufflinks transcript id}
    \item{\code{exon_number}}{number of exon for this transcript}
    \item{\code{gene_name}}{The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.}
    \item{\code{oId}}{old identifier of this gene e.g. in cufflinks output}
    \item{\code{nearest_ref}}{The reference transcript to which the class code refers, if any}
    \item{\code{class_code}}{The type of relationship between the Cufflinks transcripts and the reference transcript}
    \item{\code{tss_id}}{The tss_id associated with the object, or "-" if not a transcript/primary transcript, or if tss_id isn't present}
         }
         cuffmerge specific attributes:
         \describe{
         \item{\code{gene_id}}{Cufflinks gene id}
    \item{\code{transcript_id}}{Cufflinks transcript id}
    \item{\code{exon_number}}{number of exon for this transcript}
    \item{\code{gene_name}}{The gene_name attribute of the reference GTF record for this transcript, if present. Otherwise gene_id is used.}
    \item{\code{oId}}{old identifier of this gene e.g. in cufflinks output}
    \item{\code{contained_in}}{is this transcript contained in a nother bigger transcript}
    \item{\code{nearest_ref}}{The reference transcript to which the class code refers, if any}
    \item{\code{class_code}}{The type of relationship between the Cufflinks transcripts and the reference transcript}
    \item{\code{tss_id}}{The tss_id associated with the object, or "-" if not a transcript/primary transcript, or if tss_id isn't present}
         }
         }
         \references{CUfflinks Manual: \url{http://cufflinks.cbcb.umd.edu/manual.html}}