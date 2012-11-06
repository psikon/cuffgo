\name{mapAllGOtoData}
         \alias{mapAllGOtoData}
         \title{seach all GO ontologies for complete dataframes}
         \description{function to search in all three GO ontologies for GO terms for the complete dataframe. This function creates lists of GO terms seperated by the three major ontologies and appends them to the dataframe.}
         \usage{
         mapAllGOtoData(dataframe,genome,identifier,rowname.BP,rowname.CC,rowname.MF)
         }
         \arguments{
           \item{dataframe}{dataframe where the GO Identifier will be added}
           \item{genome}{which genome will be used}
           \item{identifier}{which type of identifiers will be used e.g. ensembl="ensgene",...} 
           \item{rowname.BP}{name of the row that will be added}
           \item{rowname.CC}{name of the row that will be added}
           \item{rowname.MF}{name of the row that will be added}
          
         }
         \examples{
         ## map to all Elements in the dataframe all Terms of the three ontologies using ensembl annotation
         mapAllGOtoData(cuffdiff, "hg19", "ensGene","BP","CC","MF")
         }
         \seealso{\link{mapGOtoData}, \link{mapAllGOtoData}, \link[goseq]{goseq}
         }