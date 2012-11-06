\name{mapGOtoElement}
         \alias{mapGOtoElement}
         \title{search GO terms for one element}
         \description{search for one element all specific GO terms of       one ontology}
         \usage{
         mapGOtoElement(element,genome,identifier,category)
         }
         \arguments{
           \item{element}{Identifier to search for GO Terms e.g. Ensembl Accession number}
           \item{genome}{genome for searching}
           \item{identifier}{type of identifier}
           \item{category}{which ontologie will be used e.g. BP - biological process, CC - cellular component, MF - molecular function}

         }
         \examples{
         ## map GO Term to one Element
         mapGOtoElement("ENSG00000000003", "hg19","ensGene","BP")
         }
         \seealso{\link{mapGOtoData}, \link{mapAllGOtoData}, \link[goseq]{goseq}
         }