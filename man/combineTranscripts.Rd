% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineTranscripts.R
\name{combineTranscripts}
\alias{combineTranscripts}
\title{A helper for countReads function}
\usage{
combineTranscripts(gList, dup, geneId = "")
}
\arguments{
\item{gList}{A GRanges object that store new exon coordinates.}

\item{dup}{A GRanges object of narrowed down reference coordinates of human genes. Only contain a target gene's protein coding exon.}

\item{geneId}{A string of characters that stores the target Ensembl gene ID.}
}
\value{
A list of two GRanges ojects.
}
\description{
A function that helps combine all transcript exons coordinates
}
\examples{
gRangesList <- GenomicRanges::GRanges()
referencesFile<- system.file("extdata", "example_refCoord.gff3",
package = "LSplicing")
refCoord <- rtracklayer::import(referencesFile)
gene <- "ENSG00000108788.11"

\dontrun{
combineTranscripts(gList = gRangesList,
                   dup = refCoord,
                   geneId = gene)
}
}
