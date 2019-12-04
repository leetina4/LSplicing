#' A helper for countReads function
#'
#' A function that helps combine all transcript exons coordinates
#'
#' @param gList A GRanges object that store new exon coordinates.
#' @param dup A GRanges object of narrowed down reference coordinates of human genes. Only contain a target gene's protein coding exon.
#' @param geneId A string of characters that stores the target Ensembl gene ID.
#'
#' @return A list of two GRanges ojects.
#'
#' @examples
#' gRangesList <- GenomicRanges::GRanges()
#' referencesFile<- system.file("extdata", "example_refCoord.gff3",
#' package = "LSplicing")
#' refCoord <- rtracklayer::import(referencesFile)
#' gene <- "ENSG00000108788.11"
#'
#' \dontrun{
#' combineTranscripts(gList = gRangesList,
#'                    dup = refCoord,
#'                    geneId = gene)
#' }
#' @import GenomicRanges
#' @importFrom IRanges IRanges

combineTranscripts <- function(gList, dup, geneId = "") {
  count <- 1
  start <- GenomicRanges::start(dup[1])
  end <- GenomicRanges::end(dup[1])
  newGranges <- GenomicRanges::GRanges()
  prev <- dup[1]
  # create New range storing new exon coordinates
  for (i in 2:length(dup)) {
    if (i == length(dup)) {
      newGranges <- GenomicRanges::union(newGranges,
                                         GenomicRanges::GRanges(dup[i]@seqnames,
                                                                ranges = IRanges::IRanges(start,
                                                                                          end, end - start + 1)))

      tmp <- data.frame("seqnames" = dup[i]@seqnames,
                        "gene_id" = geneId,
                        "exon" = count,
                        "start" = start,
                        "end" = end,
                        "width" = end - start + 1)
      if (! (nrow(merge(gList, tmp)) > 0)) {
        gList <- rbind(gList, tmp)
      }

    }
    if (GenomicRanges::start(dup[i]) > end) {
      # create New range
      newGranges <- GenomicRanges::union(newGranges,
                                         GenomicRanges::GRanges(dup[i]@seqnames,
                                                                ranges = IRanges::IRanges(start,
                                                                                          end, end - start + 1)))

      tmp <- data.frame("seqnames" = dup[i]@seqnames,
                        "gene_id" = geneId,
                        "exon" = count,
                        "start" = start,
                        "end" = end,
                        "width" = end - start + 1)
      if (! (nrow(merge(gList, tmp)) > 0)) {
        gList <- rbind(gList, tmp)
      }


      count <- count + 1
      start <- GenomicRanges::start(dup[i])
      end <- GenomicRanges::end(dup[i])
    } else {
      end <- GenomicRanges::end(dup[i])

    }
    prev <- dup[i]
  }
  result <- list(newGranges, gList)
  return(result)
}

# [END]
