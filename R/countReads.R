# countReads.R

#'
#' \code{countReads} Return a list of domains that the genes in two systems
#'    overlap along with genes that mapped to the same domains
#'
#' @param alns (GenomicAlignment object) The name of the system in interest. Is compared with system2
#'     in this function.
#' @param refCoord (GRanges) The name of the system in interest. Is compared with system1
#'     in this function.
#' @return None
#' @author {Tina Lee} (aut)
#' @examples
#' # Get Alignment data and corresponding reference coordinates and run following code
#' countReads(alns, refCoord)
#'
#' @export

countReads <- function(alns, refCoord) {
  results <- data.frame()
  gList <- GenomicRanges::GRanges()

  # Extract gene in refCo  ord
  idx <- refCoord@elementMetadata@listData[["type"]] == "gene"
  genes <- refCoord[idx]

  for (k in seq_along(alns)) {
    gene <- data.frame(union = SummarizedExperiment::assay(GenomicAlignments::summarizeOverlaps(genes, alns[k], inter.feature=FALSE)))

    # If this does not equal to 1, then skip this read
    if (sum(gene$reads) == 1) {
      (matchIdx <- which(gene[,1] > 0))
      genes[matchIdx]

      geneId <- genes[matchIdx]@elementMetadata@listData[["gene_id"]]

      dup <- sort(refCoord[refCoord@elementMetadata@listData[["gene_id"]] == geneId &
                             refCoord@elementMetadata@listData[["type"]] == "exon" &
                             refCoord@elementMetadata@listData[["transcript_type"]] == "protein_coding"])

      combined <- combineTranscripts(gList, dup, geneId)
      newGranges <- combined[[1]]
      gList <- combined[[2]]

      read <- data.frame(union = SummarizedExperiment::assay(GenomicAlignments::summarizeOverlaps(newGranges,
                                                                                                  alns[k], inter.feature=FALSE), byrow = TRUE))
      sum(read$reads)

      # Transpose the data frame
      t <- t(read)
      colnames(t) <- paste(rep("exon", ncol(t)), sep = "", 1:ncol(t))
      newDf <- cbind(data.frame("index" = k, geneId), t)

      # Assign row name
      ## use plyr::rbind.fill(x,y) to bind data frames
      results <- plyr::rbind.fill(results, newDf)
    }
  }

  # finish for loop
  # output a txt file (could use write.table())
  write.table(gList, file = "exon_coordinates.tsv", quote = FALSE, sep = "\t", col.names = NA)
  write.table(results, file = "exon_table.tsv", quote = FALSE, sep = "\t", col.names = NA)
}



# [END]
