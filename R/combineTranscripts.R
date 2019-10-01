# combineTranscripts.R

combineTranscripts <- function(gList, dup, geneId) {
  count <- 1
  start <- GenomicRanges::start(dup[1])
  end <- GenomicRanges::end(dup[1])
  newGranges <- GenomicRanges::GRanges()
  prev <- dup[1]

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
