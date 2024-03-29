#' plotReads.R
#'
#' A function that plots results exon combination of target gene (plots results of countReads function).
#'
#' @param coor A data frame that contains exon coordinates of different genes.
#' @param reads A data frame containing rows of reads. If a value = 1,
#'     then that row of read is mapped to that column of exon. If a value = 0, then that row
#'     of read is not mapped to that column of exon.
#' @param gene A string of characters representing Ensembl gene ID of the gene of interest.
#' @param alnsFile The file that store the long-read RNA alignments generated by Minimap2.
#'
#' @return A ggplot graph that represent the long-read RNAs exon combinations
#' @examples
#' alignments <- system.file("extdata", "mlx_reads.sorted.bam",
#' package = "LSplicing")
#' referenceCoord <- system.file("extdata", "example_refCoord.gff3",
#' package = "LSplicing")
#'
#' \dontrun{
#' results <- countReads(alignments, referenceCoord)
#' coor <- results[[1]]
#' r <- results[[2]]
#' plotReads(coor = coor,
#'           reads = r,
#'           gene = "ENSG00000108788.11",
#'           alns = alignments)
#' }
#' @export
#' @import ggplot2
#' @import readr
#'

plotReads <- function(coor, reads, gene, alnsFile) {
  # narrow down the exon combinations to target gene
  reads <- reads[reads$geneId == gene, ]

  # remove columns with only NA, because it means this gene does not have those exons
  reads <- reads[,colSums(is.na(reads)) == 0]

  # remove reads that were not included in the countRead result
  idx <- reads[, "index"]
  alns <- GenomicAlignments::readGAlignments(alnsFile)
  alns <- alns[idx]

  # get coordinates for correpsonding gene
  speExon <- coor[coor$gene_id == gene, ]

  numAlignments <- length(alns)

  # calculates coverage of alignments
  cvg <- readCoverage(speExon, alns)

  # Initialize a ggplot that plot the coverage of alignments
  p <- ggplot2::ggplot(cvg, aes(x=pos, y=count)) +
    expand_limits(y=c(0, 3*numAlignments)) +
    coord_fixed(ratio = 5) +
    geom_area(color="#00AFBB", fill = "#00AFBB")


  # count and plot arc (missing exons in between)
  arc <- data.frame()

  for (i in 1:(ncol(reads) - 2)) {
    for (j in 1:(ncol(reads) - 2)) {
      if (abs(i - j) > 1 & (i < j)) {
        # get columns i to j
        s <- i + 2
        e <- j + 2
        tmp <- reads[,s:e]
        tmp <- tmp[rowSums(tmp) == 2,]
        first <- noquote(paste0("exon", i))
        sec <- noquote(paste0("exon", j))
        tmp <- tmp[tmp[[first]] == 1 & tmp[[sec]] == 1, ]
        count <- nrow(tmp)
        if (count > 0) {
          arc <- rbind(arc, data.frame("start" = speExon$end[i],
                                       "end" = speExon$start[j],
                                       "count" = count))
        }
      }
    }
  }

  # plot arcs

  for (i in 1:nrow(arc)) {
    start <- arc$start[i]
    end <- arc$end[i]
    ystart <- cvg$count[cvg$pos == start]
    yend <- cvg$count[cvg$pos == end]
    scale <- (ystart + yend)/2 + ((end - start)/12)
    p <- p +
      ggplot2::geom_curve(x = start,
                          xend = end,
                          y = ystart,
                          yend = yend,
                          curvature = -0.75,
                          size = 0.005,
                          color="#00AFBB") +
      ggplot2::geom_text(x = (end - start)/2 + start,
                         y = scale,
                         color = "grey45",
                         label = arc$count[i],
                         size = 2.5)
  }


  # add horizontal lines
  hr <- data.frame()

  # get start and end coordinates of horizontal lines
  for (i in 1:(nrow(speExon) - 1)) {
    first <- noquote(paste0("exon", i))
    sec <- noquote(paste0("exon", i+1))
    tmp <- reads[reads[[first]] == 1 & reads[[sec]] == 1, ]
    count <- nrow(tmp)
    hr <- rbind(hr, data.frame("start" = speExon$end[i],
                               "end" = speExon$start[i+1],
                               "count" = count))
  }

  # plot horizontal lines ans add label
  for (i in 1:nrow(hr)) {
    start <- hr$start[i]
    end <- hr$end[i]
    p <- p +
      ggplot2::geom_segment(x = start,
                            xend = end,
                            y = 0,
                            yend = 0,
                            size = 0.17,
                            color="#00AFBB") +
      ggplot2::geom_text(data = hr,
                         x = start + ((end - start)/2),
                         y = 10,
                         color = "grey45",
                         label = hr$count[i],
                         size = 2.5)
  }

  return(p)
}



# [END]
