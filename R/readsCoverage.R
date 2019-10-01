# readCoverage.R

readCoverage <- function(coor, alns) {

  cvg <- coverage(alns)
  ans <- GRanges(cvg)

  # get certain chromosome

  a <- subset(eval(ans), seqnames == "chr17")

  # check IRanges
  first <- which(start(a) <= coor[1,"start"] & coor[1,"start"] <= end(a))
  last <- which(start(a) <= coor[nrow(coor), "end"] & coor[nrow(coor), "end"] <= end(a))

  a <- a[first:last]

  cvg <- data.frame()
  exon <- 1
  j <- 1


  while (j <= length(a) && exon <= nrow(coor)) {
    e <- coor[exon, "end"]
    s <- coor[exon, "start"]
    if (s > end(a[j])) {
      j <- j + 1
    } else if (s > start(a[j]) && e > end(a[j])) {
      tmp <- data.frame("pos"= s:end(a[j]), "count"= rep(a[j]$score, end(a[j])-s+1))
      cvg <- rbind(cvg, tmp)
      j <- j + 1
    } else if (s <= start(a[j]) && e > end(a[j])) {
      tmp <- data.frame("pos"= start(a[j]):end(a[j]), "count"= rep(a[j]$score, width(a[j])))
      cvg <- rbind(cvg, tmp)
      j <- j + 1
    } else if (s > start(a[j]) && e < end(a[j])) {
      tmp <- data.frame("pos"= s:e, "count"= rep(a[j]$score, e-s+1))
      cvg <- rbind(cvg, tmp)
      exon <- exon + 1
    } else if (s <= start(a[j]) && e < end(a[j])) {
      tmp <- data.frame("pos"= start(a[j]):e, "count"= rep(a[j]$score, e-start(a[j])+1))
      cvg <- rbind(cvg, tmp)
      exon <- exon + 1
    } else if (s <= start(a[j]) && e == end(a[j])) {
      tmp <- data.frame("pos"= start(a[j]):end(a[j]), "count"= rep(a[j]$score, width(a[j])))
      cvg <- rbind(cvg, tmp)
      j <- j + 1
      exon <- exon + 1
    } else if (s > start(a[j]) && e == end(a[j])) {
      tmp <- data.frame("pos"= s:end(a[j]), "count"= rep(a[j]$score, end(a[j])-s+1))
      cvg <- rbind(cvg, tmp)
      j <- j + 1
      exon <- exon + 1
    }
  }

  # Get intron position and assign 0 to them
  start <- coor$start[1]
  end <- coor$end[nrow(coor)]
  all <- start:end
  intron <- all[!all %in% cvg$pos]

  # add the intron positions to cvg
  tmp <- data.frame("pos"=intron, "count"= rep(0, length(intron)))
  cvg <- rbind(cvg, tmp)

  return(cvg)
}



# [END]





