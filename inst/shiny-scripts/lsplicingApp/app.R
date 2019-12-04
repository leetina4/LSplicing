# app.R
#
# Ref: https://shiny.rstudio.com/tutorial/ Shiny tutorial
#

# Load packages ----
library(shiny)
library(ggplot2)
library(GenomicAlignments)

# Load data ----
load("/Users/tinalee/Documents/BCB410/LSplicing/inst/extdata/testCountReads.RData")

# Define UI ----
myUi <- fluidPage(

  titlePanel("Splicing graph on long-rna reads and reference sequence"),
  sidebarLayout( position = "right",
    sidebarPanel(
      helpText(h3("Please input reads file and reference genome file as needed.")),
      fileInput(inputId = "readsFile",
                label = "Input reads file:"),
      fileInput(inputId = "refGenomeFile",
                label = "Input reference genome file:"),
      # Horizontal line
      tags$hr(),
      helpText(h3("Please enter target gene ID")),
      textInput(inputId = "geneID",
                label = "Gene ID:")
    ),

    mainPanel(
      textOutput(outputId = "message"),
      plotOutput("plot")
    )
  )
)

# Define Server ----
myServer <- function(input, output) {
  observe({
    if (is.null(input$readsFile) || is.null(input$refGenomeFile) || is.null(input$geneID)) {
      message <-
        sprintf("Default reads file is mlx_reads.sorted.bam in extdata of LSplicing package.
                Default reference Genome File is example_refCoord.gff3 in extdata of LSplicing package.
                Please upload reads file or transcriptsFile.")
      count <- testCountReads
      alnsFile <- system.file("extdata", "mlx_reads.sorted.bam", package = "LSplicing")
      gene <- "ENSG00000108788.11"
    } else {
      message <-
        sprintf("Uploaded reads file is %s.
                Uploaded transcripts file is %s.",
                input$readsFile$name, input$refGenomeFile$name)
      count <- countReads(input$readsFile$datapath, input$refGenomeFile$datapath)
      alnsFile <- input$readsFile$datapath
      gene <- input$geneID
    }

    output$message <- renderText({message})

    output$plot <- renderPlot({
      coor <- count[[1]]
      reads <- count[[2]]

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
        coord_fixed(ratio = 6) +
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
      print(p)
    })
  })
}

# Run APP ----
shinyApp(ui = myUi, server = myServer)

# [END]
