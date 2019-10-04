# test_plotReads.R

context("plotReads")

# ==== BEGIN SETUP AND PREPARE =================================================
#
alnsFile <- system.file("extdata", "mlx_reads.sorted.bam", package = "LSplicing")

load(system.file("extdata", "testCountReads.RData", package = "LSplicing"))
load(system.file("extdata", "testDataAlgn.Rdata", package = "LSplicing"))
coor <- testCountReads[[1]]
reads <- testCountReads[[2]]
#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(plotReads(coor, reads, "ENSG00000108788.11"),
               "argument \"alnsFile\" is missing, with no default")
})

test_that("a sample input prouces the expected output",  {
  expect_is(plotReads(coor, reads, "ENSG00000108788.11", alignments),
            "ggplot")
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#
rm(testCountReads)
rm(alignments)
# remove the output and my_index by file.remove()
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
