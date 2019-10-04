# test_countReads.R

context("countReads")

# ==== BEGIN SETUP AND PREPARE =================================================
#
alnsFile <- system.file("extdata", "mlx_reads.sorted.bam", package = "LSplicing")
referencesFile <- system.file("extdata", "example_refCoord.gff3", package = "LSplicing")

load(system.file("extdata", "testCountReads.RData", package = "LSplicing"))
load(system.file("extdata", "testDataAlgn.Rdata", package = "LSplicing"))
load(system.file("extdata", "testDataRef.Rdata", package = "LSplicing"))
#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(countReads(alnsFile),
               "argument \"referencesFile\" is missing, with no default")
  expect_error(countReads("alns.bam", referencesFile),
               "No such file, please check the path to files")
})

test_that("a sample input prouces the expected output",  {
  expect_equal(countReads(alignments, referenceCoord), testCountReads)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#
rm(testCountReads)
rm(alignments)
rm(referenceCoord)
file.remove(list.files("./", pattern = "^exon_coordinates.tsv"))
file.remove(list.files("./",pattern = "^exon_table.tsv"))
# remove the output and my_index by file.remove()
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
