# test_readCoverage.R

context("readCoverage")

# ==== BEGIN SETUP AND PREPARE =================================================
#
load(system.file("extdata", "testCountReads.RData", package = "LSplicing"))
load(system.file("extdata", "testDataAlgn.Rdata", package = "LSplicing"))
load(system.file("extdata", "testReadCovAlgn.Rdata", package = "LSplicing"))
load(system.file("extdata", "testReadCoverage.Rdata", package = "LSplicing"))
coor <- testCountReads[[1]]
coor <- coor[coor$gene_id == "ENSG00000108788.11", ]
#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(readCoverage(coor),
               "argument \"alns\" is missing, with no default")
})

test_that("a sample input prouces the expected output",  {
  expect_equal(readCoverage(coor, testReadAlns), testReadCoverage)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#
rm(testReadAlns)
rm(testReadCoverage)

# remove the output and my_index by file.remove()
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
