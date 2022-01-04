test_that("Check for correct char vector", {
 expect_identical(listPredefinedDataSets(), c("protCodingGenes","CGI","SNP","TFBS"))
})
