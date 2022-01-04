test_that("Check class of addDataSetFromHub", {
 myOGRE <- makeExampleOGREDataSet()
 myOGRE <- loadAnnotations(myOGRE)
 myOGRE <- addDataSetFromHub(myOGRE,"protCodingGenes","query")
 expect_equal(class(myOGRE),class(GRangesList()))
})
