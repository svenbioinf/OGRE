test_that("Check class of loadAnnotations", {
 myOGRE <- makeExampleOGREDataSet()
 myOGRE <- loadAnnotations(myOGRE)
 expect_equal(class(myOGRE),class(GRangesList()))
})
