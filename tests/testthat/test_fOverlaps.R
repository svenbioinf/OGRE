test_that("Check class of fOverlaps", {
 myOGRE <- makeExampleOGREDataSet()
 myOGRE <- loadAnnotations(myOGRE)
 myOGRE <- fOverlaps(myOGRE)
 expect_equal(class(myOGRE),class(GRangesList()))
})

