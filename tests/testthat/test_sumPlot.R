test_that("Check class of sumPlot", {
 myOGRE <- makeExampleOGREDataSet()
 myOGRE <- loadAnnotations(myOGRE)
 myOGRE <- fOverlaps(myOGRE)
 myOGRE <- sumPlot(myOGRE)
 expect_equal(class(myOGRE),class(GRangesList()))
})


