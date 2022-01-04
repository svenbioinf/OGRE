test_that("Check class of gvizPlot", {
 myOGRE <- makeExampleOGREDataSet()
 myOGRE <- loadAnnotations(myOGRE)
 myOGRE <- fOverlaps(myOGRE)
 myOGRE <- gvizPlot(myOGRE,query="ENSG00000142168")
 expect_equal(class(myOGRE),class(GRangesList()))
})

