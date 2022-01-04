test_that("Check class of addGRanges", {
 myOGRE <- makeExampleOGREDataSet()
 myOGRE <- loadAnnotations(myOGRE)
 myOGRE <- addGRanges(myOGRE,makeExampleGRanges(),"subject")
 expect_equal(class(myOGRE),class(GRangesList()))
})
