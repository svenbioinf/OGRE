test_that("Check class of OGREDataSetFromDir", {
 myQueryFolder <- file.path(system.file('extdata', package = 'OGRE'),"query")
 mySubjectFolder <- file.path(system.file('extdata', package = 'OGRE'),"subject")
 myOGRE <- OGREDataSetFromDir(queryFolder=myQueryFolder,subjectFolder=mySubjectFolder)
 expect_equal(class(myOGRE),class(GRangesList()))
})
