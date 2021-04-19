context("tidy_micro works")

testthat::test_that("OTU cols should be Lib names", {

  data("mrsa_clin"); data("mrsa_gen")

  expect_error( tidy_micro(mrsa_gen, tab_names = "G", mrsa_clin, library_name = "X",
                           count_summary = F) )


  mrsa_clin$Lib <- paste0(mrsa_clin$Lib, "_X")
  expect_error( tidy_micro(mrsa_gen, tab_names = "G", mrsa_clin, count_summary = F) )

})

testthat::test_that("Named lists and lists give same result", {

  data("mrsa_clin"); data("mrsa_gen")

  nmd <- tidy_micro(list("G" = mrsa_gen), clinical = mrsa_clin, count_summary = F)
  tbn <- tidy_micro(mrsa_gen, tab_names = "G", mrsa_clin, count_summary = F)

  expect_equivalent(nmd, tbn)
})

# testthat::test_that("otu_filter shoud give same reults", {
#
#   data("mrsa_clin"); data("mrsa_gen")
#
#   tidy.f <- tidy_micro(mrsa_gen, tab_names = "G", clinical = mrsa_clin,
#                        count_summary = F, filter_summary = F,
#                        prev_cutoff = 5, ra_cutoff = 1)
#
#   set <- tidy_micro(mrsa_gen, tab_names = "G", mrsa_clin,
#                      count_summary = F, filter_summary = F)
#
#   otu.f <- otu_filter(set, prev_cutoff = 5, ra_cutoff = 1)
#
#   expect_equivalent(tidy.f, otu.f)
# })





