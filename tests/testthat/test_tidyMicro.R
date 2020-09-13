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






