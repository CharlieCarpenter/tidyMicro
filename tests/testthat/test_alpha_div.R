context("alpha_div")


test_that("min_goods and min_depth stoppers work",{
  skip_on_cran()
  data("mrsa_clin"); data("mrsa_gen")
  tm <- tidy_micro(list("G" = mrsa_gen), clinical = mrsa_clin, count_summary = F)

  expect_error(alpha_div(tm, table = "G", iter = 10, min_depth = 200000))
  expect_error(alpha_div(tm, table = "G", iter = 10, min_goods = 100))
})



