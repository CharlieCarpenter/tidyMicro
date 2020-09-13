context("cor_rocky_mtn")

test_that("cor_rocky_mtn takes numeric values only", {

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  expect_error( cor_rocky_mtn(set, table = "P", bpd1), "'x' must be numeric",
                class = "dplyr_error")
  expect_error( cor_rocky_mtn(set, table = "P") )
})

test_that("cor_rocky_mtn takes numeric values only", {

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  expect_error( cor_rocky_mtn(set, table = "P", weight, method = "X"),
                "'method' must be one of: pearson, kendall, spearman")
})
