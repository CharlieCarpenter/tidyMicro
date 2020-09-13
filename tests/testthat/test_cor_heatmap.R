context("cor_heatmap")

test_that("cor_heatmap takes numeric values only", {

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  expect_error( cor_heatmap(set, table = "P", bpd1), "'x' must be numeric" )
  expect_error( cor_heatmap(set, table = "P") )
})

test_that("cor_heatmap takes numeric values only", {

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  expect_error( cor_heatmap(set, table = "P", weight, method = "X"),
                "'method' must be one of: pearson, kendall, spearman")
})
