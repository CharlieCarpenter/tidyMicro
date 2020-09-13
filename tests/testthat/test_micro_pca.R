context("micro_pca")

test_that("micro_pca needs grp_var",{
  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  expect_error( micro_pca(set, table = "P"), "'grp_var' must be specified" )
})

test_that("micro_pca distance matrix",{
  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  nn <- nrow(beta_div(set, table = "P"))
  X <- data.frame(x1 = rnorm(nn), x2 = rpois(nn, 1))

  expect_error( micro_pca(set, table = "P", dist = X, grp_var = bpd1),
                "'dist' must be a distance matrix. Usually created from 'beta_div' or a similar function" )
})

