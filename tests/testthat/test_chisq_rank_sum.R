context("micro_chisq and micro_rank_sum")

test_that("micro_chisq only takes correct modsum", {
  data(bpd_cla); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_cla), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  X <- rnorm(10)
  expect_error(
    micro_chisq(set, table = "P", grp_var = bpd1, mod = X),
    "'mod' must be output from either nb_mods or bb_mods")
})
