context("bb_bars")

test_that("bb_bars knows it is a bb_mod",{

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  nb <- nb_mods(set, table = "P", bpd1)

  expect_error( bb_bars(nb, bpd1))

})
