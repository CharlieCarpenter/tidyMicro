context('nb_bars')

test_that("nb_bars shouldn't run on bb_bars output", {
  skip_on_cran()

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  bb <- bb_mods(set, table = "P", bpd1)
  expect_error( nb_bars(bb, bpd1), "'modsum' should be the output from 'nb_mods'")
})
