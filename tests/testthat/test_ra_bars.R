context('ra_bars')

test_that("ra_bars can't have both top_taxa and RA",{
  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  expect_error(ra_bars(set, table = "P", bpd1, top_taxa = 10, RA = 10),
               "Can not aggregate based on both 'top_taxa' and 'RA'")
})
