context('beta_heatmap')

test_that('beta_heat map breaks with more than 1 variable',{

  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  bd <- beta_div(set, table = "P")
  expect_error( beta_heatmap(bp, set, bpd1, gender), "Must use one factor variable." )
})
