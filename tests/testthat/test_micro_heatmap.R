context("micro_heatmap")

test_that("micro_heatmap needs correct modsum", {
  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  nb <- nb_mods(set, table = "P", bpd1)
  nb$Model_Type <- "X"

  expect_error( micro_heatmap(nb), "'modsum' must be output from either nb_mods or bb_mods" )
})
