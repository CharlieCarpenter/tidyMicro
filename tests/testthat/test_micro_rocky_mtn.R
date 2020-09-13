context('micro_rocky_mtn')

test_that("micro_rocky_mtn requires specified covariate", {
  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  nb <- nb_mods(set, table = "P", bpd1)

  expect_error(micro_rocky_mtn(nb), "NB_RockMtn requires a model coefficient")
})

test_that("micro_rocky_mtn requires nb_ or bb_mods", {
  data(bpd_phy); data(bpd_clin)

  set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                    prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                    filter_summary = FALSE, count_summary = FALSE) %>%
    filter(day == 7)

  nb <- nb_mods(set, table = "P", bpd1)
  nb$Model_Type <- 'nn'

  expect_error(micro_rocky_mtn(nb, bpd1),
               "'modsum' must be output from either nb_mods or bb_mods")
})
