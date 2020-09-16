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

test_that("no repeated taxa in Model_Coef", {
  data("mrsa_clin"); data("mrsa_gen")

  tidy.mrsa <- tidy_micro(mrsa_gen, mrsa_clin, "Genus", count_summary = F,
                          prev_cutoff = 1, ra_cutoff = 1,
                          exclude_taxa = c("Unclassified", "Bacteria"),
                          filter_summary = F)

  nb.int <- nb_mods(tidy.mrsa, table = "Genus", Aureus_Positive*Age)
  nb.int$Model_Coef$Taxa %<>% pull.lev(6)

  expect_error(nb_bars(nb.int, Aureus_Positive*Age, top_taxa = 10, quant_style = 'discrete'),
    "Repeated Taxa names exist in model's 'Model_Coef'")
})
