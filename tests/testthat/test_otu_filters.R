context('filter_funs')


test_that('filter limits work',{
  context("cor_heatmap")

    data(bpd_phy); data(bpd_clin)

    set <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                      prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                      filter_summary = FALSE, count_summary = FALSE) %>%
      filter(day == 7)

    ## Prev
    expect_error(tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                            prev_cutoff = -5),
                 "Prevalence cutoff must be between 0 and 100")
    expect_error(tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                            prev_cutoff = 105),
                 "Prevalence cutoff must be between 0 and 100")

    expect_error(otu_filter(set, prev_cutoff = -5),
                 "Prevalence cutoff must be between 0 and 100")
    expect_error(otu_filter(set, prev_cutoff = 105),
                 "Prevalence cutoff must be between 0 and 100")

    ## Prev
    expect_error(tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                            ra_cutoff = -5),
                 "Relative abundance cutoff must be between 0 and 100")
    expect_error(tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                            ra_cutoff = 105),
                 "Relative abundance cutoff must be between 0 and 100")

    expect_error(otu_filter(set, ra_cutoff = -5),
                 "Relative abundance cutoff must be between 0 and 100")
    expect_error(otu_filter(set, ra_cutoff = 105),
                 "Relative abundance cutoff must be between 0 and 100")

    ## Exclude
    expect_error(tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                            exclude_taxa = -5),
                 "'exclue_taxa' should be a character or left as NULL")
    expect_error(tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                            exclude_taxa = TRUE),
                 "'exclue_taxa' should be a character or left as NULL")

    expect_error(otu_filter(set, exclude_taxa = -5),
                 "'exclue_taxa' should be a character or left as NULL")
    expect_error(otu_filter(set, exclude_taxa = TRUE),
                 "'exclue_taxa' should be a character or left as NULL")
})
