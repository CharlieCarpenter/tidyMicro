context('pca_3d')

test_that("warning message pops",{
  data(bpd_phy); data(bpd_clin)

  expect_warning(
    tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
               prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
               filter_summary = FALSE, count_summary = FALSE) %>%
      pca_3d(table = "P", time_var = day, subject = study_id)
  )
})

test_that("Requires 'time_var' and 'subject'",{
  data(bpd_phy); data(bpd_clin)

  ts <- tidy_micro(otu_tabs = list("P" = bpd_phy), clinical = bpd_clin,
                   prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = "Bacteria",
                   filter_summary = FALSE, count_summary = FALSE)

  expect_error(pca_3d(ts, table = "P", time_var = day),
               "Must supply a column name for subject ID")

  expect_error(pca_3d(ts, table = "P", time_var = day, subject = X),
               "subject is not a column name in supplied micro_set")

  expect_error(pca_3d(ts, table = "P", subject = study_id),
               "Must supply a column name for the time variable")

  expect_error(pca_3d(ts, table = "P", time_var = X, subject = study_id),
               "time_var is not a column name in supplied micro_set")
})
