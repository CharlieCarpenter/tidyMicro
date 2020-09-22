###############################
##
## Project: tidyMicro
##
## Purpose: Helper functions for pca_3d
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-01-23
##
## ---------------------------
## Notes:
##
##
## ---------------------------

## Helper functions to get data in the correct format
A_matricization <- function(micro_set, table, subject, y, time_var){

  l_otu <- micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::select(.data$Taxa, sub = !!rlang::enquo(subject),
                  time = !!rlang::enquo(time_var),
                  !!rlang::enquo(y)) %>%
    # Just in case there are extra levels of the time variable
    # that aren't in the data set
    dplyr::mutate(time = factor(.data$time, levels = sort(unique(.data$time))))

  # 'empty' data frame to fill
  wide_otu <- data.frame(sub = unique(l_otu$sub))
  # for loop to 'spread' taxa out by timepoint
  for(i in seq(1,length(levels(l_otu$time)))){

    ## Filter out by time point, and uniting time_taxa vars into one unique column
    ss <- l_otu %>%
      dplyr::filter(.data$time == levels(.data$time)[i]) %>%
      tidyr::unite("Time_Taxa", .data$time, .data$Taxa)

    ## Filtering out subjects that aren't consistent b/t time points
    wide_otu %<>% dplyr::filter(.data$sub %in% ss$sub)
    ss %<>% dplyr::filter(.data$sub %in% wide_otu$sub)

    ## Spreading counts into wide format
    ss %<>% tidyr::pivot_wider(names_from = "Time_Taxa", values_from = !!rlang::enquo(y))

    ## Tacking on time point counts to join by subject
    wide_otu %<>% dplyr::full_join(ss, by = "sub")

  }

  ## Checking that subjects are consistent across all time points
  if(!all(unique(l_otu$sub) %in% unique(wide_otu$sub))){
    warning("Subjects are not consistent across time points.\nOnly complete cases will be used.",
            call. = FALSE)
  }

  ## making subjects rownames
  wide_otu %>%
    tibble::column_to_rownames(var = "sub")
}
