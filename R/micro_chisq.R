#' @title Run Chi-Squared tests for each taxa
#' @name micro_chisq
#' @description Run Chi-Squared tests for presence / absence of each taxa in you data set, or each taxa that didn't converge in negative binomial models
#' @param micro_set A tidy_micro data set
#' @param table The OTU table you'd like to test
#' @param grp_var Grouping variable for chi-squared test
#' @param y Response variable for chi-squared test. Default is presence / absence (bin)
#' @param mod The output from mods if you'd like to only run on taxa that did not converge
#' @param ... Options to be passed to chisq.test
#' @details If the taxa are present or absent in every subject the chi-sqared test will not but run. The returned chi-sqared stat will either be "All Absent" or "All Present." This will be clear in the output
#' @references \code{help(chisq.test)}
#' @return A data from containing the taxa, the chi-squared statistic, and the p-value of the test.
#' @examples
#' data(cla); data(clin)
#'
#' set <- tidy_micro(otu_tabs = cla, tab_names = "Class", clinical = clin,
#' prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Chi-squared test on every taxa's presence/absence
#' set %>% micro_chisq(table = "Class", grp_var = bpd1,
#' simulate.p.value = TRUE)
#'
#' ## Chi-squared test on every taxa whose model didn't converge
#' nb_cla <- set %>% nb_mods(table = "Class", bpd1)
#'
#' micro_chisq(micro_set = set, table = "Class", grp_var = bpd1,
#' mod = nb_cla, simulate.p.value = TRUE)
#' @export
micro_chisq <- function(micro_set, table, grp_var, y = bin, mod = NULL, ...){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  ## Run on every taxa
  if(is.null(mod)){
    dat <- null.mod.chisq(micro_set, table,
                          grp_var = !!rlang::enquo(grp_var),
                          y = !!rlang::enquo(y), ...)
  } else { ## Running on taxa that don't converge

    ## Pulling taxa that didn't converge
    if("FE_Converged" %in% names(mod$RA_Summary)){
      tax <- mod$RA_Summary %>%
        dplyr::filter(!(.data$FE_Converged)) %>%
        dplyr::pull(.data$Taxa)
    }
    if("RE_Converged" %in% names(mod$RA_Summary)){
      tax <- mod$RA_Summary %>%
        dplyr::filter(!(.data$RE_Converged)) %>%
        dplyr::pull(.data$Taxa)
    }

    dat <- mod.chisq(micro_set, table, grp_var = !!rlang::enquo(grp_var),
                     y = !!rlang::enquo(y), mod = mod, tax = tax, ...)
  }

  suppressWarnings(
    dat %>% dplyr::mutate(FDR_Pval = stats::p.adjust(.data$P_Val, method = "BH")) %>%
      Taxa_ord() %>%
      dplyr::arrange(.data$Taxa) %>%
      dplyr::select(-1))
}


null.mod.chisq <- function(micro_set, table, grp_var, y, ...){
  micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    plyr::ddply(~ .data$Taxa, function(set){

      ## Pulling grouping variable and y
      grp <- set %>% dplyr::pull(!!rlang::enquo(grp_var))
      rsp <- set %>% dplyr::pull(!!rlang::enquo(y))

      ## Running Chi-Squared test when appropriate
      if(dim(table(rsp)) < 2){

        X <- ifelse(sum(set$bin) == nrow(set), "All Present", "All Absent")

        dat <- data.frame(Taxa = unique(set$Taxa),
                          Chi_squared = X,
                          P_Val = NA)
      } else{
        chi <- stats::chisq.test(x = grp, y = rsp, ...)

        dat <- data.frame(Taxa = unique(set$Taxa),
                          Chi_squared = chi$statistic,
                          P_Val = chi$p.value)
      }

      dat
    } )
}

mod.chisq <- function(micro_set, table, grp_var, y = bin, mod, tax, ...){
  micro_set %>%
    ## Filtering out table and taxa that didn't converge
    dplyr::filter(.data$Table == table, .data$Taxa %in% tax) %>%
    plyr::ddply( ~ .data$Taxa, function(set){

      ## Pulling grouping variable and y
      grp <- set %>% dplyr::pull(!!rlang::enquo(grp_var))
      rsp <- set %>% dplyr::pull(!!rlang::enquo(y))

      ## Running Chi-Squared test when appropriate
      if(dim(table(rsp)) < 2){

        X <- ifelse(sum(set$bin) == nrow(set), "All Present", "All Absent")

        dat <- data.frame(Taxa = unique(set$Taxa),
                          Chi_squared = X,
                          P_Val = NA)
      } else{
        chi <- stats::chisq.test(x = grp, y = rsp, ...)

        dat <- data.frame(Taxa = unique(set$Taxa),
                          Chi_squared = chi$statistic,
                          P_Val = chi$p.value)
      }

      dat
    } )
}


