#' @title Run rank sum tests for each taxa within an OTU table
#' @name micro_rank_sum
#' @description Runs a rank sum test for each taxa within an OTU table or each taxa that didn't converge in \code{\link[tidyMicro]{nb_mods}} or \code{\link[tidyMicro]{bb_mods}}
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param grp_var A factor variable for grouping
#' @param y A continuous response variable. Taxa relative abundance (ra) is recommended
#' @param mod The output from \code{\link[tidyMicro]{nb_mods}} or \code{\link[tidyMicro]{bb_mods}} if desired
#' @param ... Options to be passed to \code{\link[stats]{wilcox.test}} or \code{\link[stats]{kruskal.test}}
#' @references \code{\link[stats]{kruskal.test}} and \code{\link[stats]{wilcox.test}}
#' @details The grp_var must have a least 2 levels. For a 2 level factor a Mann-Whitney test will be calculated through \code{\link[stats]{wilcox.test}}, and for 3 or more levels a Kruskal-Wallis test will be run throuh \code{\link[stats]{kruskal.test}}
#' @return A data frame containing the p-value for each taxa's rank sum test.
#' @examples
#' data(bpd_cla); data(bpd_clin)
#'
#' set <- tidy_micro(otu_tabs = bpd_cla, tab_names = "Class", clinical = bpd_clin,
#' prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Rank sum test on every taxa's relative abundance
#' set %>% micro_rank_sum(table = "Class", grp_var = bpd1)
#'
#' ## Rank sum test on every taxa whose model didn't converge
#' nb_cla <- nb_mods(set, table = "Class", bpd1)
#'
#' micro_rank_sum(micro_set = set, table = "Class",
#' grp_var = bpd1, mod = nb_cla)
#' @export
micro_rank_sum <- function(micro_set, table, grp_var, y = ra, mod = NULL, ...){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")
  nlev <- nlevels(micro_set %>% dplyr::pull(!!rlang::enquo(grp_var)))

  if(!is.null(mod)){
    if(mod$Model_Type %nin% c('bb_mod', 'nb_mod')){
      stop("'mod' must be output from either nb_mods or bb_mods")
    }
    if('Model_Type' %nin% names(mod)){
      stop("'mod' must be output from either nb_mods or bb_mods")
    }
  }

  ## Kruskal-Wallis rank sum
  if(nlev > 2){
    message("Running Kruskal-Wallis rank sum test.")

    if(is.null(mod)) {
      dat <- kruskal_rank_sum_no_mod(micro_set, table,
                                     !!rlang::enquo(grp_var), !!rlang::enquo(y), ...)
    }
    else{
      dat <- kruskal_rank_sum_mod(micro_set, table,
                                  !!rlang::enquo(grp_var), !!rlang::enquo(y), mod, ...)
    }

    ## Wilcoxon for 2 groups
  } else if(nlev == 2){

    message("Running Wilcoxon rank sum test.")
    if(is.null(mod)){
      dat <- wilcox_rank_sum_no_mod(micro_set, table,
                                    !!rlang::enquo(grp_var), !!rlang::enquo(y), ...)
    }
    else{
      dat <- wilcox_rank_sum_mod(micro_set, table,
                                 !!rlang::enquo(grp_var), !!rlang::enquo(y), ...)
    }

  } else stop("grp_var must have two or more levels")

  suppressWarnings(
    dat %>% dplyr::mutate(FDR_Pval = stats::p.adjust(.data$P_Val, method = "BH")) %>%
      Taxa_ord() %>%
      dplyr::arrange(.data$Taxa) %>%
      dplyr::select(-1))
}

kruskal_rank_sum_no_mod <- function(micro_set, table, grp_var, y, ...){
  ## Run on every taxa
  micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    plyr::ddply(~ .data$Taxa, function(set){

      ## Pulling grouping variable and y
      grp <- set %>% dplyr::pull(!!rlang::enquo(grp_var))
      rsp <- set %>% dplyr::pull(!!rlang::enquo(y))

      if(length(unique(rsp)) == 1){
        X <- ifelse(all(rsp == 0), "All Absent", "All Present")

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = X)
      } else{
        ##running kruskal-wallis rank sum test
        kwt <- stats::kruskal.test(rsp ~ grp)

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = kwt$p.value)
      }

      dat
    } )
}

kruskal_rank_sum_mod <- function(micro_set, table, grp_var, y, mod, ...){
  ## Pulling taxa that didn't converge
  if("FE_Converged" %in% names(mod$RA_Summary)){
    if(all(mod$RA_Summary$FE_Converged)) stop("No tests run. All taxa models converged")

    tax <- mod$RA_Summary %>% dplyr::filter(!(.data$FE_Converged)) %>% dplyr::pull(.data$Taxa)
  }
  if("RE_Converged" %in% names(mod$RA_Summary)){
    if(all(mod$RA_Summary$RE_Converged)) stop("No tests run. All taxa models converged")

    tax <- mod$RA_Summary %>% dplyr::filter(!(.data$RE_Converged)) %>% dplyr::pull(.data$Taxa)
  }

  dat <- micro_set %>%
    ## Filtering out rank and taxa that didn't converge
    dplyr::filter(.data$Table == table, .data$Taxa %in% tax) %>%
    plyr::ddply(~ .data$Taxa, function(set){

      ## Pulling grouping variable and y
      grp <- set %>% dplyr::pull(!!rlang::enquo(grp_var))
      rsp <- set %>% dplyr::pull(!!rlang::enquo(y))

      if(length(unique(rsp)) == 1){
        X <- ifelse(all(rsp == 0), "All Absent", "All Present")

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = X)
      } else{
        ##running kruskal-wallis rank sum test
        kwt <- stats::kruskal.test(rsp ~ grp)

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = kwt$p.value)
      }

      dat
    } )

  dat
}

wilcox_rank_sum_no_mod <- function(micro_set, table, grp_var, y, ...){
  micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    plyr::ddply( ~ .data$Taxa, function(set){

      ## Pulling grouping variable and splitting based on that variable
      grp <- set %>% dplyr::pull(!!rlang::enquo(grp_var))
      spl <- split(set, grp)

      rsp <- c( dplyr::pull(spl[[1]], !!rlang::enquo(y)),
                dplyr::pull(spl[[2]], !!rlang::enquo(y)) )

      if(length(unique(rsp)) == 1){
        X <- ifelse(all(rsp == 0), "All Absent", "All Present")

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = X)
      } else{
        ## Running Wilcoxon / Mann-Whitney test
        wct <- stats::wilcox.test(spl[[1]] %>% dplyr::pull(!!rlang::enquo(y)),

                                  spl[[2]] %>% dplyr::pull(!!rlang::enquo(y)),

                                  ...) %>% broom::tidy()

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = wct$p.value)
      }

      dat
    } )
}

wilcox_rank_sum_mod <- function(micro_set, table, grp_var, y, mod, ...){
  ## Pulling taxa that didn't converge
  if("FE_Converged" %in% names(mod$RA_Summary)){
    if(all(mod$Convergent_Summary$FE_Converged)) stop("No tests run. All taxa models converged")

    tax <- mod$RA_Summary %>% dplyr::filter(!(.data$FE_Converged)) %>% dplyr::pull(.data$Taxa)
  }
  if("RE_Converged" %in% names(mod$RA_Summary)){
    if(all(mod$Convergent_Summary$RE_Converged)) stop("No tests run. All taxa models converged")

    tax <- mod$RA_Summary %>% dplyr::filter(!(.data$RE_Converged)) %>% dplyr::pull(.data$Taxa)
  }

  dat <- micro_set %>%
    ## Filtering out table and taxa that didn't converge
    dplyr::filter(.data$Table == table, .data$Taxa %in% tax) %>%
    plyr::ddply( ~ .data$Taxa, function(set){

      ## Pulling grouping variable and splitting based on that variable
      grp <- set %>% dplyr::pull(!!rlang::enquo(grp_var))
      rsp <- set %>% dplyr::pull(!!rlang::enquo(y))
      spl <- split(set, grp)

      if(length(unique(rsp)) == 1){
        X <- ifelse(all(rsp == 0), "All Absent", "All Present")

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = X)
      } else{
        ## Running Wilcoxon / Mann-Whitney test
        wct <- stats::wilcox.test(spl[[1]] %>% dplyr::pull(!!rlang::enquo(y)),
                                  spl[[2]] %>% dplyr::pull(!!rlang::enquo(y)),
                                  ...) %>% broom::tidy

        dat <- data.frame(Taxa = unique(set$Taxa),
                          P_Val = wct$p.value)
      }

    } )
  dat
}

