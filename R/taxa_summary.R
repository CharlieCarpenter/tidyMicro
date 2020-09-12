#' @title Summarize the information
#' @name taxa_summary
#' @description Give taxa summary table stratified by variables of interest and/or OTU tables
#' @param micro_set A tidy_micro data set
#' @param ... Covariates of interest
#' @param table OTU table of interest. If NULL, all tables will be used
#' @param obj The taxonomic information of interest
#' @param taxa Logical; Whether or not to stratify by taxa
#' @return A tibble containing columns of stratifying variables and several summary columns
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#'
#' otu_tabs <- list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' mutate(bpd1 = factor(bpd1))
#'
#' ## Summarize each taxa by Table
#' set %>% taxa_summary
#'
#' ## Summarize each taxa by a categorical variable of interest
#' set %>% taxa_summary(bpd1)
#'
#' ## Summarize each taxa by a categorical variable of interest within a Table
#' set %>% taxa_summary(bpd1, table = "Phylum")
#'
#' ## Summarize within group or table only
#' set %>% taxa_summary(taxa = FALSE)
#' @export
taxa_summary <- function(micro_set, ..., table = NULL, obj = ra, taxa = TRUE){

  if(!is.null(table)){
    if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")
    micro_set %<>% dplyr::filter(.data$Table == table) ## Filtering out table
  }

  ## Making "Taxa" a character so that it won't be grouped by below
  ## Option for when you don't want it
  if(!taxa) micro_set %<>% dplyr::mutate(Taxa = as.character(.data$Taxa))

  micro_set %>%
    dplyr::select(!!rlang::enquo(obj), .data$Table, cov_str(...), .data$Taxa) %>%
    dplyr::group_by_if(is.factor) %>%  ## Groups by factor covariates when given
    taxa_summarize(!!rlang::enquo(obj))
}
