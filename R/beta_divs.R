#' @title Beta Diversity Calculations for tidy_micro
#' @name beta_div
#' @description Calculate beta diversities of your tidy_micro set. This function reformats the data into the original OTU table and then feeds that into the vegdist function
#' @param micro_set A tidy_micro data set
#' @param table Table you'd like to use when calculating alpha diversity. Your lowest level is recommended
#' @param method A dissimilarity method compatible with \code{\link[vegan]{vegdist}}
#' @references \code{\link[vegan]{vegdist}}
#' @return A symmetrix distance matrix
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Bray-Curtis beta diversity
#' bray <- set %>% beta_div(table = "Family")
#'
#' ## Morisita-Horn beta diversity
#' horn <- set %>% beta_div(table = "Family", method = "horn")
#' @export
beta_div <- function(micro_set, table, method = "bray"){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  ## Getting original OTU table
  beta_div <- micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::select(.data$cts, .data$Lib, .data$Taxa) %>%
    tidyr::pivot_wider(names_from = .data$Taxa, values_from = .data$cts)

  ## Storing Library names
  Lib <- beta_div$Lib

  ## Calculating Beta divs
  beta_div %<>%
    dplyr::select(-.data$Lib) %>%
    vegan::vegdist(method = method) %>%
    as.matrix

  ## restoring rownames
  rownames(beta_div) <- Lib

  beta_div
}
