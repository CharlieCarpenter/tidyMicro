#' @title A function to run PERMANOVA on tidi_micro data sets
#' @name micro_PERMANOVA
#' @description A wrapper function to call \code{\link[vegan]{adonis2}} from the \code{vegan} package. PERMANOVA is a method for partitioning distance matrices among sources of variation and fitting linear models (e.g., factors, polynomial regression) to distance matrices; uses a permutation test with pseudo-F ratios
#' @param micro_set A tidy_micro data set
#' @param beta_div A dissimilarity matrix calculated by \code{beta_div}
#' @param method A character string indicating the method used to calculated dissimilarity
#' @param ... Covariates of interest
#' @param nperm Number of permutations
#' @details The function adonis2 is based on the principles of McArdle & Anderson (2001) and can perform sequential, marginal and overall tests. Function adonis2 also allows using additive constants or squareroot of dissimilarities to avoid negative eigenvalues
#' @references \code{\link[vegan]{vegdist}} \code{\link[vegan]{adonis2}}
#' @seealso \code{\link[vegan]{adonis}}
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Bray-Curtis beta diversity
#' bray <- set %>% beta_div(table = "Family")
#'
#' set %>% micro_PERMANOVA(bray, method = "bray", bpd1)
#' @export
micro_PERMANOVA <- function(micro_set, beta_div, method, ..., nperm = 999){
  micro_set %<>%
    dplyr::filter(.data$Lib %in% rownames(beta_div)) %>%
    dplyr::distinct(.data$Lib, .keep_all = T) %>%
    dplyr::arrange(.data$Lib)

  f <- paste("beta_div ~", suppressWarnings(adonis_formula(...)) ) %>%
    stats::as.formula()

  vegan::adonis2(f, data = micro_set, method = method, permutations = nperm)
}
