#' @title Fit beta binomial models to each taxa within an OTU table
#' @name bb_mods
#' @description Fit beta binomial models to each taxa within an OTU table through \code{\link[VGAM]{vglm}} in the \pkg{VGAM} package. Summaries for models or confidence intervals that fail to converge will not be returned, but taxa summaries will be provided in the output. Rank-Sum tests or presence/absence tests can be run on these taxa using \code{tidi_rank_sum} or \code{tidi_chisq}, respectively
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param ... Covariates of interest. Can be interactions such as Group*Age
#' @param CI_method Character indicating the type of method used for confidence interval estimation. Wald intervals are the current default. Abbreviations allowed. See \code{\link[VGAM]{confintvglm}} for more details
#' @param SS_type Type of sums of squares calculated in \code{\link[VGAM]{anova.vglm}}. Either type II (2) or type III (3) sums of squares. Type II is the default
#' @param trace Print messages of model fitting proceedure
#' @references \code{\link[VGAM]{anova.vglm}}, \code{\link[VGAM]{vglm}}, \code{\link[VGAM]{betabinomial}}
#' @details Models containing only fixed effects are fit using \code{\link[VGAM]{vglm}} in the \pkg{VGAM} package. ANOVA / ANCOVA tests are conducted using a Likelihood Ratio test
#' @return A list containing several different model components and summaries
#' \item{Convergend_Summary}{A data.frame of model summaries from convergent models. Includes the Taxa name, the model coefficient, the estimated beta, the beta's 95 percent confidence interval, Z score, p_value, false discovery rate p-value, and p-value from likelihood ratio test}
#' \item{Estimate_Summary}{A data.frame of model estimates from convergent models intended to be ready for export for publications. Includes the Taxa name, the model coefficient, the estimated Rate Ratio, the Wald 95 percent confidence interval, the Z-score, and false discovery rate p-value}
#' \item{RA_Summary}{A data.frame of taxa summaries. Includes the Taxa name, grouping variables (each factor variable in your models), sample size (n), percent of 0 counts, basic summaries of relative abundance, percentiles of relative abundance, and a logical indicator of whether or not the model converged}
#' \item{formula}{The formula used in the model}
#' \item{Model_Coef}{Model coefficients (used in plotting funcitons)}
#' \item{Model_Covs}{Model covariates (used in plotting functions)}
#' @note False Discovery Rate p-values are calculated using \code{\link[stats]{p.adjust}}. Estimated rate ratios and confidence intervals for interactions in the Estimate_Summary table include all main effects. It is not simply the exponentiated interaction beta, it is the interaction of the sum of the intercept, corresponding main effect betas, and interaction betas
#' @import VGAM
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#'   filter(day == 7) ## Only including first week
#'
#' \donttest{
#' bb_phy <- set %>%
#'
#' ## Filtering out low abundance and unclassified taxa
#' ## These models will either break or we don't care about them
#' otu_filter(prev_cutoff = 5, ra_cutoff = 0.1,
#'            exclude_taxa = c("Unclassified", "Bacteria")) %>%
#'
#' ## Beta binomial models for each Family of taxa with bpd1 as a covariate
#' bb_mods(table = "Phylum", bpd1, CI_method = "wald")
#'
#' names(bb_phy)
#' bb_phy$Estimate_Summary
#' }
#' @export
bb_mods <- function(micro_set, table, ..., CI_method = c("wald", "profile"),
                    SS_type = c(2, 3, "II", "III"), trace = FALSE){
  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

    micro_set %<>% dplyr::filter(.data$Table == table) %>% ## Filtering out by table
      dplyr::mutate(Taxa = factor(.data$Taxa, levels = unique(.data$Taxa)))
  ## Making Taxa levels the Taxa names w/in table for tapply

  ## Defining formula for glm.nb with  w/ or w/o and Offset. Offset is the default

  f <- suppressWarnings(formula_fun_bb(...))

  if(missing(CI_method)) CI_method <- "wald"
  if(grepl("[profi]", CI_method)) warning("Profile likelihood confidence intervals greatly increase computation time.\n")
  if(missing(SS_type)) SS_type <- 2
  if(SS_type %nin% c(2,3,"II","III")) stop("SS_type must be either 2, 3, 'II', or 'III' \n")

  Convergent_Models <- micro_set %>% ## Gives summaries of convergent models
    plyr::ddply(~ .data$Taxa, FE_bb, f. = f, trace = trace,
                method. = CI_method, SS_type. = SS_type, ...)

  if(nrow(Convergent_Models) == 0) stop("No taxa models converged.\n")

  Cov <- micro_set %>% dplyr::select(cov_str(...)) ## pulling out model covariates

  ## False Discovery Rate adjustment
  Convergent_Models %<>%
    dplyr::mutate(FDR_Pval = stats::p.adjust(.data$P_val, method = "BH") %>%
                    round(4))

  ## Gives Taxa summaries
  RA_Summary <- micro_set %>% N.con(.data$ra,...) %>%
    dplyr::left_join(Convergent_Models %>%
                dplyr::select(.data$Taxa, .data$FE_Converged), by = "Taxa") %>%
    dplyr::arrange(.data$FE_Converged, .data$Taxa) %>% unique

  ## TRUE / FALSE for taxa converging
  con_mod <- RA_Summary %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$Taxa, .keep_all = T) %>%
    dplyr::pull(.data$FE_Converged)

  message("\n", sum(con_mod), " taxa converged.")
  message(sum(!con_mod), " taxa did not converge.")

  result <- list(Convergent_Summary=Convergent_Models %>%  ## Convergent model summaries
                   dplyr::filter(.data$FE_Converged) %>%
                   dplyr::select(.data$Taxa, .data$Coef, .data$Beta,
                                 .data$CI, .data$LRT, .data$P_val, .data$FDR_Pval),
                 Estimate_Summary=Convergent_Models %>%
                   dplyr::filter(.data$FE_Converged) %>%
                   dplyr::select(.data$Taxa, .data$Coef, .data$OR,
                                 .data$CI_95, .data$LRT, .data$FDR_Pval) %>%
                   dplyr::filter(!grepl("(Intercept)", .data$Coef)),
                 RA_Summary=RA_Summary,  ## Nonconvergent summaries
                 formula=f,  ## The formula used so you can check it is what you intended
                 Model_Coef = Convergent_Models %>%
                   dplyr::filter(.data$FE_Converged) %>%
                   dplyr::select(.data$Taxa, .data$Coef, .data$Intercept,
                                 Estimate = .data$Beta, .data$Cov_Type),
                 Model_Covs = Cov,
                 Model_Type = "bb_mod"
  )

  result
}

## Code to relevel a factor variable (unnessessary)
# ## relevels the factor variables
# if(!is.null(ref)){
#   if(!is.character(ref)) stop("ref must be a character string of valid reference levels")
#   Cov <- cov_str(...) %>%
#     dplyr::select(micro_set,rlang::.data$.)              ## pulling out model covariates
#
#   Fac <- Cov[sapply(Cov,class) == "factor"]   ## selecting factor variables
#
#   ## for loops are bad, but this should happen quickly (no need to apply)
#   ## relevels each factor
#   for(i in seq(1,length(ref))){
#     Fac[i] <- relevel(Fac[,i], ref = ref[i])
#   }
#
#   Cov <- data.frame(Fac, Cov[sapply(Cov,class) != "factor"])
#
#   micro_set <- data.frame(
#     micro_set[which( !(names(micro_set) %in% names(Cov)) )],
#     Cov)
# } else {
#   Cov <- cov_str(...) %>%
#     dplyr::select(micro_set,.)              ## pulling out model covariates
# }


