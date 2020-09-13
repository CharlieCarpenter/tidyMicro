#' @title Fit negative binomial models to each taxa within an OTU table
#' @name nb_mods
#' @description Fit negative binomial models to each taxa within an OTU table through \code{\link[MASS]{glm.nb}} in the \pkg{MASS} package. Models can include a random effect if desired. Modesl will then be fit through \code{\link[lme4]{glmer.nb}} in the lmer package. Summaries for models or confidence intervals that fail to converge will not be returned, but taxa summaries will be provided in the output. Rank-Sum tests or presence/absence tests can be run on these taxa using \code{tidi_rank_sum} or \code{tidi_chisq}, respectively
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param ... Covariates of interest. Can be interactions such as Group*Age
#' @param Offset Logical; include subject sequencing depth as an offset for negative binomial models. This is highly recommended
#' @param ref A character vector of the desired reference levels for each factor covariate. The order of the specifed references must match the order for the corresponding covariates specified in '...'
#' @param SS_type Type of sums of squares calculated in \code{\link[car]{Anova}}. Either type II (2) or type III (3) sums of squares
#' @references \code{\link[car]{Anova}}, \code{\link[MASS]{glm.nb}}, \code{\link[lme4]{glmer.nb}}
#' @details Models containing only fixed effects are fit using \code{\link[MASS]{glm.nb}} in the \pkg{MASS} package and models containing random effects are fit using \code{\link[lme4]{glmer.nb}}. ANOVA / ANCOVA tests are conducted using a Likelihood Ratio test for fixed effects models and Chi-Squared tests for random effect models.
#' @return A list containing several different model components and summaries
#' \item{Convergend_Summary}{A data.frame of model summaries from convergent models. Includes the Taxa name, the model coefficient, the estimated beta, the beta's 95 percent confidence interval, Z score, p_value, false discovery rate p-value, and p-value from likelihood ratio test}
#' \item{Estimate_Summary}{A data.frame of model estimates from convergent models intended to be ready for export for publications. Includes the Taxa name, the model coefficient, the estimated Rate Ratio, the Wald 95 percent confidence interval, the Z-score, and false discovery rate p-value}
#' \item{RA_Summary}{A data.frame of taxa summaries. Includes the Taxa name, grouping variables (each factor variable in your models), sample size (n), percent of 0 counts, basic summaries of relative abundance, percentiles of relative abundance, and a logical indicator of whether or not the model converged}
#' \item{formula}{The formula used in the model}
#' \item{Model_Coef}{Model coefficients (used in plotting funcitons)}
#' \item{Model_Covs}{Model covariates (used in plotting functions)}
#' @note False Discovery Rate p-values are calculated using \code{\link[stats]{p.adjust}}. Estimated rate ratios and confidence intervals for interactions in the Estimate_Summary table include all main effects. It is not simply the exponentiated interaction beta, it is the interaction of the sum of the intercept, corresponding main effect betas, and interaction betas
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#'
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7)
#'
#' nb_fam <- set %>%
#' otu_filter(prev_cutoff = 5, ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' nb_mods(table = "Family", bpd1)
#'
#' names(nb_fam)
#' nb_fam$Estimate_Summary
#' @export
nb_mods <- function(micro_set, table, ..., Offset=TRUE, ref=NULL, SS_type=c(2,3,"II","III")){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  micro_set %<>% dplyr::filter(.data$Table == table) %>% ## Filtering out by table
    dplyr::mutate(Taxa = factor(.data$Taxa, levels = unique(.data$Taxa)))
  ## Making Taxa levels the Taxa names w/in table for tapply

  ## Defining formula for glm.nb with  w/ or w/o and Offset. Offset is the default

  f <- suppressWarnings(formula_fun(...))
  if(Offset) f <- paste(f, "+ offset(log(Total))")

  if(missing(SS_type)) SS_type <- 2
  if(SS_type %nin% c(2,3,"II","III")) stop("SS_type must be either 2, 3, 'II', or 'III' \n")

  Cov <- micro_set %>% dplyr::select(cov_str(...)) ## pulling out model covariates

  Convergent_Models <- micro_set %>%
    FE.con(f, SS_type, ...) ## Gives summaries of convergent models

  ## Gives Taxa summaries
  RA_Summary <- micro_set %>% N.con(.data$ra, ...) %>%
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
                   dplyr::filter(.data$FE_Converged==T) %>%
                   dplyr::select(.data$Taxa, .data$Coef, .data$Beta, .data$CI,
                                 .data$Z, .data$P_val, .data$FDR_Pval, .data$LRT),
                 Estimate_Summary=Convergent_Models %>%
                   dplyr::filter(.data$FE_Converged==T) %>%
                   dplyr::select(.data$Taxa, .data$Coef, .data$RR,
                                 .data$CI_95, .data$Z, .data$FDR_Pval) %>%
                   dplyr::filter(.data$Coef != "(Intercept)"),
                 RA_Summary=RA_Summary,  ## Nonconvergent summaries
                 formula=f,  ## The formula used so you can check it is what you intended
                 Model_Coef = Convergent_Models %>%
                   dplyr::filter(.data$FE_Converged==T) %>%
                   dplyr::select(.data$Taxa, .data$Coef, .data$Intercept,
                                 Estimate = .data$Beta, .data$Cov_Type),
                 Model_Covs = Cov,
                 Model_Type = "nb_mod"
  )

  result
}

## Code to relevel a factor variable (unnessessary)
# ## relevels the factor variables
# if(!is.null(ref)){
#   if(!is.character(ref)) stop("ref must be a character string of valid reference levels")
#   Cov <- cov_str(...) %>%
#     dplyr::select(micro_set,.)              ## pulling out model covariates
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


## Helper functions
# MK: I moved this source to the set-up code in tidy_MIBI
# source('~/Documents/Research/Current/Jed_Stuff/03_tidy_microbiome/Code/Summary_Helper.R')

# ## Fixed Effects Negative Binomial Models
# FE_mod <- function(d,f){
#
#   ## Models might converge, but have problematic fitting issues.
#   options(warn = 2) ## Changes warnings to errors for the try() statements
#
#   ## Pulling out convergent Taxa
#   conv <- unique(d$Taxa)[ !sapply(unique(d$Taxa), function(x,f){
#     try(glm.nb(as.formula(f), data = d %>%
#                  filter(Taxa == as.character(x)), maxit = 5000), silent = T) %>%
#       inherits("try-error")
#   }, f=f) ]
#
#   dc <- d %>% filter(Taxa %in% conv)
#
#   ## Repeating the process just to make sure the models converge
#   Converged <- unique(dc$Taxa)[ !sapply(unique(dc$Taxa), function(x,f){
#     try(glm.nb(as.formula(f), data = dc %>%
#                  filter(Taxa==as.character(x)), maxit=5000), silent = T) %>%
#       inherits("try-error")
#   }, f=f) ]
#
#   options(warn = 0) ## Making warnings warnings again
#
#   d %>% mutate(FE_Converged = Taxa %in% Converged) %>%
#     return()  ## Returns data filtered by table with new column. Mostly for piping into NB_summary
# }
#
# ## Random Effects Negative Binomial Models
# RE_mod <- function(d,f){
#   ## Pulling out convergent Taxa
#   conv <- unique(d$Taxa)[ !sapply(unique(d$Taxa), function(x,f){
#     try(glmer.nb(as.formula(f), data = d %>%
#                    filter(Taxa==as.character(x))),silent = T) %>%
#       inherits("try-error")
#   }, f=f) ]
#
#   dc <- d %>% filter(Taxa %in% conv)
#
#   ## Repeating the process just to make sure the models converge
#   Converged <- unique(dc$Taxa)[ !sapply(unique(dc$Taxa), function(x,f){
#     try(glmer.nb(as.formula(f), data = dc %>%
#                    filter(Taxa==as.character(x))),silent = T) %>%
#       inherits("try-error")
#   }, f=f) ]
#
#   d %>% mutate(RE_Converged = Taxa %in% Converged) %>%
#     return()  ## Returns data filtered by table with new column. Mostly for piping into NB_summary
# }

# #' @param RE A character string specifying the random effect term in the syntax of a lmer model; i.e. a random intercept for subject would be written as RE = "(1|Subject)"
## Random Effect Option
# ## Adding random effect
# f <- paste(f,"+", RE)
#
# Convergent_Models <- micro_set %>%
#   RE.con(f, SS_type, !!!rlang::quos(...), RE = RE) ## Gives summaries of convergent models
#
# ## Gives Taxa summaries
# RA_Summary <- micro_set %>% N.con(ra,!!!rlang::quos(...)) %>%
#   left_join(Convergent_Models %>%
#               dplyr::select(.data$Taxa, .data$RE_Converged), by = "Taxa") %>%
#   arrange(RE_Converged, Taxa) %>% unique
#
# ## TRUE / FALSE for taxa converging
# con_mod <- RA_Summary %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(.data$Taxa, .keep_all = T) %>%
#   dplyr::pull(.data$RE_Converged)
#
# cat("\n", sum(con_mod), " taxa converged\n", sep = "")
# cat(sum(!con_mod), "taxa did not converge\n")
#
# result <- list(Convergent_Summary = Convergent_Models %>%  ## Convergent model summaries
#                  dplyr::filter(.data$RE_Converged) %>%
#                  dplyr::select(.data$Taxa, .data$Coef, .data$Beta,
#                                .data$CI, .data$Z, .data$P_val, .data$FDR_Pval,
#                                .data$LRT),
#                Estimate_Summary = Convergent_Models %>%
#                  dplyr::filter(.data$RE_Converged) %>%
#                  dplyr::select(.data$Taxa, .data$Coef, .data$RR,
#                                .data$CI_95, .data$Z, .data$FDR_Pval),
#                RA_Summary = RA_Summary,  ## Nonconvergent summaries
#                formula = f,  ## The formula used so you can check it is what you intended
#                Model_Coef = Convergent_Models %>%
#                  dplyr::select(.data$Taxa, .data$Coef, .data$Intercept,
#                                Estimate = .data$Beta, .data$Cov_Type),
#                Model_Covs = .data$Cov
# )

