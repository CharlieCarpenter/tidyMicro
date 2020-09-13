#' @title Linear regression on alpha diversities within a micro_set
#' @name micro_alpha_reg
#' @description A simple wrapper to run standard linear regression though the \code{\link[stats]{lm}} function. Will only use alpha diversities distinct libraries (Lib) from the specified table as to not inflate the sample size
#' @param alpha_set A tidy_micro data set with alpha diversities calculated by \code{alpha_div}
#' @param table OTU table of interest
#' @param ... Covariates of interest. Can include interaction terms such as \code{Group*Age}
#' @return A data frame containing the model estimates for each alpha diversity
#' @note Be aware of your minimal sequencing depth as this will be the size of all bootstrapped resamples (rarefied).
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including first week
#'
#' \donttest{
#' set_fam_alpha <- set %>% alpha_div(table = "Family", min_depth = 5000, min_goods = 90)
#' set_fam_alpha %>% micro_alpha_reg(table = "Family", bpd1)
#' }
#' @export
micro_alpha_reg <- function(alpha_set, table, ...){

  if(table %nin% unique(alpha_set$Table)) stop("Specified table is not in supplied micro_set")

  alpha_set %<>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::distinct(.data$Lib, .keep_all = TRUE)

  alphas <- c("Sobs", "Chao1", "ShannonE", "ShannonH", "SimpsonD")
  alpha_tab <- NULL
  error_tab <- NULL

  alpha_set %<>% dplyr::distinct(.data$Lib, .keep_all = TRUE)

  for(i in seq(1,length(alphas))){

    f <- suppressWarnings(
      stats::as.formula(paste0(alphas[i], "~", adonis_formula(...)))
    )

    mod <- stats::lm(f, data = alpha_set)
    CI <- stats::confint(mod) %>% round(4)

    alpha_tab <- rbind(alpha_tab,
                       cbind(Alpha_Div = alphas[i],
                             mod %>% broom::tidy(),
                             CI_95 = paste0("(", CI[,1], ", ", CI[,2], ")"))
    )
  }

  alpha_tab %>% dplyr::rename(Coef = .data$term, Beta = .data$estimate,
                              t.stat = .data$statistic)
}


## Option to include random effect
# else{
#   if(!is.character(RE)) stop("RE must be written as a character string")
#
#   alpha_set %<>% filter(Taxa == unique(Taxa)[1])
#
#
#   for(i in seq(1,length(alphas))){
#
#
#     f <- suppressWarnings(
#       as.formula(
#         paste(
#           paste0(alphas[i], "~", adonis_formula(...)),
#           RE, sep = "+"
#         )
#       )
#     )
#
#     mod <- lme4::lmer(f, data = alpha_set)
#
#     CI <- stats::confint(mod) %>% round(4) %>%
#       as.data.frame() %>% tibble::rownames_to_column()
#
#     fe <- mod %>% broom::tidy() %>%
#       dplyr::filter(group == "fixed") %>%
#       dplyr::select(-group)
#
#     fe_CI <- CI %>% filter(rowname %in% fe$term) %>% dplyr::select(-rowname)
#     er_CI <- CI %>% filter(rowname %nin% fe$term) %>% dplyr::select(-rowname)
#
#     alpha_tab <- rbind(alpha_tab,
#                        cbind(Alpha_Div = alphas[i],
#                              fe,
#                              CI_95 = paste0("(", fe_CI[,1], ", ", fe_CI[,2], ")"))
#     )
#
#     error_tab <- rbind(error_tab,
#                        cbind(Alpha_Div = alphas[i],
#                              mod %>% broom::tidy() %>%
#                                dplyr::filter(group != "fixed") %>%
#                                dplyr::select(term, estimate, group),
#                              CI_95 = paste0("(", er_CI[,1], ", ", er_CI[,2], ")")
#                        )
#     )
#   }
#
#   alpha_tab %<>% dplyr::rename(Coef = term, Beta = estimate, t.stat = statistic)
#   error_tab %<>% dplyr::rename(Sigma = estimate)
#
#   alpha_tab <- list(FE = alpha_tab, RE = error_tab)
# }

