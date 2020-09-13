#' @title Create forest plots from negative binomial taxa models
#' @name micro_forest
#' @description Create forest plots for specified coefficients in negative binomial taxa models. Plots estimated beta coefficients and confidence intervals
#' @param modsum The output from nb_mods
#' @param ... The covariate you'd like to plot. Must be in the models created by nb_mods
#' @param main The title for your plot
#' @param ylab The label for the y-axis; default is "Taxa"
#' @param xlab The label for the x-axis; default is output from function "TeX"
#' @param subtitle The plot subtitle
#' @param legend_title The title of the plot's legend
#' @param legend_labs The names of the elements within the legend
#' @return Returns a ggplot that you can add geoms to if you'd like
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Creating negative binomial models on filtered tidi_micro set
#' nb_fam <- set %>%
#' otu_filter(prev_cutoff = 5, ra_cutoff = 0.1,
#' exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' nb_mods(table = "Family", bpd1)
#'
#' nb_fam %>% micro_forest(bpd1)
#' @export
micro_forest <- function(modsum, ..., main, ylab, xlab, subtitle,
                         legend_title, legend_labs){

  if(modsum$Model_Type %nin% c('bb_mod', 'nb_mod')){
    stop("'mod' must be output from either nb_mods or bb_mods")
  }

  CC <- modsum$Convergent_Summary %>%
    dplyr::mutate(CI.u = stringr::str_sub(.data$CI, start = 2L,
                                          end = stringr::str_length(.data$CI)-1) %>%
             stringr::str_split(", ") %>% purrr::map(2) %>% unlist %>% as.numeric,
           CI.l = stringr::str_sub(.data$CI, start = 2L,
                                   end = stringr::str_length(.data$CI)-1) %>%
             stringr::str_split(", ") %>% purrr::map(1) %>% unlist %>% as.numeric
           )

  if(any(is.na(CC$CI.u)) | any(is.na(CC$CI.l))){
    warning("Some confidence intervals contain NAs and will be removed.")
    CC %<>% dplyr::filter(!is.na(.data$CI.u), !is.na(.data$CI.l))
  }

  if(missing(main)) main <- NULL; if(missing(ylab)) ylab <- "Taxa"
  if(missing(xlab)) xlab <- latex2exp::TeX("$\\hat{\\beta}$ with 95% CI")
  if(missing(subtitle)) subtitle <- NULL
  if(missing(legend_title)) legend_title <- "Model_Covs"

  if(any((CC$CI.u - CC$CI.l) > 10)){
    warning("At least one Confidence interval width is over 10. These estimates may be unstable.")
  }

  if(!missing(...)){
    if(length(cov_str(...)) == 1) CC %<>% dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...)))
    if(length(cov_str(...)) == 2) {
      CC %<>% dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                              stringr::str_detect(.data$Coef, cov_str(...)[2]))
    }
  }

  gg <- ggplot2::ggplot(CC, ggplot2::aes(y = .data$Beta,
                                         x = interaction(.data$Coef,.data$Taxa) %>% sort,
                        col = .data$Coef, ymin = .data$CI.l, ymax = .data$CI.u)) +
    ggplot2::scale_x_discrete(labels = CC$Taxa) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(yintercept=0, lty=3) +
    ggplot2::scale_y_continuous(breaks = sort(round(c(seq(min(CC$CI.l), max(CC$CI.u), length.out=5), 0)))) +
    ggplot2::coord_flip() +
    ggplot2::labs(y = xlab, x = ylab, title = main, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::guides(col = ggplot2::guide_legend(title = legend_title))

    if(!missing(legend_labs)){
      gg <- gg + ggplot2::scale_colour_discrete(labels = legend_labs)
    }

  gg
}
