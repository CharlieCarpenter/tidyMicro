#' @title Create heatmaps of estiamted coefficients from negative binomial models
#' @name micro_heatmap
#' @description A function to create heatmaps of estimated beta coeffients from each model fit by nb_mods
#' @param modsum The output from nb_mods
#' @param low_grad The low gradient colors for the coefficient magnitude. Will be fed into scale_fill_gradient
#' @param high_grad The high gradient colors for the coefficient magnitude. Will be fed into scale_fill_gradient
#' @param mid_grad The medium gradient colors for the coefficient magnitude. Will be fed into scale_fill_gradient
#' @param midpoint Midpoint for coefficient magnitude in legend
#' @param top_taxa Only plot X taxa with the largest magnitude beta coefficients
#' @param low_lim Lower limits of the fill gradient. Will default to the largest magnitude effect size
#' @param high_lim Upper limits of the fill gradient. Will default to the largest magnitude effect size
#' @param mute_cols Mute the colors of the fill gradients
#' @param alpha Mark beta coefficient cells with p-values below this cutoff
#' @param dot_size size of marker in cells
#' @param dot_shape shape of marker in cells
#' @param main Plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param subtitle Plot label
#' @param xaxis Labels for the x-axis ticks
#' @param legend_title Title of figure legend
#' @param caption plot caption to be displayed at the bottom of plot
#' @details The output will give gray columns if there are missing values in the supplied continuous variable
#' @return Returns a ggplot that you can add geoms to if you'd like
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Creating negative binomial models on filtered tidy_micro set
#' nb_fam <- set %>%
#' mutate(bpd1 = factor(bpd1)) %>% ## making bpd1 a factor
#' otu_filter(ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' nb_mods(table = "Family", bpd1)
#'
#' nb_fam %>% micro_heatmap
#' @export
micro_heatmap <- function(modsum, low_grad, high_grad, mid_grad, midpoint = 0, top_taxa = 10,
                       low_lim, high_lim, mute_cols = T, alpha = 0.05, dot_size = 2,
                       dot_shape = 8, main = NULL, xlab = NULL, ylab = NULL, subtitle = NULL,
                       xaxis = NULL, legend_title = NULL, caption = NULL){

  if(modsum$Model_Type %nin% c('bb_mod', 'nb_mod')){
    stop("'modsum' must be output from either nb_mods or bb_mods")
  }

  if(is.null(xaxis)) xaxis <- modsum$Model_Coef %>%
      dplyr::filter(!grepl("(Intercept)",.data$Coef)) %>%
      dplyr::distinct(.data$Coef) %>% purrr::simplify()

  if(missing(low_grad)) low_grad <- "blue"; if(missing(high_grad)) high_grad <- "red"
  if(missing(mid_grad)) mid_grad <- "beige"
  if(is.null(legend_title)) legend_title <- latex2exp::TeX("$\\hat{\\beta}$")
  if(is.null(caption)) caption <- latex2exp::TeX("'*' denotes significance at $\\alpha$ = 0.05")

  ## selecting the Taxa, model Coefs, and RR, and arranging them by RR within Taxa
  CC <- modsum$Convergent_Summary %>%
    dplyr::select(.data$Taxa, .data$Coef, .data$Beta, .data$FDR_Pval) %>%
    dplyr::filter(!grepl("(Intercept)", .data$Coef))

  if(missing(low_lim)) low_lim <- -max(abs(max(CC$Beta)),
                                       abs(min(CC$Beta)))
  if(missing(high_lim)) high_lim <- max(abs(max(CC$Beta)),
                                        abs(min(CC$Beta)))

  ## Organizing by strongest and most consistent effect size
  ar <- CC %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarise(max = max(.data$Beta), avg = mean(.data$Beta)) %>%
    dplyr::left_join(CC, by = "Taxa") %>%
    unique %>%
    dplyr::arrange(dplyr::desc(.data$max), dplyr::desc(.data$avg))

  ## Pulling out top taxa by effect size and plotting
  CC %>% dplyr::filter(.data$Taxa %in% unique(ar$Taxa)[seq(1,top_taxa)]) %>%

    ggplot2::ggplot(ggplot2::aes(.data$Coef, .data$Taxa)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Beta)) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient2(limits = c(low_lim, high_lim),
                         low = low_grad,
                         high = high_grad,
                         mid = mid_grad,
                         midpoint = midpoint) +
    ggplot2::labs(title = main, x = xlab, y = ylab,
                  subtitle = subtitle, caption = caption,
                  fill = legend_title) +
    ggplot2::scale_x_discrete(labels = xaxis) +
    ggplot2::geom_point(ggplot2::aes(size=ifelse(.data$FDR_Pval < alpha, "dot", "no_dot")), shape = 8) +
    ggplot2::scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
}

