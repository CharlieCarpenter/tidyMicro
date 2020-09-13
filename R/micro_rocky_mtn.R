#' @title Create Rocky Mountain plots from negative binomial taxa models
#' @name micro_rocky_mtn
#' @description Display the magnitude of log p-values for each of the taxa in \code{nb_mods} as vertical bars next to each other along the x-axis. The direction of the bars will be determined by the direction of the estimated relationship. The taxa will be color coded by the phylum they belong to, and taxa that have FRD adjusted p-values below your desired significance cutoff for the specified covariate will be labeled
#' @param modsum The output from nb_mods
#' @param ... The covariate you'd like to plot. Must be in the models created by nb_mods
#' @param main Plot title
#' @param ylab y-axis labels
#' @param subtitle Plot subtitle
#' @param pval_lines Logical; include horizonal dashed lines at corresponding p-values
#' @param pval_text Logical; label the y-axis with corresponding p-values
#' @param sig_text Logical; label the taxa with p-values below specified alpha
#' @param facet_labels Labels for different facets if covariate has more than 1 beta coefficient
#' @param alpha Significance cutoff
#' @param lwd Line width for pval_lines
#' @param lty Line type for pval_lines
#' @return A ggplot you can add geoms to if you'd like
#' @author Charlie Carpenter, Rachel Johnson, Dan Frank
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#'
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Creating negative binomial models on filtered tidy_micro set
#' nb_fam <- set %>%
#' otu_filter(ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' nb_mods(table = "Family", bpd1)
#'
#' nb_fam %>% micro_rocky_mtn(bpd1)
#' @export
micro_rocky_mtn <- function(modsum, ..., main = NULL,
                         ylab = NULL, subtitle = NULL,
                         pval_lines = TRUE, pval_text = TRUE, sig_text = TRUE,
                         facet_labels = NULL, alpha = 0.05, lwd = 2, lty = 1){

  if('Model_Type' %nin% names(modsum)){
    stop("'modsum' must be output from either nb_mods or bb_mods")
  }
  if(modsum$Model_Type %nin% c('nb_mod', 'bb_mod') ){
    stop("'modsum' must be output from either nb_mods or bb_mods")
  }

  if(missing(...)) stop("NB_RockMtn requires a model coefficient")
  if(is.null(ylab)) ylab <- expression("log"[10]*" p-value")

  CC <- modsum$Convergent_Summary %>%
    dplyr::filter(.data$Coef != "(Intercept)", stringr::str_detect(.data$Coef, cov_str(...))) %>%
    dplyr::mutate(FDR_Pval = ifelse(.data$FDR_Pval == 0.0000, 0.0005, .data$FDR_Pval),
           phyl = sapply(.data$Taxa, phy_fun)) %>%
    phyl_ord()

  CC %<>%
    dplyr::mutate(Coef = factor(.data$Coef, levels = unique(.data$Coef)),
           rnum = seq(1,nrow(CC)),
           ys = rep(0,nrow(CC)),
           ye = ifelse(.data$Beta < 0, log(.data$FDR_Pval, base = 10), -log(.data$FDR_Pval, base = 10)))

  pline <- c(-1,1)*sapply(rep(c(0.01, 0.05, 0.1), 2), log, base = 10)

  gg <- ggplot2::ggplot(CC, ggplot2::aes(x = .data$rnum, y = .data$ye))

  if(pval_lines){
    gg <- gg + ggplot2::geom_hline(yintercept = pline, linetype = "dotdash")
  }

  if(pval_text){
    gg <- gg +
      ggplot2::scale_y_continuous(breaks = pline,
                                  labels = c("0.01","0.05","0.10", "0.01","0.05","0.10"))
  } else{
    gg <- gg +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }

  gg <- gg +
    ggplot2::geom_segment(ggplot2::aes(x = .data$rnum, xend = .data$rnum,
                                      y = .data$ys, yend = .data$ye, colour = .data$phyl),
                          size = lwd) +
    ggplot2::theme_bw() +
    ggplot2::labs(colour = "Phylum", y = ylab, title = main, subtitle = subtitle,
                  x = NULL) +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank())

  if(length(unique(CC$Coef)) > 1L){
    if(!is.null(facet_labels) & is.null(names(facet_labels))){
      names(facet_labels) <- unique(CC$Coef)
    }

    gg <- gg + ggplot2::facet_wrap( ~ Coef, labeller = ggplot2::labeller(Coef = facet_labels))
  }

  if(sig_text){
   gg <- gg + ggrepel::geom_label_repel(data = CC %>% dplyr::filter(.data$FDR_Pval < alpha),
     ggplot2::aes(label = .data$Taxa), show.legend = F)
  }

  gg +
    ggplot2::geom_hline(yintercept = 0)
}
