#' @title Create Rocky Mountain plots from taxa relative abundance correlations
#' @name cor_rocky_mtn
#' @description Calculate the correlation between the relative abundance of each taxa within a specified table and a continuous variable of interest. Correlation is calculated by \code{\link[stats]{cor}}. By default, Kendall's correlation is used to account for the prevalence of ties that often occur (lots of 0s)
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param x Continuous variable of interest
#' @param y The taxa information. The centered log ratio (clr) is recommended.
#' @param method Correlation type; must be supported by \code{\link[stats]{cor}}. By default it is "spearman" to use with clr. If you'd like to use taxa ra, it is recommend you switch to Kendall's correlation to account for the large number of ties common in taxa ra (lots of 0s)
#' @param main Plot title
#' @param xlab Lable for x-axis
#' @param ylab Label for y-axis
#' @param subtitle Plot subtitle
#' @param cut_lines Add lines for p-value cutoffs
#' @param line_text Label p-value cut-offs
#' @param sig_text Label taxa with correlations greater than \code{cor_label} in magnitude
#' @param lwd line width for cut_lines
#' @param cor_label Cutoff for correlations to be labeled
#' @param breaks Where to place cut_lines along y-axis
#'
#' @author Charlie Carpenter, Dan Frank
#' @return A ggplot you can add geoms to if you'd like
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' set %>% cor_rocky_mtn(table = "Family", weight, cor_label = 0.3)
#' @export
cor_rocky_mtn <- function(micro_set, table, x, y = clr,
                          method = c("pearson", "kendall", "spearman"),
                          main = NULL, xlab = NULL, ylab = NULL, subtitle = NULL,
                          cut_lines = TRUE, line_text = TRUE,
                          sig_text = TRUE, lwd = 1, cor_label = 0.5,
                          breaks = c(-0.6, -0.5, -0.3, 0.3, 0.5, 0.6)){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  if(missing(method)) method <- "spearman"
  if(method %nin% c("pearson", "kendall", "spearman")){
    stop("'method' must be one of: pearson, kendall, spearman")
  }

  if(is.null(ylab)) ylab = paste(method, "correlations") %>% stringr::str_to_title()

  cor_set <- micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarise(CC = stats::cor(!!rlang::enquo(x), !!rlang::enquo(y),
                                     method = method, use = "pairwise.complete.obs"))

  if(table != "Phylum"){
    cor_set %<>%
      dplyr::filter(!is.na(.data$CC)) %>%
      dplyr::mutate(phyl = sapply(.data$Taxa, phy_fun)) %>%
      phyl_ord() %>% ## Make a Phylum variable for color
      dplyr::arrange(.data$phyl) %>%
      dplyr::mutate(ys = 0) %>% ## start of y segment
      tibble::rowid_to_column(var = "pos") ## Position along x-axis

    gg <- cor_set %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$pos, y = .data$CC)) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$pos, xend =.data$pos,
                                         y = .data$ys, yend = .data$CC,
                                         colour = .data$phyl),
                            data = cor_set)
  }else{
    cor_set %<>%
      dplyr::filter(!is.na(.data$CC)) %>%
      Taxa_ord() %>%
      dplyr::arrange(.data$Taxa) %>%
      dplyr::mutate(ys = 0 ) %>%  ## start of y segment
      tibble::rowid_to_column(var = "pos") ## Position along x-axis

    gg <- cor_set %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$pos, y = .data$CC)) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$pos, xend = .data$pos,
                                         y = .data$ys, yend = .data$CC,
                                         colour = .data$Taxa),
                            data = cor_set)
  }

  if(cut_lines){
    gg <- gg + ggplot2::geom_hline(yintercept = breaks, linetype = "dotdash")
  }

  if(line_text){
    gg <- gg +
      ggplot2::scale_y_continuous(breaks = breaks,
                                  labels = as.character(breaks))
  }else{
    gg <- gg +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }

  if(sig_text){
    gg <- gg + ggrepel::geom_label_repel(data = cor_set %>%
                                           dplyr::filter(abs(.data$CC) >= cor_label),
                                         ggplot2::aes(label = .data$Taxa),
                                         show.legend = F)
  }

  gg +
    ggplot2::theme_bw() +
    ggplot2::labs(colour = "Phylum", y = ylab, x = xlab, title = main, subtitle = subtitle) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0)
}
