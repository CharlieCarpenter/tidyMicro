#' @title Create heatmaps of the supplied dissimilarity matrices
#' @name beta_heatmap
#' @description Create heatmaps of the supplied dissimilarity matrices ordered by supplied grouping variables
#' @param beta_div A dissimilarity matrix calculated by \code{beta_div}
#' @param micro_set A tidy_micro data set
#' @param ... Variables for ordering
#' @param low_grad Colors for the corelation magnitude. Will be fed into scale_fill_gradient
#' @param high_grad Colors for the corelation magnitude. Will be fed into scale_fill_gradient
#' @param main Plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param subtitle Plot label
#' @param natural_order Keep order of axes in the conventional order for dissimilarity matrices
#' @param legend_title Title for the legend
#' @return Returns a ggplot that you can add geoms to if you'd like
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
#' bray %>% beta_heatmap(micro_set = set, bpd1)
#' @export
beta_heatmap <- function(beta_div, micro_set, ..., low_grad, high_grad,
                         main = NULL, xlab = NULL, ylab = NULL, subtitle = NULL,
                         natural_order = TRUE,
                         legend_title = "Dissimilarity"){

  if(!missing(low_grad) & missing(high_grad)) stop("Must specify both low_grad and high_grad.")
  if(missing(low_grad) & !missing(high_grad)) stop("Must specify both low_grad and high_grad.")

  micro_set %<>%
    dplyr::distinct(.data$Lib, .keep_all = TRUE) %>% ## unique subjects from micro_set
    dplyr::select(.data$Lib, ...)

  if(ncol(micro_set) > 2) stop("Must use one factor variable for `...` argument.")
  if(class(micro_set[,2]) %nin% c("character", "factor")) stop("Must use one factor variable for `...` argument.")

  ## Ensuring colnames are the Lib names for gather step
  colnames(beta_div) <- micro_set$Lib
  ## For reording leves
  Var <- micro_set[, 2]

  CC <- beta_div %>%
    as.data.frame %>% ## Getting beta diversity file into format for ggplot
    tibble::rownames_to_column(var = "Lib") %>%
    dplyr::full_join(micro_set, by = "Lib") %>%
    tidyr::unite("II", ..., .data$Lib, sep = ":_:") %>%
    ## Supposed to be a 'unique' sep so we can locate it later
    ## Will break if ':_:' is in someones Libs or variable name
    tidyr::pivot_longer(-.data$II, names_to = "Lib2", values_to = "val") %>%
    dplyr::full_join(micro_set %>% dplyr::rename(Lib2 = .data$Lib), by = "Lib2") %>%
    tidyr::unite("III", ..., .data$Lib2, sep = ":_:")

  if(natural_order){
    CC %<>%
      dplyr::mutate(II = factor(.data$II),
                    II = factor(.data$II,
                                levels = rev(levels(.data$II)))
      ) ## gets the same subject beta divs on the main diagonal
    Var_y <- factor(Var, levels = rev(levels(Var)))

    gg <- ggplot2::ggplot(CC, ggplot2::aes(.data$II, .data$III, fill = .data$val)) +
      ggplot2::geom_tile() +
      ggplot2::labs(x = xlab, y = ylab, subtitle = subtitle, title = main, fill = legend_title) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)) +
      ggplot2::scale_y_discrete(labels = rep(levels(Var_y), table(Var_y))) +
      ggplot2::scale_x_discrete(labels = rep(levels(Var), table(Var)))
  } else {
    gg <- ggplot2::ggplot(CC, ggplot2::aes(.data$II, .data$III, fill = .data$val)) +
      ggplot2::geom_tile() +
      ggplot2::labs(x = xlab, y = ylab, subtitle = subtitle, title = main, fill = legend_title) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)) +
      ggplot2::scale_y_discrete(labels = rep(levels(Var), table(Var))) +
      ggplot2::scale_x_discrete(labels = rep(levels(Var), table(Var)))
  }
  if(!missing(low_grad) & !missing(high_grad)){
    gg <- gg + ggplot2::scale_fill_gradient(low = low_grad, high = high_grad)
  }
  gg
}

