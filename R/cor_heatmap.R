#' @title Create correlation heatmaps of taxa and another continuous variable
#' @name cor_heatmap
#' @description Calculated the correlation between a specified continuous variable and some taxa measure. Correlation type and taxa measure (count, relative abundance, etc.) can be specified by the user but is "spearman" and relative abundance, respectively, by default
#' @param micro_set A tidy_micro data set
#' @param table The OTU table
#' @param ... Continuous variables of interest
#' @param y The taxa information: cts, ra, etc. The centered log ratio (clr) is recommended.
#' @param low_grad Colors for the corelation magnitude. Will be fed into scale_fill_gradient
#' @param high_grad Colors for the corelation magnitude. Will be fed into scale_fill_gradient
#' @param method Correlation type; must be supported by \code{\link[stats]{cor}}. By default it is "spearman" to use with clr. If you'd like to use taxa ra, it is recommend you switch to Kendall's correlation to account for the large number of ties common in taxa ra (lots of 0s)
#' @param main Plot title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param subtitle Plot label
#' @param legend_title Title for the legend
#' @details The output will give gray columns if there are missing values in the supplied continuous variable
#' @return Returns a ggplot that you can add geoms to if you'd like
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' set %>% cor_heatmap(table = "Class", gestational_age, weight)
#' @export
cor_heatmap <- function(micro_set, table, ..., y = clr,
                        method = c("pearson", "kendall", "spearman"),
                        main = NULL, xlab = NULL, ylab = NULL, subtitle = NULL,
                        legend_title = NULL, low_grad, high_grad){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  if(missing(method)) method <- "spearman"
  if(method %nin% c("pearson", "kendall", "spearman")){
    stop("method must be one of: pearson, kendall, spearman")
  }

  if(is.null(legend_title)){
    legend_title <- paste0("spearman \n", "correlation") %>% stringr::str_to_title()
  }

  ## Calculates correlation of all '...' put in, and gathers it for headmap format
  gg <-  micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::select(.data$Taxa, !!rlang::enquo(y), ...) %>%
    ## Correlations for each selected var in ...
    plyr::ddply(~ .data$Taxa, .fun = function(set, taxa_info, .method, ...){
      .y <- set %>% dplyr::pull(!!rlang::enquo(taxa_info))
      set %>% dplyr::select(...) %>% 
      apply(2, stats::cor, y = .y,
            use = "pairwise.complete.obs", method = .method)
    }, taxa_info = !!rlang::enquo(y), .method = method, ...) %>%
    dplyr::rename(Taxa = .data$`.data$Taxa`) %>%
    ## making long format
    tidyr::pivot_longer(-.data$Taxa, names_to = "Var", values_to = "corr") %>%

    ## Making pase ggplot
    ggplot2::ggplot(ggplot2::aes(.data$Var, .data$Taxa)) +
    ggplot2::geom_tile(ggplot2::aes(fill = as.numeric(.data$corr))) +
    ggplot2::theme_bw()

  if(!missing(low_grad) & missing(high_grad)) stop("Must specify both low_grad and high_grad.")
  if(missing(low_grad) & !missing(high_grad)) stop("Must specify both low_grad and high_grad.")

  if(!missing(low_grad) & !missing(high_grad)){
    gg <- gg + ggplot2::scale_fill_gradient(low = low_grad, high = high_grad)
  }

  gg +
    ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = legend_title))
}

