#' @title Calculate and plot principle components
#' @name micro_pca
#' @description Principle components are calculated on the centerted log ratio tranformation of the OTU table using the \code{\link[stats]{prcomp}} function from the \code{\link{stats}} package. Scaling the OTU table to a unit variance is the default option, and recommended, but this can be changed using scaled = F.
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param dist A distance matrix, such as a beta diversity. If supplied a PCoA plot will be returned
#' @param grp_var Categorical grouping variable
#' @param y Value to calculate principle components or coordinates on. Default is centered log ratio (recommended)
#' @param scale Logical. Indicating whether the variables should be scaled to have unit variance before the analysis takes place
#' @param axes_arrows Logical. Plot component axes arrows
#' @param ellipse Logical. Plot normal data ellipses by groups
#' @param ellipse.prob Numeric.
#' @param main Plot title
#' @param subtitle Plot subtitle
#' @param legend_title Legend title
#' @details PCA calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix, not by using eigen on the covariance matrix. This is generally the preferred method for numerical accuracy. Calculations are accomplished through the \code{\link[stats]{prcomp}} function, and the plot is created through internal code based on the ggbiplot function \url{https://github.com/vqv/ggbiplot}.
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.
#'
#' Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) Multivariate Analysis, London: Academic Press.
#'
#' Venables, W. N. and B. D. Ripley (2002) Modern Applied Statistics with S, Springer-Verlag.
#'
#' Vincent Q. Vu (2011). ggbiplot: A ggplot2 based biplot. \url{https://github.com/vqv/ggbiplot}
#' @return A ggplot you can add geoms to if you'd like
#' @examples
#' data(mrsa_gen); data(mrsa_clin)
#'
#' set <- tidy_micro(otu_tabs = mrsa_gen, tab_names = "Genus", clinical = mrsa_clin)
#'
#' ## PCA Plot
#' set %>% micro_pca(table = "Genus", grp_var = Aureus_Positive)
#'
#' ## PCoA Plot (Recommended for p > n)
#'
#' bray_beta <- set %>% beta_div(table = "Genus")
#' micro_pca(set, dist = bray_beta, grp_var = Aureus_Positive, ellipse = TRUE)
#'
#' @export
micro_pca <- function(micro_set, table = NULL, dist = NULL, grp_var, y = clr,
                      scale = TRUE, axes_arrows = F,
                      ellipse = FALSE, ellipse.prob = 0.68,
                      main = NULL, subtitle = NULL, legend_title = NULL){

  if(missing(grp_var)) stop("'grp_var' must be specified")

  if(is.null(dist)){
    ## PCA using
    if(is.null(table)) stop("Please supply a table with your tidy_micro set")

    if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

    dd <- micro_set %>%
      dplyr::filter(.data$Table == table) %>% ## Filtering out table we want
      dplyr::select(.data$Lib, .data$Taxa, !!rlang::enquo(y), !!rlang::enquo(grp_var)) %>% ## getting impartant vars
      tidyr::spread(.data$Taxa, !!rlang::enquo(y)) ## Formatting for PCA

    ## % of variation explained by first two components
    grp_var <- dd %>% dplyr::pull(!!rlang::enquo(grp_var))

    dd %<>%
      dplyr::select(-c(1,2)) %>%  ## Removing Lib and var
      stats::prcomp(scale = scale)

    ve <- round(summary(dd)$importance[2,c(1,2)], 3) * 100

    gg <- dd %>%   ## Principle components
      micro_biplot(var.axes = axes_arrows, ## Annoying arrows
                   groups = grp_var, ## pulling var for grouping
                   ellipse = ellipse, ellipse.prob = ellipse.prob
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, subtitle = subtitle,
                    x = paste0("Principle Component 1 (", ve[1], "%)"),
                    y = paste0("Principle Component 2 (", ve[2], "%)"), colour = legend_title)

    message("PCA plot created")
  } else {
    if(!is.matrix(dist)) stop("'dist' must be a distance matrix. Usually created from 'beta_div' or a similar function")

    ## Not the most elegant but it works...
    ## pulling grouping var for color
    grp_var <- micro_set %>%
      dplyr::filter(.data$Lib %in% rownames(dist)) %>% ## Only including Libs from dist
      dplyr::distinct(.data$Lib, .keep_all = TRUE) %>% ## removing repeated libs
      dplyr::pull(!!rlang::enquo(grp_var)) ## Pulling grouping var

    dd <- dist %>%
      stats::prcomp(scale = scale) ## Principle coordinants

    ve <- round(summary(dd)$importance[2,c(1,2)], 3) * 100

    gg <- dd %>%
      micro_biplot(var.axes = axes_arrows, ## Annoying arrows
                   groups = grp_var, ## pulling var for grouping
                   ellipse = ellipse, ellipse.prob = ellipse.prob
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, subtitle = subtitle,
                    x = paste0("Principle Component 1 (", ve[1], "%)"),
                    y = paste0("Principle Component 2 (", ve[2], "%)"), colour = legend_title)

    message("PCoA plot created")
  }

  gg
}
