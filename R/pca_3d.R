#' @title Create 3d PCA plots
#' @name pca_3d
#' @description Create three dimensional PCA plots from longitudinal data or multiple omics data sets.
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param time_var The time point variable column name in your tidi_MIBI set
#' @param subject The subject variable column name in your tidi_MIBI set
#' @param y Value to calculate principle components or coordinates on. Default is centered log ratio (recommended)
#' @param plot_scores Plot the scores instead of the principle components
#' @param dist_method Dissimilartiy method to be calculated by \code{\link[vegan]{vegdist}}. Euclidean by default
#' @param type "PCA" for principle components or "PCoA" to calculated dissimilarity matrix using \code{\link[vegan]{vegdist}}
#' @param plot_scores Plot the scores instead of the principle components
#' @param pch Plotting "character", i.e. symbol to use.
#' @param cex.axis Options for \code{\link[scatterplot3d]{scatterplot3d}}
#' @param cex.lab Options for \code{\link[scatterplot3d]{scatterplot3d}}
#' @param cex Options for \code{\link[scatterplot3d]{scatterplot3d}}
#' @param main Plot title
#' @param subtitle Plot subtitle
#' @param scalewt Logical; center and scale OTU table, recommended
#' @param print.legend Logical; print plot legend
#' @param legend.title Title for plot legend. Ignored if \code{print.legend = FALSE}
#' @param legend.position 'x' argument in \code{\link[graphics]{legend}}
#' @details Requires that you have separate columns for subject ID and time point. Data must be complete across time points. The function will automatically filter out incomplete cases with a warning message.
#'
#' When type = "PCoA" the component matrices must be specified prior to the optimization. This is handled automatically.
#'
#' @references \code{\link[vegan]{vegdist}}
#' @author Charlie Carpenter, Kayla Williamson
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#'
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin)
#'
#' set %>% pca_3d(table = "Family", time_var = day, subject = study_id, legend.title = "Day")
#' @importFrom scatterplot3d scatterplot3d
#' @export
pca_3d <- function(micro_set, table, time_var, subject, y = clr,
                   dist_method = "euclidean", type = "PCoA", plot_scores = FALSE,
                   pch = 16, cex.axis = 1, cex.lab = 1, cex = 1,
                   main = NULL, subtitle = NULL, scalewt = TRUE,
                   print.legend = TRUE,
                   legend.title = "Time Points", legend.position = 'right'){

#  modes = c("AC","BA","CB"),
#  #' @param modes Components of the data to focus on: time, subjects, bacteria, etc. "AC" by default

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")
  modes <- "ACplot"

  if(missing(subject)){
    stop("Must supply a column name for subject ID")
  }
  if(cov_str(!!rlang::enquo(subject)) %nin% names(micro_set)){
    stop("subject is not a column name in supplied micro_set")
  }

  if(missing(time_var)){
    stop("Must supply a column name for the time variable")
  }
  if(cov_str(!!rlang::enquo(time_var)) %nin% names(micro_set)){
    stop("time_var is not a column name in supplied micro_set")
  }

  ## Making the wide otu format that is needed for U_matrices
  wide_otu <- micro_set %>%
    A_matricization(table = table, subject = !!rlang::enquo(subject),
                    y = !!rlang::enquo(y),
                    time_var = !!rlang::enquo(time_var))

  ## Number of timepoints, subjects, and taxa
  time_v <- micro_set %>% dplyr::pull(!!rlang::enquo(time_var)) %>% unique %>% sort
  n_time <- length(time_v)
  n_sub <- nrow(wide_otu); n_taxa <- ncol(wide_otu)/n_time

  message("Found ", n_time, " time points and ", n_sub, " subjects with complete cases.")
  ## Pulls the maximum number of components if not specified
  ## Get better plots when you don't do this
  n_compA <- n_sub; n_compB <- n_taxa; n_compC <- n_time

  if(type %nin% c("PCA", "PCoA")) stop("'type' must be either 'PCA' or 'PCoA'")

  ## Scales the OTU counts which is recommended by KW Thesis
  if(scalewt) wide_otu %<>% ade4::scalewt()

  if(type == "PCA"){
    ## Principle Components (start = 0)
    T3P <- ThreeWay::T3func(wide_otu, n = n_sub, m = n_taxa, p = n_time,
                            r1 = n_compA, r2 = n_compB, r3 = n_compC,
                            start = 0, conv = 0.0000000000000002)
    T3P <- T3_plots(r1 = n_compA, r2 = n_compB, r3 = n_compC,
                    T3P, modes, plot_scores)

    pv <- round(100* T3P$B_eig[1:3]/sum(T3P$B_eig), 2)
    T3P <- T3P[modes] %>% as.data.frame
  } else if(type == "PCoA"){

    T3 <- U_matrices(wide_otu, method = dist_method, n_compA, n_compB, n_compC,
                     n_sub, n_taxa, n_time)

    ## Principle Coordinant (we define the component matrices, start = 2)
    T3P <- ThreeWay::T3func(wide_otu, n = n_sub, m = n_taxa, p = n_time,
                            r1 = n_compA, r2 = n_compB, r3 = n_compC,
                            start = 2, conv = 0.0000000000000002,
                            T3$A,T3$B,T3$C,T3$G)
    T3P <- T3_plots(r1 = n_compA, r2 = n_compB, r3 = n_compC,
                    T3P, modes, plot_scores)

    pv <- round(100* T3P$B_eig[1:3]/sum(T3P$B_eig), 2)
    T3P <- T3P[modes] %>% as.data.frame
  }

  if(modes == "ACplot"){
    cc <- factor(rep(seq(1,n_time), each = n_sub))
    lb <- paste0("Three Mode Component ", 1:3, " (",  pv, ")%")
    scatterplot3d::scatterplot3d(T3P[,1], T3P[,2], T3P[,3],
                                 color = cc, pch=pch,
                                 xlab = lb[1], ylab = lb[2], zlab = lb[3],
                                 cex.axis = cex.axis, cex.lab = cex.lab,
                                 main = main, sub = subtitle)

    if(print.legend){
      graphics::legend(legend.position, legend = time_v, cex = cex,
             title = legend.title, pch = pch, col = levels(cc))
    }
  }
}


  # if(modes == "BAplot"){
  #   cc <- factor(rep(seq(1,n_taxa), each = n_sub))
  #   cols. <-
  #
  #   scatterplot3d::scatterplot3d(T3P[,1], T3P[,2], T3P[,3],
  #                                color = cc, pch=pch,
  #                                xlab = "C1", ylab = "C2", zlab = "C3",
  #                                cex.axis = cex.axis, cex.lab = cex.lab,
  #                                main = main, sub = subtitle)
  #
  #   if(print.legend){
  #     ll <- micro_set %>%
  #       dplyr::filter(.data$Table == table) %>%
  #       dplyr::pull(.data$Taxa) %>% unique %>% pull.lev(4)
  #
  #     plot.new()
  #     graphics::legend(legend.position, legend = ll, cex = cex,
  #            title = legend.title, pch = pch, col = topo.colors(length(levels(cc))))
  #   }
  # }
  #
  # if(modes == "CBplot"){
  #   cc <- factor(rep(seq(1, n_taxa), each = n_time))
  #   scatterplot3d::scatterplot3d(T3P[,1], T3P[,2], T3P[,3],
  #                                color = cc, pch=pch,
  #                                xlab = "A1", ylab = "A2",zlab = "A3",
  #                                cex.axis = cex.axis, cex.lab = cex.lab,
  #                                main = main, sub = subtitle)
  #
  #   if(print.legend){
  #     plot.new()
  #     ll <- micro_set %>% dplyr::pull(study_id) %>% unique
  #     graphics::legend(legend.position, legend = ll, cex = cex,
  #            title = legend.title, pch = pch, col = levels(cc), ...)
  #   }
  # }

# ## Stops for incorrect dimensions
# if(n_compA > n_sub){
#   stop("n_compA can't be greater than the number of subjects in every time point.")
# }
# if(n_compB > n_taxa){
#   stop("n_compB can't be greater than the number of unique taxa in the specified table.")
# }
# if(n_compC > n_time){
#   stop("n_compC can't be greater than the number of time points in your data set.")
# }
# if(modes == "CB" & n_time < 3){
#   stop("Need at least 3 time points for CB modes.")
# }

# if(modes %nin% c("BA", "AC", "CB")) {
#   stop("'modes' must be one of 'BA', 'AC', or 'CB'. The function will plot based on the non-specified components")
# }
