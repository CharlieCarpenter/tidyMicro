#' @title Create 3d PCA plots
#' @name pca_3d
#' @description Create three dimensional PCA plots from longitudinal data or multiple omics data sets.
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param time_var The time point variable column name in your tidi_MIBI set
#' @param subject The subject variable column name in your tidi_MIBI set
#' @param y Value to calculate principle components or coordinates on. Default is centered log ratio (recommended)
#' @param modes Components of the data to focus on: time, subjects, bacteria, etc. "AC" by default
#' @param plot_scores Plot the scores instead of the principle components
#' @param dist_method Dissimilartiy method to be calculated by \code{\link[vegan]{vegdist}}. Euclidean by default
#' @param type "PCA" for principle components or "PCoA" to calculated dissimilarity matrix using \code{\link[vegan]{vegdist}}
#' @param plot_scores Plot the scores instead of the principle components
#' @param n_compA The number of components along first axis. See details
#' @param n_compB The number of components along second axis. See details
#' @param n_compC The number of components along third axis. See details
#' @param pch Plotting "character", i.e. symbol to use.
#' @param cex.axis Options for \code{\link[scatterplot3d]{scatterplot3d}}
#' @param cex.lab Options for \code{\link[scatterplot3d]{scatterplot3d}}
#' @param main Plot title
#' @param subtitle Plot subtitle
#' @param scalewt Logical; center and scale OTU table, recommended
#' @param print.legend Logical; print plot legend
#' @param legend.title Title for plot legend. Ignored if \code{print.legend = FALSE}
#' @param legend.position 'x' argument in \code{\link[graphics]{legend}}
#' @details Requires that you have columns for subject name and time point. Data must be complete across time points. The function will filter out inconsistent subjects
#'
#' When type = "PCoA" the component matrices must be specified prior to the optimization. This is handled automatically.
#'
#' If n_compA, n_compB, and n_compC aren't specified they will default to the number of complete subjects, the number of taxa, and the number of time points, respectively. This slows down performance slightly, but will not change the results.
#' @references \code{\link[vegan]{vegdist}}
#' @author Charlie Carpenter, Kayla Williamson
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#'
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin)
#'
#' set %>% pca_3d(table = "Family", time_var = day, subject = study_id)
#' @export
pca_3d <- function(micro_set, table, time_var, subject, y = clr,
                   modes = c("AC","BA","CB"), dist_method = "euclidean",
                   type = "PCoA", plot_scores = FALSE,
                   n_compA, n_compB, n_compC, pch = 16, cex.axis = 1, cex.lab = 1,
                   main = NULL, subtitle = NULL, scalewt = TRUE,
                   print.legend = FALSE, legend.title = NULL, legend.position = NULL){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")
  if(missing(modes)) modes <- "AC"; if(is.null(legend.position)) legend.position <- "right"

  ## Making the wide otu format that is needed for U_matrices
  wide_otu <- micro_set %>%
    A_matricization(table = table, subject = !!rlang::enquo(subject),
                    y = !!rlang::enquo(y),
                    time_var = !!rlang::enquo(time_var))

  ## Number of timepoints, subjects, and taxa
  n_time <- micro_set %>% dplyr::pull(!!rlang::enquo(time_var)) %>% unique %>% length
  n_sub <- nrow(wide_otu); n_taxa <- ncol(wide_otu)/n_time

  message("Found ", n_time, " time points and ", n_sub, " subjects with complete cases.\n")
  ## Pulls the maximum number of components if not specified
  ## Get better plots when you don't do this
  if(missing(n_compA)) n_compA <- n_sub
  if(missing(n_compB)) n_compB <- n_taxa
  if(missing(n_compC)) n_compC <- n_time
  ## Stops for incorrect dimensions
  if(n_compA > n_sub){
    stop("n_compA can't be greater than the number of subjects in every time point.")
  }
  if(n_compB > n_taxa){
    stop("n_compB can't be greater than the number of unique taxa in the specified table.")
  }
  if(n_compC > n_time){
    stop("n_compC can't be greater than the number of time points in your data set.")
  }
  if(modes == "CB" & n_time < 3){
    stop("Need at least 3 time points for CB modes.")
  }

  if(type %nin% c("PCA", "PCoA")) stop("'type' must be either 'PCA' or 'PCoA'")

  ## Scales the OTU counts which is recommended by KW Thesis
  if(scalewt) wide_otu %<>% ade4::scalewt()

  if(modes %nin% c("BA", "AC", "CB")) {
    stop("'modes' must be one of 'BA', 'AC', or 'CB'. The function will plot based on the non-specified components")
  }

  modes <- paste0(modes, "plot")

  if(type == "PCA"){
    ## Principle Components (start = 0)
    T3P <- ThreeWay::T3func(wide_otu, n = n_sub, m = n_taxa, p = n_time,
                            r1 = n_compA, r2 = n_compB, r3 = n_compC,
                            start = 0, conv = 0.0000000000000002)
    T3P <- T3_plots(r1 = n_compA, r2 = n_compB, r3 = n_compC,
                    T3P, modes, plot_scores = F)[modes] %>%
      as.data.frame
  } else if(type == "PCoA"){

    T3 <- U_matrices(wide_otu, method = dist_method, n_compA, n_compB, n_compC,
                     n_sub, n_taxa, n_time)

    ## Principle Coordinant (we define the component matrices, start = 2)
    T3P <- ThreeWay::T3func(wide_otu, n = n_sub, m = n_taxa, p = n_time,
                            r1 = n_compA, r2 = n_compB, r3 = n_compC,
                            start = 2, conv = 0.0000000000000002,
                            T3$A,T3$B,T3$C,T3$G)
    T3P <- T3_plots(r1 = n_compA, r2 = n_compB, r3 = n_compC,
                    T3P, modes, plot_scores)[modes] %>%
      as.data.frame
  }

  if(modes == "ACplot"){
    cc <- factor(rep(seq(1,n_time), each = n_sub))
    scatterplot3d::scatterplot3d(T3P[,1], T3P[,2], T3P[,3],
                                 color = cc, pch=pch,
                                 xlab = "B1", ylab = "B2", zlab = "B3",
                                 cex.axis = cex.axis, cex.lab = cex.lab,
                                 main = main, sub = subtitle)
    
    if(print.legend){
      ll <- micro_set %>% dplyr::pull(!!rlang::enquo(time_var)) %>% unique %>% sort
      graphics::legend(legend.position, legend = ll, 
             title = legend.title, pch = pch, col = levels(cc))
    } 
  }

  if(modes == "BAplot"){
    cc <- factor(rep(seq(1,n_taxa), each = n_sub))
    scatterplot3d::scatterplot3d(T3P[,1], T3P[,2], T3P[,3],
                                 color = cc, pch=pch, 
                                 xlab = "C1", ylab = "C2", zlab = "C3",
                                 cex.axis = cex.axis, cex.lab = cex.lab,
                                 main = main, sub = subtitle)
  
    if(print.legend){
      ll <- micro_set %>% 
        dplyr::filter(Table == table) %>% 
        dplyr::pull(Taxa) %>% unique
      
      graphics::legend(legend.position, legend = ll, 
             title = legend.title, pch = pch, col = levels(cc))
    }
  }

  if(modes == "CBplot"){
    cc <- factor(rep(seq(1, n_taxa), each = n_time))
    scatterplot3d::scatterplot3d(T3P[,1], T3P[,2], T3P[,3],
                                 color = cc, pch=pch, 
                                 xlab = "A1", ylab = "A2",zlab = "A3",
                                 cex.axis = cex.axis, cex.lab = cex.lab,
                                 main = main, sub = subtitle)
    
    if(print.legend){
      ll <- micro_set %>% dplyr::pull(!!rlang::enquo(subject)) %>% unique
      graphics::legend(legend.position, legend = ll, 
             title = legend.title, pch = pch, col = levels(cc))
    } 
  }
}

