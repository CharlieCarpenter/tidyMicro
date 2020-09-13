#' @title Create Three Mode PCA and PCoA plots
#' @name three_mode
#' @description Three Mode Principal Components, an ordination method that can take into account repeated measure of subjects. These methods have also been extended to other common ecological distance metrics for Three Mode Principal Coordinate Analysis
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param group A categorical variable to color by
#' @param time_var The time point variable column name in your tidi_MIBI set
#' @param subject The subject variable column name in your tidi_MIBI set
#' @param y Value to calculate principle components or coordinates on. Default is centered log ratio (recommended)
#' @param plot_scores Plot the scores instead of the principle components
#' @param main Plot title
#' @param subtitle Plot subtitle
#' @param legend_title Plot legend title
#' @param scalewt Logical; center and scale OTU table, recommended
#' @details Requires that you have columns for subject name and time point. Data must be complete across time points. The function will filter out inconsistent subjects
#'
#' If n_compA, n_compB, and n_compC aren't specified they will default to the number of complete subjects, the number of taxa, and the number of time points, respectively. This slows down performance slightly, but will not change the results.
#' @return A ggplot you can add geoms to if you'd like
#' @author Charlie Carpenter, Kayla Williamson
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#'
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin)
#'
#' set %>% three_mode(table = "Family", group = bpd1, time_var = day, subject = study_id)
#' @export
three_mode <- function(micro_set, table, group, time_var, subject, y = clr,
                       plot_scores = F, main = NULL, subtitle = NULL,
                       legend_title = NULL, scalewt = TRUE){

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  if(missing(subject)) stop("subject variable missing with no default.")
  if(missing(time_var)) stop("time_var variable missing with no default.")
  modes <- "AC"

  wide_otu <- micro_set %>%
    A_matricization(table = table, subject = !!rlang::enquo(subject),
                  y = !!rlang::enquo(y),
                  time_var = !!rlang::enquo(time_var))

  ## Number of timepoints, subjects, and taxa
  n_time <- micro_set %>% dplyr::pull(!!rlang::enquo(time_var)) %>% unique %>% length
  n_sub <- nrow(wide_otu); n_taxa <- ncol(wide_otu)/n_time

  message("Found ", n_time, " time points and ", n_sub, " subjects with complete cases.\n")
  ## Pulls the maximum number of components
  n_compA <- nrow(wide_otu); n_compB <- ncol(wide_otu)/n_time; n_compC <- n_time

  ## Scales the OTU counts which is recommended by KW Thesis
  if(scalewt) wide_otu %<>% ade4::scalewt()

  ## Pulling out important var and arranging by Subject
  meta <- micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::select(sub = !!rlang::enquo(subject), !!rlang::enquo(group)) %>%
    ## Getting distinct rows
    dplyr::distinct() %>%
    ## Only including the subjects in every time point
    dplyr::filter(.data$sub %in% rownames(wide_otu)) %>%
    dplyr::arrange(.data$sub)

  modes <- paste0(modes, "plot")

  Tmod <- ThreeWay::T3func(wide_otu, n = n_sub, m = n_taxa,
                   p = n_time, r1 = n_compA, r2 = n_compB, r3 = n_compC,
                   start = 0,conv = 0.0000000000000002)
  Tmod.p <- T3_plots(n_compA, n_compB, n_compC,
                     Tmod, modes, plot_scores)

  ## x,y Labs
  p.var <- paste0(100*round(Tmod.p$B_eig[1:2]/sum(Tmod.p$B_eig), 4), "%")
  ll <- paste0("Principle Component ", c("B1","B2"), " (", p.var, ")")

  Tmod.p <- Tmod.p[modes] %>% as.data.frame
  gg <- suppressWarnings(data.frame(x = Tmod.p[,1], y = Tmod.p[,2],
                                    sub = rep(rownames(wide_otu), each = n_time)) %>%
                           dplyr::full_join(meta, by = "sub") %>%
                           ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y,
                                                        col = !!rlang::enquo(group))) +
                           ggplot2::geom_point() +
                           ggplot2::labs(title = main,
                                         subtitle = subtitle, x = ll[1], y = ll[2])
  )

  if(!is.null(legend_title)) gg <- gg + ggplot2::labs(col = legend_title)

  gg + ggplot2::theme_bw()
}


# if(modes == "BAplot"){
#   xlab = "C1"; ylab = "C2"
#
#   gg <- data.frame(x = Tmod.p[,1], y = Tmod.p[,2],
#                    sub = rep(rownames(wide_otu), each = n_sub*n_taxa)) %>%
#     dplyr::full_join(meta, by = "sub") %>%
#     ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y,
#                                  col = !!rlang::enquo(group))) +
#     ggplot2::geom_point() +
#     ggplot2:: labs(title = main, subtitle = subtitle, x = xlab, y = ylab)
# }
#
# if(modes == "CBplot"){
#   xlab = "A1"; ylab = "A2"
#
#   gg <- data.frame(x = Tmod.p[,1], y = Tmod.p[,2],
#                    sub = rep(rownames(wide_otu), each = n_taxa*n_time)) %>%
#     dplyr::full_join(meta, by = "sub") %>%
#     ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y,
#                                  col = !!rlang::enquo(group))) +
#     ggplot2::geom_point() +
#     ggplot2::labs(title = main, subtitle = subtitle, x = xlab, y = ylab)
# }
#

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

# if(modes %nin% c("BA", "AC", "CB")) {
#   stop("'modes' must be one of 'BA', 'AC', or 'CB'. The function will plot based on the non-specified components")
# }

