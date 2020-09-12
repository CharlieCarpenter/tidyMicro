#' @title Function to make boxplots of taxa counts or relative abundance
#' @name taxa_boxplot
#' @description A function to make boxplots of one specified taxa relative abundance with the option to stratify by a factor variable
#' @param micro_set A tidy_micro data set
#' @param taxa A character string. The name of the taxa of interest
#' @param ... The factor variable you'd like to stratify by
#' @param y The taxa information
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param main Plot title
#' @param subtitle Subtitle for the plot
#' @param legend_title Title of plot legend
#' @return A ggplot that you can add geoms to if you'd like
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' set %>%
#' taxa_boxplot("Firmicutes/Bacilli/Bacillales/Staphylococcaceae", bpd1)
#' @export
taxa_boxplot <- function(micro_set, taxa, ..., y = ra, xlab = NULL, ylab = NULL,
                         main = NULL, subtitle = NULL, legend_title = NULL){

  if(!is.character(taxa)) stop("taxa must be a character string")
  if(taxa %nin% micro_set$Taxa) stop(paste(taxa, "not in supplied tidy_micro set. \n"))
  if(is.null(ylab)) ylab <- paste(taxa,
                                  stringr::str_to_upper(cov_str(!!rlang::enquo(y))),
                                  sep = "\n" )

  micro_set %>%
    dplyr::filter(.data$Taxa == taxa) %>%
    dplyr::group_by(!!!rlang::quos(...)) %>%
    dplyr::mutate(box_groups = interaction(!!!rlang::quos(...), sep = ":")) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$box_groups, y = !!rlang::enquo(y),
                                 fill = .data$box_groups)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(x = xlab, y = ylab, title = main, subtitle = subtitle, fill = legend_title) +
    ggplot2::theme_bw()

}
