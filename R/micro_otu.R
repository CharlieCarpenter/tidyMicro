#' @title Extract OTU table from tidyMicro set
#' @name micro_otu
#' @description A simple wrapper to extract an OTU table from a tidyMicro set
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param taxa_info The taxa info to pull
#' @return A tibble
#' @examples
#' data(mrsa_gen); data(mrsa_clin)
#'
#' ## Creating tidyMicro set
#' set <- tidy_micro(otu_tabs = mrsa_gen, tab_names = "Genus", clinical = mrsa_clin)
#'
#' ## Filtering out unwanted OTUs
#' filt.set <- otu_filter(set, prev_cutoff = 1, ra_cutoff = 1, filter_summary = FALSE)
#'
#' ## Extract filtered OTU table
#' filt.otu.cts <- micro_otu(filt.set, table = "Genus")
#'
#' ## Extract filtered relative abundances table
#' filt.otu.ra <- micro_otu(filt.set, table = "Genus", taxa_info = ra)
#'
#' @export
micro_otu <- function(micro_set, table, taxa_info = cts){
  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::select(.data$Lib, .data$Taxa, !!enquo(taxa_info)) %>%
    dplyr::filter(!is.na(.data$Taxa)) %>%
    tidyr::pivot_wider(names_from = .data$Lib, values_from = !!enquo(taxa_info))
}
