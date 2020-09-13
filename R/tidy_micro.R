#' @title A function to merge multiple OTU tables and clinical data into a "tidy" format
#' @name tidy_micro
#' @description A function to take any number of OTU tables (or other sequencing data tables), calculate taxa prevalence, relative abundance, and a CLR transformation, and finally merges clinical data
#' @param otu_tabs A single table or list of metagenomic sequencing data. Tables should have a first column of OTU Names and following columns of OTU counts. Column names should be sequencing library names
#' @param clinical Sequencing level clinical data. Must have a column with unique names for library (sequencing ID)
#' @param tab_names names for otu_tabs. These will become the "Tables" column. It is also an option to simply name the OTU tables in the list supplied to otu_tabs
#' @param prev_cutoff A prevalence cutoff where *X* percent of libraries must have this taxa or it will be included in the "Other" category
#' @param ra_cutoff A relative abundance (RA) cutoff where at least one library must have a RA above the cutoff or the taxa will be included in the "Other" category
#' @param exclude_taxa A character vector used to specify any taxa that you would like to included in the "Other" category. Taxa specified will be included in "Other" for every OTU table provided
#' @param library_name The column name containing sequencing library names. Should match with column names of supplied OTU tables (after first column)
#' @param complete_clinical Logical; only include columns from OTU tables who's library name is in clinical data
#' @param filter_summary Logical; print out summaries of filtering steps. Ignored \code{prev_cutoff}, \code{ra_cutoff}, and \code{exclude_taxa} are all left as default values
#' @param count_summary Logical: print out summary of unique library names and sequencing depth
#' @details Column names of the OTU tables must be the same for each table, and these should be the the library names inside of your clinical. Please see the \link{vignette} for a detailed description.
#'
#' The CLR transformation adds (1 / sequencing depth) to each OTU count for each library before centering and log transforming in order to avoid issues with 0 counts.
#'
#' The list of OTU tables are split, manipulated, and stacked into a data frame using the \code{\link[plyr]{ldply}} function from the \pkg{plyr} package. Names of OTU tables supplied will be the name of their "Table" in the final tidy_micro set
#' @return A data.frame in the tidy_micro format
#' @author Charlie Carpenter
#' @import magrittr tidyverse
#' @importFrom rlang .data
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#'
#' ## Multiple OTU tables with named list
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin)
#'
#' ## Multiple OTU tables with unnamed list
#' unnamed_tabs <- list(bpd_phy, bpd_cla, bpd_ord, bpd_fam)
#' set <- tidy_micro(otu_tabs = unnamed_tabs,
#' tab_names = c("Phylum", "Class", "Order", "Family"), clinical = bpd_clin)
#'
#' ## Single OTU table
#' set <- tidy_micro(otu_tabs = bpd_cla, clinical = bpd_clin, tab_names = "Class")
#'
#' ## Filtering out low abundance or uninteresting taxa right away
#' ## WARNING: Only do this if you do not want to calculate alpha diversities with this tidy_micro set
#'
#' filter_set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin,
#'               prev_cutoff = 5, ## 5% of libraries must have this bug, or it is filtered
#'               ra_cutoff = 1, ## At least 1 libraries must have RA of 1, or it is filtered
#'               exclude_taxa = c("Unclassified", "Bacteria") ## Unclassified taxa we don't want
#'               )
#' @export
tidy_micro <- function(otu_tabs, clinical, tab_names,
                       prev_cutoff = 0.0, ra_cutoff = 0.0, exclude_taxa = NULL,
                       library_name = "Lib", complete_clinical = TRUE, filter_summary = TRUE,
                       count_summary = TRUE){

  if(library_name %nin% names(clinical)) stop("Must provide 'library_name' matching column name of sequencing IDs from clinical data")

  names(clinical)[names(clinical) == library_name] <- "Lib"
  if(library_name != "Lib") warning(paste0("Renaming '", library_name, "' to 'Lib'") )

  if(is.data.frame(otu_tabs) | is.matrix(otu_tabs)){ ## One OTU table
    if(length(tab_names) != 1) {
      stop("tab_names must have length of 1 if only one OTU table supplied.")}

    l_otu <- suppressWarnings( otu_tabs %>%
                                 mul_otu_long(.meta = clinical) %>%
                                 dplyr::mutate(Table = tab_names)
    )
  } else if(is.list(otu_tabs)){ ## Multiple OTU tables

    if(missing(tab_names) & is.null(names(otu_tabs))){
      stop("otu_tabs must have names, or names must be supplied through tab_names.")
    }

    if(is.null(names(otu_tabs))) names(otu_tabs) <- tab_names

    l_otu <- suppressWarnings(otu_tabs %>%
                                plyr::ldply(mul_otu_long, .meta = clinical, .id = "Table"))
  } else stop("otu_tabs must be a single data.frame/matrix or list of data.frames/matrices.")

  l_otu %<>%
    dplyr::distinct() %>%
    otu_filter(prev_cutoff = prev_cutoff, ra_cutoff = ra_cutoff,
               exclude_taxa = exclude_taxa, filter_summary = filter_summary) %>%
    dplyr::select(.data$Table, dplyr::everything())

  ## Checking if all sequence libraries have clinical data and visa versa
  if(complete_clinical){
    if(!all(l_otu$Lib %in% clinical$Lib) | !all(clinical$Lib %in% l_otu$Lib)){
      warning("Warning: Some libraries in OTU table and clinical data don't match. Only libraries contained in each will be kept.\n")

      l_otu %<>% dplyr::filter(.data$Lib %in% clinical$Lib)
    }
  }

  if(count_summary){
    ## Counts subjects
    message("Contains ",length(unique(l_otu$Lib))," libraries from OTU files.\n")

    ss <- l_otu %>%
      dplyr::distinct(.data$Lib, .keep_all = T) %>%
      dplyr::pull(.data$Total) %>%
      summary

    ## Printing out a summary of the sequencing depth
    message("Summary of sequencing depth:"); print(ss)
  }

  l_otu %<>%
    dplyr::mutate(Table = factor(.data$Table), Taxa = factor(.data$Taxa))

  l_otu[names(l_otu) != ".data$Table"]
}
