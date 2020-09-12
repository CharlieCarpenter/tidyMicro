#' @title A function to aggregate low prevalence, abundance, or unwanted taxa together
#' @name otu_filter
#' @description Will take a tidi_micro set and aggregate the raw counts of taxa with a low prevalence and/or abundance into a new "Other" taxa. Can also find specific taxa you'd like to include in the "Other" taxa counts. Once the counts are aggregated taxa relative abundance, centered log ratio (CLR) transformations, and presence will be recalculated. This recalculation will only change the "Other" category
#' @param micro_set A tidy_micro data set
#' @param prev_cutoff Minimum percent of subjects with OTU counts above 0
#' @param ra_cutoff At leat one subject must have RA above this subject
#' @param exclude_taxa A character vector of OTU names that you would like filter into your "Other" category
#' @param filter_summary Logical; print out summaries of filtering steps
#' @details \eqn{\frac{1}{Total}}{1/Total} will be added to each taxa count for CLR tranformations in order to avoid issues with log(0)
#' @return Returns a tidy_micro set
#' @author Charlie Carpenter and Dan Frank
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#'
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' filter_set <- set %>%
#' otu_filter(prev_cutoff = 5, ## 5% of subjects must have this bug, or it is filtered
#'   ra_cutoff = 1, ## At least 1 subject must have RA of 1, or it is filtered
#'   exclude_taxa = c("Unclassified", "Bacteria") ## Unclassified taxa we don't want
#' )
#' @export
otu_filter <- function(micro_set, prev_cutoff = 0, ra_cutoff = 0, exclude_taxa = NULL,
                       filter_summary = T){

  ## Table, Lib, Taxa name, meta data left over
  met <- micro_set %>% dplyr::select(-dplyr::matches("Table"), -dplyr::matches("Taxa"),
                                    -dplyr::matches("Total"), -dplyr::matches("bin"),
                                    -dplyr::matches("cts"), -dplyr::matches("clr"),
                                    -dplyr::matches("ra")
  ) %>%
    dplyr::distinct(.data$Lib, .keep_all = T)

  ## Applies function to each table
  micro_set %>%
    plyr::ddply(~ .data$Table, function(set, meta,
                                               prev_cutoff, ra_cutoff,
                                               exclude_taxa){

      ## Pulls out Lib, Taxa, and counts and spreads to original otu structure
      otu <- set %>%
        dplyr::select(.data$Lib, .data$Taxa, .data$cts) %>%
        dplyr::filter(!is.na(.data$Taxa)) %>%
        tidyr::pivot_wider(names_from = .data$Taxa, values_from = .data$cts)

      Lib <- otu$Lib; total <- rowSums(otu[,-1]); rr <- unique(as.character(set$Table))

      ## applies all filters specified
      if(filter_summary){
        otu.f <- mul_filter(otu[,-1], prev_cutoff, ra_cutoff, total = total,
                            ex = exclude_taxa, rr = rr)
      } else{
        otu.f <- suppressMessages(
          mul_filter(otu[,-1], prev_cutoff, ra_cutoff, total = total,
                     ex = exclude_taxa, rr = rr)
        )
      }

      ## recalculates ra, clr, and bin (just to recalculate the "Other")
      ra <- my_ra(otu.f, total) %>% dplyr::mutate(Lib = Lib)
      clr <- my_clr(otu.f, total) %>% dplyr::mutate(Lib = Lib)
      bin <- my_bin(otu.f, total) %>% dplyr::mutate(Lib = Lib)

      ## Recalculating totals (for "Other" category) and remaking cts
      Tot <- data.frame(Total = rowSums(otu.f), Lib = Lib)
      otu.f$Lib <- Lib ## filtered cts

      ## Melted data
      m.ra <- ra %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "ra")
      m.clr <- clr %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "clr")
      m.bin <- bin %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "bin")
      m.cts <- otu.f %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "cts")

      ## Joining all counts and transformations
      long_OTU <- suppressWarnings(dplyr::left_join(m.bin, m.cts, by = c("Lib", "Taxa")) %>%
                                     dplyr::left_join(m.clr, by = c("Lib", "Taxa")) %>%
                                     dplyr::left_join(m.ra, by = c("Lib", "Taxa")) %>%
                                     dplyr::left_join(Tot, by = "Lib") %>%
                                     dplyr::mutate(Table = rr) %>%
                                     dplyr::left_join(meta, by = "Lib")
      )

      long_OTU
    }, meta = met, prev_cutoff, ra_cutoff, exclude_taxa) %>%

    ## Reordering to original order
    dplyr::filter(!is.na(.data$Total)) %>%
    dplyr::select(.data$Table, .data$Lib, .data$Taxa, .data$Total,
                  .data$bin, .data$cts, .data$clr, .data$ra, dplyr::everything(),
                  -dplyr::matches(".data$Total"))
}
