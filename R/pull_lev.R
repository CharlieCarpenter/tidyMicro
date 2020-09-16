#' @title Shorten long taxa names
#' @name pull_lev
#' @description \code{pull_lev} is designed to work with Taxa names that contain multiple taxonomy levels. The function splits Taxa names using \code{\link[stringr]{str_split}}, then grabs the string in position \code{l} or higher. If nothing exists at or past position \code{l}, the function will add "Unclassified" to the last position. The function skips over any Taxa named "Other".
#' @param Taxa Taxa names with taxonomic rank separated by some character
#' @param l Numeric. Position / taxonomic rank. Best if it is the highest rank in Taxa names
#' @param sep Character to split
#' @return A character string
#' @examples
#' data(bpd_fam); data(bpd_clin)
#' set <- tidy_micro(otu_tabs = bpd_fam, clinical = bpd_clin, tab_names = "Family",
#' prev_cutoff = 5, ra_cutoff = 1, exclude_taxa = "Bacteria",
#' count_summary = FALSE, filter_summary = FALSE) %>%
#' filter(day == 7)
#'
#' ## Taxa names
#' pull_lev(unique(set$Taxa), 4)
#'
#' ## Most useful to clean up plots ##
#'
#' # Hard to read!
#' (ch <- cor_heatmap(set, table = "Family", weight, gestational_age))
#'
#' # Much better
#' ch$data$Taxa %<>% pull_lev(4)
#' ch
#' @export
pull_lev <- function(Taxa, l, sep = "/"){
  Taxa %>% stringr::str_split(sep) %>%
    lapply(function(t, l.){
      if(length(t) >= l) t <- t[length(t)]
      else if(any(t == "Other")) t <- t
      else if(length(t) < l) t <- paste("Unclassified", t[length(t)])

      t
    }, l. = l) %>%
    unlist
}
