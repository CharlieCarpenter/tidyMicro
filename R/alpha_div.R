#' @title Alpha Diversity Calculations for tidy_micro
#' @name alpha_div
#' @description A wrapper function to calculate Sobs, Choa1, Goods, Shannon's diversity and evenness, and Simpson's diversity and evenness alpha diversities for your micro_set. Estimates are calculated based on rarefied bootstrapped samples
#' @param micro_set A tidy_micro data set
#' @param table OTU table of interest
#' @param iter The number of bootstrap resamples used for estimation
#' @param min_depth Filter out libraries with sequencing depth (Total) below min_depth
#' @param min_goods Filter out libraries Good's coverage below min_goods
#' @details If you have multiple otu tables, you can specify the table you'd like to use to calculate your alpha diversities using the \code{table} option. We highly recommend using the lowest taxonomic rank available to calculate your alpha diversity. If you would like to calculate alpha diversities for each otu table in your micro_set, you can leave the \code{table} option as \code{NULL} and the function will calculate the alpha diversity for each table. The function will append the estimated alpha diversities to the tidy_micro supplied. The alpha diversity columns will be just before your clinical data. Since alpha diversity is estimated for each individual library (Lib), it will be repeated within each taxa block.
#' @return A tidy_micro set with alpha diversity columns added in to the left of clinical data
#' @note Be aware of your minimal sequencing depth as this will be the size of all bootstrapped resamples (rarefied).
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#'
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' \donttest{
#' ## calculate alpha diversity for every table
#' set_alpha <- set %>% alpha_div(min_depth = 5000, min_goods = 90)
#'
#' ## calculate alpha diversity for a specific table
#' set_fam_alpha <- set %>% alpha_div(table = "Family", min_depth = 5000, min_goods = 90)
#' }
#' @export
alpha_div <- function(micro_set, table = NULL, iter = 100, min_depth = 0, min_goods = 0){

  if(min(micro_set$Total) < 10000){
    warning(paste0("The minimum Library size is ", min(micro_set$Total), ".\n"))
  }

  if(!is.null(table)){ ## Specific table
    if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")
    micro_set %<>%
      dplyr::filter(.data$Table == table) %>%
      dplyr::mutate(Table = factor(.data$Table, levels = unique(.data$Table)))
  }

  ## Calculating alpha diversities
  tidy_alpha <- micro_set %>%
    ## Filtering out based on minimum sequencing depth
    dplyr::filter(.data$Total > min_depth) %>%
    dplyr::left_join(goods(micro_set, iter), by = c("Table", "Lib")) %>%
    ## Filtering out based on minimum Good's coverage
    dplyr::filter(.data$Goods > min_goods) %>%
    dplyr::left_join(sobs(micro_set, iter), by = c("Table", "Lib")) %>%
    dplyr::left_join(chao1(micro_set, iter), by = c("Table", "Lib")) %>%
    dplyr::left_join(shannonH(micro_set, iter), by = c("Table", "Lib")) %>%
    dplyr::left_join(shannonE(micro_set, iter), by = c("Table", "Lib")) %>%
    dplyr::left_join(simpsonD(micro_set, iter), by = c("Table", "Lib")) %>%
    dplyr::left_join(simpsonE(micro_set, iter), by = c("Table", "Lib")) %>%
    dplyr::select(.data$Table, .data$Lib, .data$Taxa, .data$Total,
                  .data$bin, .data$cts, .data$clr, .data$ra,
                  .data$Goods, .data$Sobs, .data$Chao1, .data$ShannonE,
                  .data$ShannonH, .data$SimpsonD, .data$SimpsonE, dplyr::everything()) %>%
    dplyr::distinct()

  ## Messages about dropped libraries if filters were applied
  if(min_depth > 0 | min_goods > 0){
      dropped_lib <- unique(micro_set$Lib)[unique(micro_set$Lib) %nin% unique(tidy_alpha$Lib)]

      if(length(dropped_lib) > 0){
        message("Dropped ", length(dropped_lib), " libraries based on sequencing depth or Good's coverage.\n")

        if(length(dropped_lib) == 1) message("Library ", dropped_lib, " was dropped.")

        if(length(dropped_lib) == 2){
          message("Libraries ", dropped_lib[1]," and ", dropped_lib[2], " were dropped.")
        }
        if(length(dropped_lib) > 2){
          drp <- dropped_lib[1]
          for(i in 2:(length(dropped_lib)-1)) drp <- paste(drp, dropped_lib[i], sep=", ")
          message("Libraries ", drp,", and ", dropped_lib[length(dropped_lib)], " were dropped.")
        }
      } else message("No libraries dropped based on sequencing depth or Good's coverage.")
    }

  tidy_alpha
}
