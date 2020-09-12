#' @title Function to make stacked bar charts of taxa relative abundance
#' @name ra_bars
#' @description A function to make stacked bar charts of taxa relative abuncance with the choice to stratify by a variable of interest
#' @param micro_set A tidy_micro data set
#' @param table OTU table you'd like to use when calculating alpha diversity. Your lowest level is recommended
#' @param ... A categorical variable by which you'd like to stratify your relative abundances
#' @param top_taxa Only plot X taxa with the highest relative abundance. The rest will be aggregated into an "Other" category.
#' @param RA Only plot taxa with a relative abundance higher than X. The rest will be aggregated into an "Other" category.
#' @param specific_taxa Plot this specific taxa even if it doesn't meet the top_taxa or RA requirements
#' @param main Plot title
#' @param subtitle Subtitle for the plot
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xaxis Labels for the x-axis ticks. Most useful for categorical variables and defaults to the levels
#' @param lines Logical; Add outlines around the different taxa colors in the stacked bar charts
#' @return Returns a ggplot that you can add geoms to if you'd like
#' @examples
#' data(phy); data(cla); data(ord); data(fam); data(clin)
#'
#' otu_tabs = list(Phylum = phy, Class = cla, Order = ord, Family = fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Full cohort abundance
#' set %>%
#' ra_bars(table = "Family", top_taxa = 10)
#'
#' ## Stratified by variable of interest
#' set %>%
#' ra_bars(table = "Family", bpd1, top_taxa = 10)
#' @export
ra_bars <- function(micro_set, table,  ..., top_taxa = 0, RA = 0, specific_taxa,
                    main, subtitle, ylab, xlab, xaxis, lines = TRUE){
  if(missing(specific_taxa)) specific_taxa <- NULL;
  if(missing(xlab)) xlab <- NULL; if(missing(main)) main <- NULL
  if(missing(subtitle)) subtitle <- NULL
  if(missing(ylab)) ylab <- "Relative Abundance (%)"

  if(table %nin% unique(micro_set$Table)) stop("Specified table is not in supplied micro_set")

  ## No grouping variable, single bar
  if(missing(...)){
    gg <- ra_bar_no_grp(micro_set, table, top_taxa, RA, specific_taxa,
                        main, subtitle, ylab, xlab, xaxis, lines)
  } else{ ## Grouping variable
    gg <- ra_bar_grp(micro_set = micro_set, table = table, ...,
                     top_taxa = top_taxa, RA = RA, specific_taxa = specific_taxa,
                     main = main, subtitle = subtitle, ylab = ylab,
                     xlab = xlab, xaxis = xaxis, lines = lines)
  }

    gg <- gg +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, subtitle = subtitle, x = xlab, y = ylab,
                    fill = "Taxa")

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

  gg
}

ra_bar_no_grp <- function(micro_set, table, top_taxa, RA, specific_taxa,
                            main, subtitle, ylab, xlab, xaxis, lines){
  arr <- micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarise(perc = sum(.data$cts) / sum(unique(.data$Total))*100) %>%
    dplyr::mutate(x = "X") %>%
    dplyr::arrange(dplyr::desc(.data$perc))

  if(RA < 0 | RA > 100) stop("RA must be between 0 and 100")
  if(top_taxa > 0 & RA > 0) stop("Can not aggregate based on both 'top_taxa' and 'RA'")

  ## Aggregating taxa
  if(top_taxa > 0 | RA > 0){
    if(top_taxa > 0){
      arr$ord <- ifelse(arr$Taxa %in% unique(arr$Taxa)[seq(1,top_taxa)], "Top", "Other")
      arr$ord <- ifelse(arr$Taxa %in% specific_taxa, "Top", arr$ord)
    } else if(RA > 0){
      arr$ord <- ifelse(arr$perc >= RA, "Top", "Other")
      arr$ord <- ifelse(arr$Taxa %in% specific_taxa, "Top", arr$ord)
    }

    ar <- rbind(
      arr %>%
        dplyr::filter(.data$ord == "Top") %>%
        as.data.frame,

      cbind(Taxa = "Other",
            arr %>% dplyr::filter(.data$ord == "Other") %>%
              dplyr::summarise(perc = sum(.data$perc))) %>%
        dplyr::mutate(x = "X", ord = "Other") %>%
        as.data.frame
    )

    if(missing(xaxis)) xaxis <- "Cohort Relative Abundance"

    gg <- ar %>%
      Taxa_ord() %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$perc))

    ## Not aggregating taxa
  } else{

    if(missing(xaxis)) xaxis <- "Cohort Relative Abundance"

    gg <- arr %>%
      Taxa_ord() %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$perc))
  }

  gg + ggplot2::scale_x_discrete(labels = xaxis)
}

ra_bar_grp <- function(micro_set, table,  ..., top_taxa, RA, specific_taxa,
                       main, subtitle, ylab, xlab, xaxis, lines){
  arr <- micro_set %>%
    dplyr::filter(.data$Table == table) %>%
    dplyr::group_by(!!!rlang::quos(...)) %>%
    dplyr::summarise(tot = sum(.data$cts)) %>%
    dplyr::left_join(micro_set %>%
                       dplyr::filter(.data$Table == table) %>%
                       dplyr::group_by(.data$Taxa,!!!rlang::quos(...)) %>%
                       dplyr::summarise(count = sum(.data$cts)), by = by_fun(!!!rlang::quos(...))) %>%
    dplyr::group_by(.data$Taxa, !!!rlang::quos(...)) %>%
    dplyr::summarise(perc = (.data$count/.data$tot)*100) %>%
    dplyr::arrange(dplyr::desc(.data$perc))

  if(RA < 0 | RA > 100) stop("RA must be between 0 and 100")
  if(top_taxa > 0 & RA > 0) stop("Can not aggregate based on both 'top_taxa' and 'RA'")

  ## Aggregating taxa
  if(top_taxa > 0 | RA > 0){
    if(top_taxa > 0){
      arr$ord <- ifelse(arr$Taxa %in% unique(arr$Taxa)[seq(1,top_taxa)], "Top", "Other")
      arr$ord <- ifelse(arr$Taxa %in% specific_taxa, "Top", arr$ord)
    } else if(RA > 0){
      arr$ord <- ifelse(arr$perc >= RA, "Top", "Other")
      arr$ord <- ifelse(arr$Taxa %in% specific_taxa, "Top", arr$ord)
    }
    ar <- rbind(
      arr %>%
        dplyr::filter(.data$ord=="Top") %>%
        dplyr::select(.data$Taxa, !!!rlang::quos(...), .data$perc) %>%
        as.data.frame,

      cbind(Taxa = "Other",
            arr %>% dplyr::filter(.data$ord == "Other") %>%
              dplyr::group_by(!!!rlang::quos(...)) %>%
              dplyr::summarize(perc = sum(.data$perc)) %>%
              as.data.frame
      )
    )
    ar$plot_grp <- ar %>%
      dplyr::select(!!!rlang::quos(...)) %>%
      interaction

    gg <- ar %>%
      Taxa_ord() %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$plot_grp, y = .data$perc))

  } else{ ## No Aggregation
    if(length(cov_str(...)) == 1){
      names(arr)[2] <- "plot_grp"

      arr %<>%
        dplyr::filter(!is.na(.data$plot_grp)) %>%
        dplyr::mutate(plot_grp = factor(.data$plot_grp))

    } else if(length(cov_str(...)) > 1){
      arr %<>%
        dplyr::mutate(plot_grp = interaction(!!!rlang::quos(...)))
    }

    gg <- arr %>%
      dplyr::ungroup() %>%
      Taxa_ord() %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$plot_grp, y = .data$perc))
  }
  if(missing(xaxis)){
    xaxis <- micro_set %>%
      dplyr::select(!!!rlang::quos(...)) %>%
      interaction %>%
      levels
  }
  gg + ggplot2::scale_x_discrete(labels = xaxis)
}



