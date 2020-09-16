#' @title Create stacked bar charts based on beta binomial model estimates
#' @name bb_bars
#' @description bb_bars takes the output from bb_mods and creates stacked bar charts of the estimated relative abundance for each taxa. The benefit of modeling each taxa before created stacked bar charts is the ability to control for potential confounders. The function will facet wrap interaction terms. Currently, only quant_style = "discrete" can be used for an interaction between two quantitative variables
#' @param modsum The output from bb_mods
#' @param ... The covariate you'd like to plot. Can be an interaction term or main effect, but must be in the models created by bb_mods
#' @param range The range you'd like to plot over for a quantitative variable. Will default to the first and third quartiles
#' @param quant_style "continuous" will plot over the entire range specified; "discrete" will plot only the endpoints of the range specified. "continuous" by default. This option is ignored without a quantitative variable
#' @param top_taxa Only plot X taxa with the highest relative abundance. The rest will be aggregated into an "Other" category
#' @param RA Only plot taxa with a relative abundance higher than X. The rest will be aggregated into an "Other" category
#' @param specific_taxa Character; Plot these specific taxa even if it doesn't meet the top_taxa or RA requirements
#' @param lines Logical; Add outlines around the different taxa colors in the stacked bar charts
#' @param xaxis Labels for the x-axis ticks. Most useful for categorical variables and defaults to the levels of the variable
#' @param main Plot title
#' @param subtitle Subtitle for the plot
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param facet_labels Labels for the facets created for interaction terms
#' @param facet_layout Rearrange the facets created for interaction terms
#' @return Returns a ggplot that you can add geoms to if you'd like
#' @examples
#' data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)
#' otu_tabs = list(Phylum = bpd_phy, Class = bpd_cla,
#' Order = bpd_ord, Family = bpd_fam)
#' set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) %>%
#' filter(day == 7) ## Only including the first week
#'
#' ## Creating beta binomial models on filtered tidy_micro set
#' \donttest{
#' bb_phy <- set %>%
#' otu_filter(ra_cutoff = 0.1, exclude_taxa = c("Unclassified", "Bacteria")) %>%
#' bb_mods(table = "Phylum", bpd1)
#'
#' bb_phy %>%
#' bb_bars(bpd1, top_taxa = 4, xlab = "BPD Severity")
#' }
#' @export
bb_bars <- function(modsum, ..., range, quant_style = c("continuous", "discrete"),
                    top_taxa = 0, RA = 0, specific_taxa, lines = TRUE,
                    xaxis, main, subtitle, xlab, ylab = "Relative Abundance (%)",
                    facet_labels, facet_layout = 1){

  if(!is.double(top_taxa)) stop("top_taxa must be an integer.")
  if(RA < 0 | RA > 100) stop("RA must be between 0 and 100.")
  if(top_taxa > 0 & RA > 0) stop("Can not plot when both top_taxa and RA are greater than 0.")
  if(top_taxa > length(unique(modsum$Convergent_Summary$Taxa))){
    stop("top_taxa must be equal to or less than the number of convergent taxa.")
  }
  if(modsum$Model_Type != "bb_mod") stop("'modsum' should be the output from 'bb_mods'")

  cc <- nb_type(modsum, ...)

  if(length(cc) == 0) stop("Variable specified for bar charts is not in original model.")

  ## Making leaving labels blank optional
  if(missing(main)) main <- NULL ; if(missing(xlab)) xlab <- NULL
  if(missing(specific_taxa)) specific_taxa <- NULL ; if(missing(subtitle)) subtitle <- NULL
  if(missing(facet_labels)) facet_labels <- NULL
  if(facet_layout %nin% c(1,2)) stop("facet_layout must be either 1 or 2")

  tc <- modsum$Model_Coef %>% dplyr::add_count(Taxa)
  if(length(unique(tc$n)) > 1) stop("Repeated Taxa names exist in model's 'Model_Coef'")

  if(cc == "categ") {
    if(missing(xaxis)) xaxis <- sort(unique(modsum$Model_Covs[,cov_str(...)]))

    bb_bars_categ(modsum, ..., top_taxa = top_taxa, lines = lines,
                  RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle, cc = cc)
  } else if(cc == "quant"){

    if(missing(quant_style)) quant_style <- "continuous"
    ## Setting range to 1st and 3rd quartile if not specified
    if(missing(range)) {
      message("Range not specified and will be set to the 1st and 3rd quartile.\n")
      range <- modsum$Model_Covs[,cov_str(...)] %>%
       stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)
    }
    if(!is.numeric(range) | length(range) != 2) stop("range must be a numeric vector of length 2")
    if(missing(xaxis)) xaxis <- as.character(range)

    bb_bars_quant(modsum, ..., quant_style = quant_style, range = range, top_taxa = top_taxa,
                  lines = lines, RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle, cc = cc)
  } else if(cc == "c*c.int"){
    if(missing(xaxis)) xaxis <- sort(unique(modsum$Model_Covs[,cov_str(...)[1]]))

    bb_bars_ccint(modsum, ..., top_taxa = top_taxa,
                  RA = RA, specific_taxa = specific_taxa, main = main, lines = lines,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle,
                  facet_labels = facet_labels, facet_layout = facet_layout, cc = cc)

  } else if(cc == "q*c.int"){
    if(missing(quant_style)) quant_style <- "continuous"
    ## Setting range to 1st and 3rd quartile if not specified
    if(missing(range)) {
      message("Range not specified and will be set to the 1st and 3rd quartile.\n")
      range <- modsum$Model_Covs[,cov_str(...)] %>%
        dplyr::select_if(is.numeric) %>% ## Pull the numeric var
        purrr::simplify() %>% ## Simplify down into a vector
       stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)
    }
    if(!is.numeric(range) | length(range) != 2) stop("range must be a numeric vector of length 2")
    if(missing(xaxis)) xaxis <- as.character(range)

    bb_bars_qcint(modsum, ..., range = range, quant_style = quant_style, top_taxa = top_taxa, lines = lines,
                  RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle,
                  facet_labels = facet_labels, facet_layout = facet_layout, cc = cc)

  } else if(cc == "q*q.int"){

    if(missing(range)){
      message("Ranges not specified and will be set to the 1st and 3rd quartile.\n")

      r1 <- modsum$Model_Covs[,cov_str(...)] %>%
        dplyr::select(1) %>% ## Pull the numeric var
        purrr::simplify() %>% ## Simplify down into a vector
       stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)

      r2 <- modsum$Model_Covs[,cov_str(...)] %>%
        dplyr::select(2) %>% ## Pull the numeric var
        purrr::simplify() %>% ## Simplify down into a vector
       stats::quantile(probs = c(0.25, 0.75), na.rm = TRUE)

      range <- rbind(r1,r2)
    }

    if(all(dim(range) != 2) | !is.matrix(range)) stop("range must be a 2x2 matrix")

    if(missing(xaxis)){ ## xaxis labels when missing
      if(facet_layout == 1) xaxis <- as.character(round(range[1,], 2))
      if(facet_layout == 2) xaxis <- as.character(round(range[2,], 2))
    }

    if(is.null(facet_labels)){ ## facet labels when missing
      if(facet_layout == 1) facet_labels <- as.character(round(range[2,], 2))
      if(facet_layout == 2) facet_labels <- as.character(round(range[1,], 2))
    }
    names(facet_labels) <- c("Low","High")

    bb_bars_qqint(modsum, ..., range = range, top_taxa = top_taxa, lines = lines,
                  RA = RA, specific_taxa = specific_taxa, main = main,
                  xaxis = xaxis, xlab = xlab, ylab = ylab, subtitle = subtitle,
                  facet_labels = facet_labels, facet_layout = facet_layout, cc = cc)
  }
}
