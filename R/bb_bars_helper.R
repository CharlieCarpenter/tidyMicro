###############################
##
## Project: tidyMicro
##
## Purpose: bb_bars helper functions
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-01-29
##
## ---------------------------
## Notes:
##
##
## ---------------------------


bb_bars_categ <- function(modsum, ..., top_taxa, RA, specific_taxa,
                          xaxis, xlab, ylab, main, subtitle, lines, cc){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs,...)
  }

  msum <- modsum$Model_Coef %>%
    dplyr::filter(.data$Cov_Type == "Intercept" |
                    stringr::str_detect(.data$Coef, cov_str(...))) %>%
    plyr::ddply(~ .data$Taxa, bb_categ_est) %>%
    dplyr::mutate(Coef = factor(.data$Coef, levels = unique(.data$Coef)))%>%
    Taxa_ord()

  ord <- msum %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarize(mean.bb = mean(.data$bb.avg)) %>%
    dplyr::arrange(dplyr::desc(.data$mean.bb))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style = NULL, ...) %>%
      bb_categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab,
                    subtitle=subtitle, lines = lines)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$bb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style = NULL, ...) %>%
      bb_categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab,
                    subtitle=subtitle, lines = lines)
  } else {

    SS <- msum %>%
      dplyr::group_by(.data$Coef) %>%
      dplyr::summarize(ss = sum(.data$bb.avg))

    if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                       "and will be rescaled to equal 100."),
                                 call. = FALSE)
    if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                      "and will be rescaled to equal 100."),
                                call. = FALSE)
    msum %>%
      plyr::ddply(~ .data$Coef,
                  function(set){
                    set %<>% dplyr::mutate(bb.avg = 100*.data$bb.avg/sum(.data$bb.avg))
                  }) %>%
      bb_categ_bars(xaxis=xaxis, main=main, xlab=xlab,
                    ylab=ylab, subtitle=subtitle, lines = lines)
  }
}

bb_bars_quant <- function(modsum, ..., range, quant_style, top_taxa, RA,
                          specific_taxa, main, xaxis, xlab, ylab, subtitle, lines, cc){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs, ...)
  }

  if(quant_style == "continuous"){
    RR <- seq(range[1], range[2], length.out = 10)

    msum <- modsum$Model_Coef %>%
      dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...)))

    ## Estimates for taxa along RR range
    Est <- exp(tcrossprod(msum$Estimate, RR) + msum$Intercept) /
      (1 + exp(tcrossprod(msum$Estimate, RR) + msum$Intercept)) * 100
    colnames(Est) <- RR

    msum %<>%
      cbind(Est) %>%
      dplyr::select(.data$Taxa, 6:15) %>%
      tidyr::pivot_longer(-.data$Taxa, names_to = "Quant",
                          values_to = "Value") %>%
      Taxa_ord()

    ord <- msum %>%
      dplyr::group_by(.data$Taxa) %>%
      dplyr::summarise(bb.avg = mean(.data$Value)) %>%
      dplyr::arrange(dplyr::desc(.data$bb.avg))
  }
  if(quant_style == "discrete"){
    msum <- modsum$Model_Coef %>%
      dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...))) %>%
      dplyr::mutate(bb.avg = exp(.data$Intercept + .data$Estimate*mean(.data$range))/ (1+exp(.data$Intercept + .data$Estimate*mean(.data$range)))*100,
                    bb.low = exp(.data$Intercept + .data$Estimate*.data$range[1] + log(100))/(1+exp(.data$Intercept + .data$Estimate*.data$range[1] + log(100)))*100,
                    bb.up = exp(.data$Intercept + .data$Estimate*.data$range[2] + log(100))/(1+exp(.data$Intercept + .data$Estimate*.data$range[2] + log(100)))*100
      ) %>%
      Taxa_ord()

    ord <- msum %>%
      dplyr::arrange(dplyr::desc(.data$bb.avg))
  }

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_quant_bars(range = range, quant_style = quant_style, xaxis = xaxis, main = main,
                    xlab = xlab, ylab = ylab, subtitle=subtitle, lines = lines)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$bb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum,cc,quant_style,...) %>%
      bb_quant_bars(range = range, quant_style = quant_style, xaxis = xaxis, main = main,
                    xlab = xlab, ylab = ylab, subtitle=subtitle, lines = lines)
  } else
    msum %>%
    bb_quant_bars(range = range, quant_style = quant_style, xaxis = xaxis, main = main,
                  xlab = xlab, ylab = ylab, subtitle=subtitle, lines = lines)

}

bb_bars_ccint <- function(modsum, ..., top_taxa, RA, specific_taxa, xaxis, lines,
                          xlab, ylab, main, subtitle, facet_labels, facet_layout, cc){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs, ...)
  }

  l1 <- levels(modsum$Model_Covs[, cov_str(...)[1]]) ; l2 <- levels(modsum$Model_Covs[, cov_str(...)[2]])
  if(any(l1%in%l2)) stop("Factor variables can not share levels.")


  msum <- suppressWarnings(
    modsum$Model_Coef %>%
      dplyr::filter(.data$Cov_Type == "Intercept" |
                      stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                      stringr::str_detect(.data$Coef, cov_str(...)[2]) ) %>%
      plyr::ddply( ~ .data$Taxa, bb_ccint_est, l1, l2)
  ) %>%
    Taxa_ord()

  ord <- msum %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarise(bb.avg = mean(.data$bb.avg)) %>%
    dplyr::arrange(dplyr::desc(.data$bb.avg))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style = NULL, ...) %>%
      bb_ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = NULL, facet_layout = facet_layout, ...)

  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$bb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style = NULL, ...) %>%
      bb_ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = facet_labels, facet_layout = facet_layout, ...)
  } else
    msum %>%
    bb_ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                  lines = lines, facet_labels = facet_labels, facet_layout = facet_layout, ...)
}

bb_bars_qcint <- function(modsum, ..., range, quant_style, top_taxa, lines, RA,
                          specific_taxa, main, xaxis, xlab, ylab, subtitle, facet_labels,
                          facet_layout, cc){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs, ...)
  }

  msum <- modsum$Model_Coef %>%
    dplyr::filter(.data$Coef == "(Intercept)" |
                    stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                    stringr::str_detect(.data$Coef, cov_str(...)[2])) %>%
    plyr::ddply( ~ .data$Taxa, bb_qcint_est,
                modsum$Model_Covs, quant_style, range, ...) %>%
    Taxa_ord()

  if(quant_style == "continuous"){
    ord <- msum %>%
      dplyr::group_by(.data$Taxa) %>%
      dplyr::summarise(bb.avg = mean(.data$Value)) %>%
      dplyr::arrange(dplyr::desc(.data$bb.avg))
  }

  if(quant_style == "discrete"){

    ord <- msum %>%
      dplyr::mutate(nb.avg = tapply(.data$Est, .data$Taxa, mean) %>%
                      rep(each = nrow(msum)/length(unique(msum$Taxa)))) %>%
      dplyr::arrange(dplyr::desc(.data$nb.avg))
  }

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_qcint_bars(m_cov = modsum$Model_Covs,
                    quant_style = quant_style, xaxis = xaxis, main = main, xlab = xlab,
                    ylab = ylab, subtitle = subtitle, lines = lines,
                    facet_labels = facet_labels, facet_layout = facet_layout, ...)

  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_qcint_bars(m_cov = modsum$Model_Covs,
                    quant_style = quant_style, xaxis = xaxis, main = main, xlab = xlab,
                    ylab = ylab, subtitle = subtitle, lines = lines,
                    facet_labels = facet_labels, facet_layout = facet_layout, ...)
  } else
    msum %>%
    bb_qcint_bars(m_cov = modsum$Model_Covs,
                  quant_style = quant_style, xaxis = xaxis, main = main, xlab = xlab,
                  ylab = ylab, subtitle = subtitle, lines = lines,
                  facet_labels = facet_labels, facet_layout = facet_layout, ...)
}

bb_bars_qqint <- function(modsum, range, ..., top_taxa, RA, specific_taxa, xaxis, lines,
                          xlab, ylab, main, subtitle, facet_labels, facet_layout, cc){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs,...)
  }

  msum <- modsum$Model_Coef %>%
    dplyr::filter(.data$Cov_Type == "Intercept" |
                    stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                    stringr::str_detect(.data$Coef, cov_str(...)[2])) %>%
    plyr::ddply(~ .data$Taxa, bb_qqint_est, range, ...) %>%
    Taxa_ord()

  ord <- msum %>%
    dplyr::arrange(dplyr::desc(.data$nb.avg))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style = NULL, ...) %>%
      bb_qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)

  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style = NULL, ...) %>%
      bb_qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)

  } else
    msum %>%
    bb_qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                  lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)
}
