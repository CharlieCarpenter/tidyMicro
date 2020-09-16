###############################
##
## Project: tidyMicro
##
## Purpose: Helper functions for bb_bars
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2019-10-21
##
## ---------------------------
## Notes:
##
##
## ---------------------------
## HEADER ----

# Helpful functions

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
    dplyr::mutate(Coef = factor(.data$Coef,
                                levels = unique(.data$Coef))) %>%
    Taxa_ord()

  ord <- msum %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarize(mean.bb = mean(.data$bb.avg)) %>%
    dplyr::arrange(dplyr::desc(.data$mean.bb))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab, subtitle=subtitle, lines = lines)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$bb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab,
                    subtitle=subtitle, lines = lines)
  } else {

    SS <- msum %>% dplyr::group_by(Coef) %>% dplyr::summarize(ss = sum(bb.avg))
    if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                       "and will be rescaled to equal 100."),
                                 call. = FALSE)
    if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                      "and will be rescaled to equal 100."),
                                call. = FALSE)
    msum %>%
      plyr::ddply(~ .data$Coef,
                  .fun = function(set){
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
      dplyr::filter(stringr::str_detect(Coef, cov_str(...)))

    ## Estimates for taxa along RR range
    Est <- exp(tcrossprod(msum$Estimate, RR) + msum$Intercept) /
      (1 + exp(tcrossprod(msum$Estimate, RR) + msum$Intercept)) * 100
    colnames(Est) <- RR

    msum %<>%
      cbind(Est) %>%
      dplyr::select(Taxa, 6:15) %>%
      tidyr::gather(Quant, Value, -Taxa) %>%
      Taxa_ord()

    ord <- msum %>%
      dplyr::group_by(Taxa) %>%
      dplyr::summarise(bb.avg = mean(Value)) %>%
      dplyr::arrange(dplyr::desc(bb.avg))
  }
  if(quant_style == "discrete"){
    msum <- modsum$Model_Coef %>%
      dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...))) %>%
      dplyr::mutate(bb.avg = exp(.data$Intercept + .data$Estimate*mean(range))/
                      (1+exp(.data$Intercept + .data$Estimate*mean(range)))*100,

                    bb.low = exp(.data$Intercept + .data$Estimate*range[1] )/
                      (1+exp(.data$Intercept + .data$Estimate*range[1]))*100,

                    bb.up = exp(.data$Intercept + .data$Estimate*range[2])/
                      (1+exp(.data$Intercept + .data$Estimate*range[2]))*100
      ) %>%
      Taxa_ord()

    ord <- msum %>%
      dplyr::arrange(dplyr::desc(bb.avg))
  }

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum,cc,quant_style, ...) %>%
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
      dplyr::filter(Cov_Type == "Intercept" |
                      stringr::str_detect(Coef, cov_str(...)[1]) |
                      stringr::str_detect(Coef, cov_str(...)[2]) ) %>%
      plyr::ddply(~ .data$Taxa, bb_ccint_est, l1, l2)
  ) %>%
    Taxa_ord()

  ord <- msum %>%
    dplyr::group_by(Taxa) %>%
    dplyr::summarise(bb.avg = mean(bb.avg)) %>%
    dplyr::arrange(desc(bb.avg))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = NULL, facet_layout = facet_layout, ...)

  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$bb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
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
    dplyr::filter(Coef == "(Intercept)" |
                    stringr::str_detect(Coef, cov_str(...)[1]) |
                    stringr::str_detect(Coef, cov_str(...)[2])) %>%
    plyr::ddply(~ .data$Taxa, bb_qcint_est, modsum$Model_Covs, quant_style, range, ...) %>%
    Taxa_ord()

  if(quant_style == "continuous"){
    ord <- msum %>%
      dplyr::group_by(Taxa) %>%
      dplyr::summarise(bb.avg = mean(Value)) %>%
      dplyr::arrange(desc(bb.avg))
  }

  if(quant_style == "discrete"){

    ord <- msum %>%
      dplyr::mutate(nb.avg = tapply(.data$Est, .data$Taxa, mean) %>%
                      rep(each = nrow(msum)/length(unique(.data$Taxa)))) %>%
      dplyr::arrange(desc(.data$nb.avg))
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
    dplyr::filter(Cov_Type == "Intercept" |
                    stringr::str_detect(Coef, cov_str(...)[1]) |
                    stringr::str_detect(Coef, cov_str(...)[2])) %>%
    plyr::ddply(~ .data$Taxa, bb_qqint_est, range, ...) %>%
    Taxa_ord()

  ord <- msum %>%
    dplyr::arrange(dplyr::desc(nb.avg))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)

  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    bb_agg_bars(msum, cc, quant_style, ...) %>%
      bb_qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                    lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)

  } else
    msum %>%
    bb_qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                  lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)
}

#### Aggregation functions ####

bb_agg_bars <- function(msum, cc, quant_style, ...){
  if(cc == "categ") {
    agg <- bb_agg_categ(msum, ...)
  } else if(cc == "quant"){
    agg <- bb_agg_quant(msum, quant_style, ...)
  } else if(cc == "c*c.int"){
    agg <- bb_agg_ccint(msum, ...)
  } else if(cc == "q*c.int"){
    agg <- bb_agg_qcint(msum, quant_style, ...)
  } else if(cc == "q*q.int"){
    agg <- bb_agg_qqint(msum, ...)
  }
  agg
}

bb_agg_categ <- function(msum, ...){
  agg <- rbind(
    msum %>%
      dplyr::filter(.data$ord=="Top") %>%
      dplyr::select(.data$Taxa, .data$Coef, .data$bb.avg),

    cbind(Taxa = "Other",
          msum %>% dplyr::filter(.data$ord == "Other") %>%
            dplyr::group_by(.data$Coef) %>%
            dplyr::summarize(bb.avg = sum(.data$bb.avg))
    )
  )
  SS <- agg %>% dplyr::group_by(.data$Coef) %>%
    dplyr::summarize(ss = sum(.data$bb.avg))

  if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                     "and will be rescaled to equal 100."),
                               call. = FALSE)
  if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                    "and will be rescaled to equal 100."),
                              call. = FALSE)

  agg %<>% plyr::ddply(~ .data$Coef,
                       function(set){
                         set %<>%
                           dplyr::mutate(bb.avg = 100*.data$bb.avg/sum(.data$bb.avg))
                       })
  agg
}

bb_agg_quant <- function(msum, quant_style, ...){
  if(quant_style == "discrete"){
    agg <- rbind(
      msum %>%
        dplyr::filter(.data$ord=="Top") %>%
        dplyr::select(.data$Taxa, .data$Coef,
                      .data$bb.low, .data$bb.up),

      cbind(Taxa = "Other",
            msum %>% dplyr::filter(.data$ord == "Other") %>%
              dplyr::group_by(.data$Coef) %>%
              dplyr::summarize(bb.low = sum(.data$bb.low), bb.up = sum(.data$bb.up))
      )
    )

    if(sum(agg$bb.up) > 105) warning(paste("Maximum column sum is", round(sum(agg$bb.up)),
                                           "and will be rescaled to equal 100."),
                                     call. = FALSE)
    if(sum(agg$bb.low) < 95) warning(paste("Minimum column sum is", round(sum(agg$bb.low)),
                                           "and will be rescaled to equal 100."),
                                     call. = FALSE)

    agg %<>%
      dplyr::mutate(bb.low = 100*.data$bb.low/sum(.data$bb.low),
                    bb.up = 100*.data$bb.up/sum(.data$bb.up))
  }
  if(quant_style == "continuous"){
    agg <- rbind(
      msum %>%
        dplyr::filter(.data$ord == "Top") %>%
        dplyr::select(.data$Taxa, .data$Quant, .data$Value),

      cbind(Taxa = "Other",
            msum %>% dplyr::filter(.data$ord == "Other") %>%
              dplyr::group_by(.data$Quant) %>%
              dplyr::summarise(Value = sum(.data$Value))
      )
    )

    SS <- agg %>%
      dplyr::group_by(.data$Quant) %>%
      dplyr::summarise(ss = sum(.data$Value))

    if(any(SS$ss > 105)) warning(paste("Maximum taxa sum is", round(max(SS$ss), 2),
                                       "and will be rescaled to equal 100."),
                                 call. = FALSE)
    if(any(SS$ss < 95)) warning(paste("Minimum taxa sum is", round(min(SS$ss), 2),
                                      "and will be rescaled to equal 100."),
                                call. = FALSE)

    agg %<>%
      plyr::ddply(~ .data$Quant,
                  function(set) set %<>% dplyr::mutate(Value = 100*.data$Value/sum(.data$Value))
      )
  }
  agg
}

bb_agg_ccint <- function(msum, ...){
  agg <- rbind(
    msum %>%
      dplyr::filter(.data$ord=="Top") %>%
      dplyr::select(.data$Taxa, .data$Coef, .data$bb.avg),

    cbind(Taxa = "Other",
          msum %>% dplyr::filter(.data$ord == "Other") %>%
            dplyr::group_by(.data$Coef) %>%
            dplyr::summarize(bb.avg = sum(.data$bb.avg))
    )
  ) %>%
    tidyr::separate(.data$Coef, cov_str(...), "[.]")
  SS <- agg %>%
    dplyr::group_by_(cov_str(...)[1], cov_str(...)[2]) %>%
    dplyr::summarise(ss = sum(.data$bb.avg))

  if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                     "and will be rescaled to equal 100."),
                               call. = FALSE)
  if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                    "and will be rescaled to equal 100."),
                              call. = FALSE)

  agg %<>%
    dplyr::group_by_(cov_str(...)[1], cov_str(...)[2]) %>%
    dplyr::mutate(bb.avg = 100*.data$bb.avg/sum(.data$bb.avg))

  agg
}

bb_agg_qcint <- function(msum, quant_style, ...){
  if(quant_style == "continuous"){
    oo <- msum %>% dplyr::filter(.data$ord == "Other") %>%
      dplyr::group_by(.data$Quant, .data$Fac) %>%
      dplyr::summarise(Value = sum(.data$Value))

    agg <- rbind(
      msum %>%
        dplyr::filter(.data$ord=="Top") %>%
        dplyr::select(.data$Taxa, .data$Quant, .data$Fac, .data$Value),

      cbind(Taxa = rep("Other", nrow(oo)), oo)
    )

    SS <- agg %>%
      dplyr::group_by(.data$Quant, .data$Fac) %>%
      dplyr::summarise(ss = sum(.data$Value))

    if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                       "and will be rescaled to equal 100."),
                                 call. = FALSE)
    if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                      "and will be rescaled to equal 100."),
                                call. = FALSE)

    agg %<>%
      plyr::ddply( ~ .data$Quant + .data$Fac,
                   .fun = function(set){
                     set %<>%
                       dplyr::mutate(Value = 100*.data$Value/sum(.data$Value),
                                     Fac = factor(.data$Fac))
                   }
      )
  }

  if(quant_style == "discrete"){
    oo <- msum %>% dplyr::filter(.data$ord == "Other") %>%
      dplyr::group_by(.data$Fac, .data$HL) %>%
      dplyr::summarize(Est = sum(.data$Est))

    agg <- rbind(
      msum %>%
        dplyr::filter(.data$ord=="Top") %>%
        dplyr::select(.data$Taxa, .data$Fac, .data$HL, .data$Est),

      cbind(Taxa = rep("Other", nrow(oo)), oo)
    )

    SS <- agg %>%
      dplyr::group_by(.data$Fac,.data$HL) %>%
      dplyr::summarise(ss = sum(.data$Est))

    if(sum(agg$nb.up) > 105) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                           "and will be rescaled to equal 100."),
                                     call. = FALSE)
    if(sum(agg$nb.low) < 95) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                           "and will be rescaled to equal 100."),
                                     call. = FALSE)

    agg %<>%
      dplyr::group_by(.data$Fac,.data$HL) %>%
      dplyr::mutate(Est = 100*.data$Est/sum(.data$Est))
  }
  agg
}

bb_agg_qqint <- function(msum, ...){

  oo <- msum %>% dplyr::filter(.data$ord == "Other") %>%
    dplyr::group_by(.data$Fac1, .data$Fac2) %>%
    dplyr::summarize(Est = sum(.data$Est))

  agg <- rbind(
    msum %>%
      dplyr::filter(.data$ord=="Top") %>%
      dplyr::select(.data$Taxa, .data$Fac1,
                    .data$Fac2, .data$Est),

    cbind(Taxa = rep("Other", nrow(oo)), oo)
  )

  SS <- agg %>%
    dplyr::group_by(.data$Fac1,.data$Fac2) %>%
    dplyr::summarise(ss = sum(.data$Est))

  if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                     "and will be rescaled to equal 100."),
                               call. = FALSE)
  if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                    "and will be rescaled to equal 100."),
                              call. = FALSE)

  agg %<>%
    dplyr::group_by(.data$Fac1,.data$Fac2) %>%
    dplyr::mutate(Est = 100*.data$Est/sum(.data$Est)) %>%
    dplyr::ungroup()

  agg$Fac1 <- factor(agg$Fac1, levels = c("Low", "High"))
  agg$Fac2 <- factor(agg$Fac2, levels = c("Low", "High"))

  agg
}

#### Estimation ####

bb_categ_est <- function(msum){
  msum %>%
    dplyr::mutate(bb.avg = ifelse(.data$Cov_Type == "Intercept", exp(.data$Intercept)/
                                    (1+exp(.data$Intercept)),

                                  exp(.data$Intercept + .data$Estimate)/
                                    (1+exp(.data$Intercept + .data$Estimate)) ) * 100
    )
}

bb_ccint_est <- function(msum, l1, l2){

  l1.fe <- numeric(0) ; l2.fe <- numeric(0) ; l.int <- numeric(0)

  # Switching variables if l2 has more levels so for loops work
  if(length(l2) > length(l1)){
    L1 <- l2; L2 <- l1
    l1 <- L1; l2 <- L2
  }

  for(i in 2:length(l1)){
    l1.fe[i] <- unique(msum$Intercept) +
      msum$Estimate[stringr::str_detect(msum$Coef, l1[i]) & msum$Cov_Type != "c*c.int"]
  }

  for(i in 2:length(l2)){
    l2.fe[i] <- unique(msum$Intercept) +
      msum$Estimate[stringr::str_detect(msum$Coef, l2[i]) & msum$Cov_Type != "c*c.int"]
  }

  for(i in 2:length(l2)){
    for(j in 2:length(l1)){
      l.int[j] <- l1.fe[j] +
        msum$Estimate[stringr::str_detect(msum$Coef, l2[i]) & msum$Cov_Type == "categ"] +
        msum$Estimate[stringr::str_detect(msum$Coef, l1[j]) & msum$Cov_Type == "c*c.int"]

    }
  }

  betas <- c(unique(msum$Intercept),c(l1.fe,l2.fe,l.int)[!is.na(c(l1.fe,l2.fe,l.int))])

  data.frame(bb.avg = exp(betas)/(1+exp(betas)) * 100,
    Coef = levels(interaction(l1,l2)),
    Cov_Type = msum$Cov_Type)
}

bb_qcint_est <- function(msum, m_cov, quant_style, range, ...){

  if(quant_style == "continuous"){
    if(is.numeric(m_cov[,cov_str(...)[1]])) {
      levs <- sort(unique(m_cov[,cov_str(...)[2]]))
      ll <- paste(cov_str(...)[2], levs, sep = "")
    } else {
      levs <- sort(unique(m_cov[,cov_str(...)[1]]))
      ll <- paste(cov_str(...)[1], levs, sep = "")
    }

    RR <- seq(range[1], range[2], length.out = 10)

    ## The estimate from the reference level along the range of RR
    Est <- data.frame(Quant = RR,
                      l0 = exp( unique(msum$Intercept) +
                                  msum$Estimate[msum$Cov_Type == "quant"]*RR ) /
                        (1 + exp( unique(msum$Intercept) +
                                   msum$Estimate[msum$Cov_Type == "quant"]*RR ))
    )

    for(i in 2:length(ll)){ ## starting at the second level (one below reference)

      eb <- exp(
        unique(msum$Intercept) + ## interaction intercept
          msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "categ"] +

          ## interaction slope
          (msum$Estimate[msum$Cov_Type == "quant"] +
             msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "q*c.int"])*RR)

      Est[, i+1 ] <- eb / (1+eb)
    }

    Est %<>%
      tidyr::pivot_longer(-.data$Quant, names_to = "Fac", values_to = "Value") %>%
      dplyr::mutate(Value = .data$Value * 100)
  }

  if(quant_style == "discrete"){
    if(is.numeric(m_cov[,cov_str(...)[1]])) {
      levs <- sort(unique(m_cov[,cov_str(...)[2]]))
      ll <- paste(cov_str(...)[2], levs, sep = "")
    } else {
      levs <- sort(unique(m_cov[,cov_str(...)[1]]))
      ll <- paste(cov_str(...)[1], levs, sep = "")
    }

    low <- exp( unique(msum$Intercept) +
                  msum$Estimate[msum$Cov_Type == "quant"]*range[1]) /
      (1 +  exp( unique(msum$Intercept) +
                    msum$Estimate[msum$Cov_Type == "quant"]*range[1]))

    high <- exp( unique(msum$Intercept) +
                   msum$Estimate[msum$Cov_Type == "quant"]*range[2]) /
      (1 +  exp( unique(msum$Intercept) +
                   msum$Estimate[msum$Cov_Type == "quant"]*range[2]))

    for(i in 2:length(ll)){ ## starting at the second level (one below reference)
     eb <- exp(
        unique(msum$Intercept) +
          msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "categ"] +
          (msum$Estimate[msum$Cov_Type == "quant"] +
             msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "q*c.int"]
          )*range[1] )

      low[i] <- eb / (1+eb)
    }

    for(i in 2:length(ll)){ ## starting at the second level (one below reference)

      eb <- exp(
        unique(msum$Intercept) +
          msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "categ"] +
          (msum$Estimate[msum$Cov_Type == "quant"] +
             msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "q*c.int"]
          )*range[2] )

      high[i] <- eb / (1+eb)
    }

    Est <- data.frame(
      Est = c(low,high) * 100,
      Fac = factor(levs),
      HL = rep(c("Low","high"), each = length(levs))
    )
  }

  Est
}

bb_qqint_est <- function(msum, range, ...){

  LL <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,1] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1]*range[2,1])

  LL <- LL / (1+LL)

  LH <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,2] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1]*range[2,2])

  LH <- LH / (1+LH)

  HL <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,1] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2]*range[2,1])

  HL <- HL / (1+HL)

  HH <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,2] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2]*range[2,2])

  HH <- HH / (1+HH)

  data.frame(Taxa = msum$Taxa,
             Est = c(LL,LH,HL,HH) * 100,
             Fac1 = rep(c("Low","High"), each = 2),
             Fac2 = rep(c("Low","High"), times = 2),
             nb.avg = mean( c(LL,LH,HL,HH) *100 )
  )
}

#### Plotting ####

bb_categ_bars <- function(msum, xaxis, main, xlab, ylab, subtitle, lines){

  ## Making sure that "Other" is at the bottom of the legend
  msum %<>%
    dplyr::mutate(Taxa = factor(.data$Taxa,
                         levels = c(as.character(sort(unique(.data$Taxa)[unique(.data$Taxa) != "Other"])),
                                    as.character(unique(.data$Taxa)[unique(.data$Taxa) == "Other"])))
    )

  gg <-  ggplot2::ggplot(msum, ggplot2::aes(.data$Coef,.data$bb.avg)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
    ggplot2::scale_x_discrete(labels = xaxis)

  if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
  else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))
  gg
}

bb_quant_bars <- function(msum, range, quant_style, xaxis, main,
                          xlab, ylab, subtitle, lines){

  if(quant_style == "continuous"){
    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x = as.numeric(.data$Quant),
                                   y = .data$Value, fill = .data$Taxa)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle)

    if(lines) gg <- gg + ggplot2::geom_area(colour = "black")
    else gg <- gg + ggplot2::geom_area()
  }

  if(quant_style == "discrete"){
    gg <- msum %>%
      dplyr::select(.data$Taxa, .data$bb.low, .data$bb.up) %>%
      tidyr::pivot_longer(-.data$Taxa, names_to = "varilabe", values_to = "value") %>%
      ggplot2::ggplot(ggplot2::aes(.data$variable, .data$value)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
      ggplot2::scale_x_discrete(labels = xaxis)

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))
  }
  gg
}

bb_ccint_bars <- function(msum, xaxis, main, xlab, ylab, lines,
                          subtitle, facet_labels, facet_layout, ...){
  Cov <- parse(text = cov_str(...)[1])
  Fac <- parse(text = cov_str(...)[2])

  gg <- ggplot2::ggplot(msum, ggplot2::aes(eval(Cov), .data$bb.avg)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
    ggplot2::scale_x_discrete(labels = xaxis)

  if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
  else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

  if(facet_layout == 1){
    gg + ggplot2::facet_grid(. ~ eval(Fac),
                             labeller = ggplot2::labeller(Fac = facet_labels))
  } else{
    gg + ggplot2::facet_grid(eval(Fac) ~ .,
                             labeller = ggplot2::labeller(Fac = facet_labels)) +
      ggplot2::coord_flip()
  }
}

bb_qcint_bars <- function(msum, m_cov, quant_style, xaxis, main, xlab, ylab,
                          subtitle, facet_labels, facet_layout, lines, ...){

  if(is.null(facet_labels)){
    if(is.factor(m_cov[,cov_str(...)[1]])){
      facet_labels <- levels(m_cov[,cov_str(...)[1]])
    } else facet_labels <- levels(m_cov[,cov_str(...)[2]])
  }

  names(facet_labels) <- levels(msum$Fac)

  if(quant_style == "continuous"){

    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x = as.numeric(.data$Quant),
                                   y = .data$Value, fill = .data$Taxa)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle)

    if(lines) gg <- gg + ggplot2::geom_area(colour = "black")
    else gg <- gg + ggplot2::geom_area()

    if(facet_layout == 1){
      gg <- gg + ggplot2::facet_grid(. ~ eval(Fac),
                               labeller = ggplot2::labeller(Fac = facet_labels))
    } else{
      gg <- gg + ggplot2::facet_grid(eval(Fac) ~ .,
                               labeller = ggplot2::labeller(Fac = facet_labels)) +
        ggplot2::coord_flip()
    }
  }


  if(quant_style == "discrete"){
    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$HL, y=.data$Est)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
      ggplot2::scale_x_discrete(labels = xaxis)

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

    if(facet_layout == 1){
      gg <- gg + ggplot2::facet_grid(. ~ eval(Fac),
                                     labeller = ggplot2::labeller(Fac = facet_labels))
    } else{
      gg <- gg + ggplot2::facet_grid(eval(Fac) ~ .,
                                     labeller = ggplot2::labeller(Fac = facet_labels)) +
        ggplot2::coord_flip()
    }
  }

  gg
}

bb_qqint_bars <- function(msum, xaxis, main, xlab, ylab, lines,
                          subtitle, facet_labels, facet_layout){

  if(facet_layout == 1){

    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Fac1, y=.data$Est)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
      ggplot2::scale_x_discrete(labels = xaxis) +
      ggplot2::facet_grid(. ~ Fac2, labeller = ggplot2::labeller(Fac2 = facet_labels))

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

  }

  if(facet_layout == 2){

    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Fac2, y=.data$Est)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle) +
      ggplot2::scale_x_discrete(labels = xaxis) +
      ggplot2::facet_grid(. ~ Fac1, labeller = ggplot2::labeller(Fac1 = facet_labels))

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))
  }

  gg
}

