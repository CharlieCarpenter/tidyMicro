###############################
##
## Project: tidyMicro
##
## Purpose: nb_bars helper functions
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

nb_bars_categ <- function(modsum, ..., cc, top_taxa, RA, specific_taxa,
                          xaxis, xlab, ylab, main, subtitle, lines){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef %<>% plyr::ddply(~ .data$Taxa, quant_cont, modsum$Model_Covs, ...)
  }

  msum <- modsum$Model_Coef %>%
    dplyr::filter(.data$Cov_Type == "Intercept" |
                    stringr::str_detect(.data$Coef, cov_str(...))) %>%
    plyr::ddply(~ .data$Taxa, categ_est) %>%
    Taxa_ord() %>%
    dplyr::mutate(Coef = factor(.data$Coef, levels = unique(.data$Coef)))

  ord <- msum %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarize(mean.nb = mean(.data$nb.avg)) %>%
    dplyr::arrange(dplyr::desc(.data$mean.nb))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style = NULL, ...) %>%
      categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab, subtitle=subtitle, lines=lines)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style = NULL, ...) %>%
      categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab, subtitle=subtitle, lines = lines)
  } else {
    msum %>%
      categ_bars(xaxis=xaxis, main=main, xlab=xlab, ylab=ylab, subtitle=subtitle, lines = lines)
  }
}

nb_bars_quant <- function(modsum, ..., cc, range, quant_style, top_taxa, RA,
                          specific_taxa, main, xaxis, xlab, ylab, subtitle, lines){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef %<>% plyr::ddply(~ .data$Taxa, quant_cont, modsum$Model_Covs,...)
  }

  if(quant_style == "continuous"){
    RR <- seq(range[1], range[2], length.out = 10)

    msum <- modsum$Model_Coef %>%
      dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...)))

    ## Estimates for taxa along RR range
    Est <- exp(tcrossprod(msum$Estimate, RR) + msum$Intercept + log(100))
    colnames(Est) <- RR

    msum %<>%
      cbind(Est) %>%
      dplyr::select(.data$Taxa, 6:15) %>%
      tidyr::pivot_longer(-.data$Taxa, names_to = "Quant", values_to = "Value") %>%
      Taxa_ord()

    ord <- msum %>%
      dplyr::group_by(.data$Taxa) %>%
      dplyr::summarise(nb.avg = mean(.data$Value)) %>%
      dplyr::arrange(dplyr::desc(.data$nb.avg))
  }

  if(quant_style == "discrete"){
    msum <- modsum$Model_Coef %>%
      dplyr::filter(stringr::str_detect(.data$Coef, cov_str(...))) %>%
      dplyr::mutate(nb.avg = exp(.data$Intercept + .data$Estimate*mean(.data$range) + log(100)),
                    nb.low = exp(.data$Intercept + .data$Estimate*.data$range[1] + log(100)),
                    nb.up = exp(.data$Intercept + .data$Estimate*.data$range[2] + log(100))
      )

    ord <- msum %>%
      dplyr::arrange(dplyr::desc(.data$nb.avg))
  }

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style, ...) %>%
      quant_bars(range = range, quant_style = quant_style, xaxis = xaxis, main = main,
                 xlab = xlab, ylab = ylab, subtitle=subtitle, lines = lines)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style, ...) %>%
      quant_bars(range = range, quant_style = quant_style, xaxis = xaxis, main = main,
                 xlab = xlab, ylab = ylab, subtitle=subtitle, lines = lines)
  } else
    msum %>%
    quant_bars(range = range, quant_style = quant_style, xaxis = xaxis, main = main,
               xlab = xlab, ylab = ylab, subtitle=subtitle, lines = lines)

}

nb_bars_ccint <- function(modsum, ..., cc, top_taxa, RA, specific_taxa, xaxis, lines,
                          xlab, ylab, main, subtitle, facet_labels, facet_layout){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs, ...)
  }

  l1 <- levels(modsum$Model_Covs[, cov_str(...)[1]]) ; l2 <- levels(modsum$Model_Covs[, cov_str(...)[2]])
  if( sum((l1%in%l2)) > 0) stop("Factor variables can not share levels.")


  msum <- suppressWarnings(
    modsum$Model_Coef %>%
      dplyr::filter(.data$Cov_Type == "Intercept" |
                      stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                      stringr::str_detect(.data$Coef, cov_str(...)[2]) ) %>%
      Taxa_ord() %>%
      plyr::ddply(~ .data$Taxa, ccint_est, l1, l2, ...)
  ) %>% dplyr::rename(Taxa = .data$`.data$Taxa`)

  ord <- msum %>%
    dplyr::group_by(.data$Taxa) %>%
    dplyr::summarise(nb.avg = mean(.data$nb.avg)) %>%
    dplyr::arrange(dplyr::desc(.data$nb.avg))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style = NULL, ...) %>%
      ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                 lines = lines, facet_labels = facet_labels, facet_layout = facet_layout, ...)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style = NULL, ...) %>%
      ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                 lines = lines, facet_labels = facet_labels, facet_layout = facet_layout, ...)
  } else
    msum %>%
    ccint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
               lines = lines, facet_labels = facet_labels, facet_layout = facet_layout, ...)
}

nb_bars_qcint <- function(modsum, ..., cc, range, quant_style, top_taxa, RA, specific_taxa, xaxis,
                          xlab, ylab, main, subtitle, facet_labels, facet_layout, lines){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs, ...)
  }

  msum <- modsum$Model_Coef %>%
    dplyr::filter(.data$Coef == "(Intercept)" |
                    stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                    stringr::str_detect(.data$Coef, cov_str(...)[2])) %>%
    Taxa_ord() %>%
    plyr::ddply(~ .data$Taxa, qcint_est, modsum$Model_Covs,
                quant_style, range, ...) %>%
    dplyr::rename(Taxa = .data$`.data$Taxa`)

  if(quant_style == "continuous"){
    ord <- msum %>%
      dplyr::group_by(.data$Taxa) %>%
      dplyr::summarise(nb.avg = mean(.data$Value)) %>%
      dplyr::arrange(dplyr::desc(.data$nb.avg))
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

    agg_bars(msum, cc, quant_style, ...) %>%
      qcint_bars(m_cov = modsum$Model_Covs,
                 quant_style = quant_style, xaxis = xaxis, main = main, xlab = xlab,
                 ylab = ylab, subtitle = subtitle, lines = lines,
                 facet_labels = facet_labels, facet_layout = facet_layout, ...)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style, ...) %>%
      qcint_bars(m_cov = modsum$Model_Covs,
                 quant_style = quant_style, xaxis = xaxis, main = main, xlab = xlab,
                 ylab = ylab, subtitle = subtitle, lines = lines,
                 facet_labels = facet_labels, facet_layout = facet_layout, !!!rlang::quos(...))
  } else
    msum %>%
    qcint_bars(m_cov = modsum$Model_Covs,
               quant_style = quant_style, xaxis = xaxis, main = main, xlab = xlab,
               ylab = ylab, subtitle = subtitle, lines = lines,
               facet_labels = facet_labels, facet_layout = facet_layout, !!!rlang::quos(...))
}

nb_bars_qqint <- function(modsum, ..., cc, range, top_taxa, RA, specific_taxa, xaxis, lines,
                          xlab, ylab, main, subtitle, facet_labels, facet_layout){

  if("quant" %in% modsum$Model_Coef$Cov_Type){
    modsum$Model_Coef <- plyr::ddply(modsum$Model_Coef,
                                     ~ .data$Taxa, quant_cont, modsum$Model_Covs,...)
  }

  msum <- modsum$Model_Coef %>%
    dplyr::filter(.data$Coef == "(Intercept)" |
                    stringr::str_detect(.data$Coef, cov_str(...)[1]) |
                    stringr::str_detect(.data$Coef, cov_str(...)[2])) %>%
    Taxa_ord() %>%
    plyr::ddply(~ .data$Taxa, qqint_est, range, ...) %>%
    dplyr::rename(Taxa = .data$`.data$Taxa`)


  ord <- msum %>%
    dplyr::arrange(dplyr::desc(.data$nb.avg))

  if(top_taxa > 0){
    if("Other" %in% ord$Taxa[seq(1,top_taxa)]) top_taxa <- top_taxa + 1

    msum$ord <- ifelse(msum$Taxa %in% unique(ord$Taxa)[seq(1,top_taxa)] & msum$Taxa != "Other", "Top", "Other")
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style = NULL, ...) %>%
      qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                 lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)
  } else if(RA > 0 & RA < 100){
    msum$ord <- ifelse(msum$Taxa == "Other", "Other",
                       ifelse(msum$nb.avg >= RA, "Top", "Other"))
    msum$ord <- ifelse(msum$Taxa %in% specific_taxa, "Top", msum$ord)

    agg_bars(msum, cc, quant_style = NULL, ...) %>%
      qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
                 lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)
  } else
    msum %>%
    qqint_bars(xaxis = xaxis, main = main, xlab = xlab, ylab = ylab, subtitle=subtitle,
               lines = lines, facet_labels = facet_labels, facet_layout = facet_layout)
}

categ_bars <- function(msum, xaxis, main, xlab, ylab, subtitle, lines){

  gg <- msum %>%
    ggplot2::ggplot(ggplot2::aes(.data$Coef,.data$nb.avg)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa") +
    ggplot2::scale_x_discrete(labels = xaxis)

  if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
  else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))
  gg
}

quant_bars <- function(msum, range, quant_style, xaxis, main,
                       xlab, ylab, subtitle, lines){

  if(quant_style == "continuous"){
    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x = as.numeric(.data$Quant),
                                   y = .data$Value, fill = .data$Taxa)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa")

    if(lines) gg <- gg + ggplot2::geom_area(colour = "black")
    else gg <- gg + ggplot2::geom_area()
  } else if(quant_style == "discrete"){
    gg <- msum %>%
      dplyr::select(.data$Taxa, .data$nb.low, .data$nb.up) %>%
      tidyr::pivot_longer(-.data$Taxa, names_to = "variable", values_to = "value") %>%
      ggplot2::ggplot(ggplot2::aes(.data$variable,.data$value)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa") +
      ggplot2::scale_x_discrete(labels = xaxis)

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))
  }
  gg
}

ccint_bars <- function(msum, xaxis, main, xlab, ylab, lines,
                       subtitle, facet_labels, facet_layout, ...){

  Cov <- parse(text = cov_str(...)[1])
  Fac <- parse(text = cov_str(...)[2])

  gg <- msum %>%
    ggplot2::ggplot(ggplot2::aes(eval(Cov), .data$nb.avg)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa") +
    ggplot2::scale_x_discrete(labels = xaxis)

  if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
  else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

  if(!is.null(facet_labels) & is.null(names(facet_labels))){
    names(facet_labels) <- unique(msum$Fac)
  }

  if(facet_layout == 1){
    gg + ggplot2::facet_grid(. ~ eval(Fac),
                             labeller = ggplot2::labeller(Fac = facet_labels))
  } else {
    gg + ggplot2::facet_grid(eval(Fac) ~ .,
                             labeller = ggplot2::labeller(Fac = facet_labels)) +
      ggplot2::coord_flip()
  }
}

qcint_bars <- function(msum, m_cov, quant_style, xaxis, main, xlab, ylab,
                       subtitle, facet_labels, facet_layout, lines, ...){

  if(is.null(facet_labels)){
    if(is.numeric(m_cov[,cov_str(...)[1]])){
      facet_labels <- sort(unique(m_cov[,cov_str(...)[2]]))
    } else facet_labels <- sort(unique(m_cov[,cov_str(...)[1]]))
  }

  names(facet_labels) <- levels(msum$Fac)

  if(quant_style == "continuous"){
    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x = as.numeric(.data$Quant),
                          y = .data$Value, fill = .data$Taxa)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa")

    if(lines) gg <- gg + ggplot2::geom_area(colour = "black")
    else gg <- gg + ggplot2::geom_area()

    if(!is.null(facet_labels) & is.null(names(facet_labels))){
      names(facet_labels) <- unique(msum$Fac)
    }

    if(facet_layout == 1){
      gg <- gg + ggplot2::facet_grid(. ~ Fac,
                                     labeller = ggplot2::labeller(Fac = facet_labels))
    } else {
      gg <- gg + ggplot2::facet_grid(Fac ~ .,
                                     labeller = ggplot2::labeller(Fac = facet_labels)) +
        ggplot2::coord_flip()
    }
  }

  if(quant_style == "discrete"){
    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$HL, y=.data$Est)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa") +
      ggplot2::scale_x_discrete(labels = xaxis)

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

    if(facet_layout == 1){
      gg <- gg + ggplot2::facet_grid(. ~ Fac,
                                     labeller = ggplot2::labeller(Fac = facet_labels))
    }
    else{
      gg <- gg + ggplot2::facet_grid(Fac ~ .,
                                     labeller = ggplot2::labeller(Fac = facet_labels)) +
        ggplot2::coord_flip()
    }
  }

  gg
}

qqint_bars <- function(msum, xaxis, main, xlab, ylab, lines,
                       subtitle, facet_labels, facet_layout){

  if(!is.null(facet_labels) & is.null(names(facet_labels))){
    names(facet_labels) <- unique(msum$Fac)
  }

  if(facet_layout == 1){

    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Fac1, y=.data$Est)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa") +
      ggplot2::scale_x_discrete(labels = xaxis) +
      ggplot2::facet_grid(. ~ Fac2, labeller = ggplot2::labeller(Fac2 = facet_labels))

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))

  }

  if(facet_layout == 2){

    gg <- msum %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$Fac2, y=.data$Est)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main, x = xlab, y = ylab, subtitle = subtitle, fill = "Taxa") +
      ggplot2::scale_x_discrete(labels = xaxis) +
      ggplot2::facet_grid(. ~ Fac1, labeller = ggplot2::labeller(Fac1 = facet_labels))

    if(lines) gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa), colour = "black")
    else gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$Taxa))
  }

  gg
}

by_fun <- function(...){
  rlang::quos(...) %>%
    rlang::splice() %>%
    unlist %>%
    as.character %>%
    stringr::str_replace_all(pattern = "~", replacement = "") %>%
    return()
}

nb_type <- function(modsum, ...){
  Cov <- rlang::quos(...) %>%
    rlang::splice() %>%
    unlist %>%
    as.character %>%
    stringr::str_remove("~")  ## Makes '...' into a character string

  if(length(Cov) < 1) stop("nb_Bars requires a model covariate. If looking for stacked bar charts of raw counts use the function ra_bars.")
  if(length(Cov) > 1) stop("nb_Bars requires the use of only one model covariate. That covariate can be an interaction term such as 'Group*Age'")

  cofs <- modsum$Model_Coef[modsum$Model_Coef$Taxa == unique(modsum$Model_Coef$Taxa)[1],]

  if(length(cov_str(...)) == 1){
    cc <- cofs$Cov_Type[stringr::str_detect(cofs$Coef, cov_str(...)) &
                          !stringr::str_detect(cofs$Coef, ":")] %>%
      unique
  }

  if(length(cov_str(...)) == 2){
    cc <- cofs$Cov_Type[stringr::str_detect(cofs$Coef, cov_str(...)[1]) &
                          stringr::str_detect(cofs$Coef, ":")] %>%
      unique
  }

  cc
}

agg_bars <- function(msum, cc, quant_style, ...){
  if(cc == "categ") {
    agg <- nb_agg_categ(msum, ...)
  } else if(cc == "quant"){
    agg <- nb_agg_quant(msum, quant_style, ...)
  } else if(cc == "c*c.int"){
    agg <- nb_agg_ccint(msum, ...)
  } else if(cc == "q*c.int"){
    agg <- nb_agg_qcint(msum, quant_style, ...)
  } else if(cc == "q*q.int"){
    agg <- nb_agg_qqint(msum, ...)
  }
  agg
}

nb_agg_categ <- function(msum, ...){
  agg <- rbind(
    msum %>%
      dplyr::filter(.data$ord=="Top") %>%
      dplyr::select(.data$Taxa, .data$Coef, .data$nb.avg),

    cbind(Taxa = "Other",
          msum %>% dplyr::filter(.data$ord == "Other") %>%
            dplyr::group_by(.data$Coef) %>%
            dplyr::summarize(nb.avg = sum(.data$nb.avg))
    )
  )
  SS <- agg %>% dplyr::group_by(.data$Coef) %>%
    dplyr::summarize(ss = sum(.data$nb.avg))

  if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                     "and will be rescaled to equal 100."),
                               call. = FALSE)
  if(any(SS$ss < 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                    "and will be rescaled to equal 100."),
                              call. = FALSE)

  agg %<>% plyr::ddply(~ .data$Coef,
                       function(set){
                         set %<>% dplyr::mutate(nb.avg = 100*.data$nb.avg/sum(.data$nb.avg))
                       })
  agg
}

nb_agg_quant <- function(msum, quant_style, ...){
  if(quant_style == "discrete"){
    agg <- rbind(
      msum %>%
        dplyr::filter(.data$ord=="Top") %>%
        dplyr::select(.data$Taxa, .data$Coef,
                      .data$nb.low, .data$nb.up),

      cbind(Taxa = "Other",
            msum %>% dplyr::filter(.data$ord == "Other") %>%
              dplyr::group_by(.data$Coef) %>%
              dplyr::summarize(nb.low = sum(.data$nb.low), nb.up = sum(.data$nb.up))
      )
    )

    if(sum(agg$nb.up) > 105) warning(paste("Maximum column sum is", round(sum(agg$nb.up)),
                                           "and will be rescaled to equal 100."),
                                     call. = FALSE)
    if(sum(agg$nb.low) < 95) warning(paste("Minimum column sum is", round(sum(agg$nb.low)),
                                           "and will be rescaled to equal 100."),
                                     call. = FALSE)

    agg %<>%
      dplyr::mutate(nb.low = 100*.data$nb.low/sum(.data$nb.low),
                    nb.up = 100*.data$nb.up/sum(.data$nb.up))
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

nb_agg_ccint <- function(msum, ...){
  agg <- rbind(
    msum %>%
      dplyr::filter(.data$ord=="Top") %>%
      dplyr::select(.data$Taxa, .data$Coef, .data$nb.avg),

    cbind(Taxa = "Other",
          msum %>% dplyr::filter(.data$ord == "Other") %>%
            dplyr::group_by(.data$Coef) %>%
            dplyr::summarize(nb.avg = sum(.data$nb.avg))
    )
  ) %>%
    tidyr::separate(.data$Coef, cov_str(...), "[.]")
  SS <- agg %>%
    dplyr::group_by_(cov_str(...)[1], cov_str(...)[2]) %>%
    dplyr::summarise(ss = sum(.data$nb.avg))

  if(any(SS$ss > 105)) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                         "and will be rescaled to equal 100."),
                                   call. = FALSE)
  if(any(SS$ss > 95)) warning(paste("Minimum column sum is", round(min(SS$ss)),
                                         "and will be rescaled to equal 100."),
                                   call. = FALSE)

  agg %<>%
    dplyr::group_by_(cov_str(...)[1], cov_str(...)[2]) %>%
    dplyr::mutate(nb.avg = 100*.data$nb.avg/sum(.data$nb.avg))

  agg
}

nb_agg_qcint <- function(msum, quant_style, ...){
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
      plyr::ddply(~ .data$Quant + .data$Fac,
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

      cbind(Taxa = "Other", oo)
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

nb_agg_qqint <- function(msum, ...){

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

  if(sum(agg$nb.up) > 105) warning(paste("Maximum column sum is", round(max(SS$ss)),
                                         "and will be rescaled to equal 100."),
                                   call. = FALSE)
  if(sum(agg$nb.low) < 95) warning(paste("Minimum column sum is", round(min(SS$ss)),
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

quant_cont <- function(m_coef, m_cov, ...){
  q <- m_cov[ ,m_coef$Coef[m_coef$Cov_Type == "quant" &
                             !(m_coef$Coef %in% cov_str(...))] %>% as.character]

  if(is.data.frame(q)){
    qm <- apply(q, 2, mean, na.rm = T) %>% as.vector
  } else qm <- mean(q, na.rm = TRUE) %>% as.vector


  m_coef$Intercept <- unique(m_coef$Intercept) +
    sum(m_coef$Estimate[m_coef$Cov_Type == "quant" &
                        !(m_coef$Coef %in% cov_str(...))]*qm)

  m_coef
}

categ_est <- function(msum){
  msum %>%
    dplyr::mutate(nb.avg = ifelse(.data$Cov_Type == "Intercept", exp(.data$Intercept + log(100)),
                                  exp(.data$Intercept + .data$Estimate + log(100)))
    )
}

ccint_est <- function(msum, l1, l2, ...){

  l1.fe <- numeric(0) ; l2.fe <- numeric(0) ; l.int <- numeric(0)

  for(i in seq(2,length(l1))){
    l1.fe[i] <- unique(msum$Intercept) +
      msum$Estimate[stringr::str_detect(msum$Coef, l1[i]) & msum$Cov_Type != "c*c.int"]
  }

  for(i in seq(2,length(l2))){
    l2.fe[i] <- unique(msum$Intercept) +
      msum$Estimate[stringr::str_detect(msum$Coef, l2[i]) & msum$Cov_Type != "c*c.int"]
  }

  if(length(l2) <= length(l1)){
    for(i in seq(2,length(l2))){
      for(j in seq(2,length(l1))){
        l.int[j] <- l1.fe[j] +
          msum$Estimate[stringr::str_detect(msum$Coef, l2[i]) & msum$Cov_Type == "categ"] +
          msum$Estimate[stringr::str_detect(msum$Coef, l1[j]) & msum$Cov_Type == "c*c.int"]

      }
    }
  }

  if(length(l1) < length(l2)){
    for(i in seq(2,length(l1))){
      for(j in seq(2,length(l2))){
        l.int[j] <- l2.fe[j] +
          msum$Estimate[stringr::str_detect(msum$Coef, l1[i]) & msum$Cov_Type == "categ"] +
          msum$Estimate[stringr::str_detect(msum$Coef, l2[j]) & msum$Cov_Type == "c*c.int"]

      }
    }
  }

  data.frame(nb.avg = exp(
    c(unique(msum$Intercept), c(l1.fe,l2.fe,l.int)[!is.na(c(l1.fe,l2.fe,l.int))]) + log(100)),
    Coef = levels(interaction(l1,l2)),
    Cov_Type = msum$Cov_Type)

}

qcint_est <- function(msum, m_cov, quant_style, range, ...){

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
                                  msum$Estimate[msum$Cov_Type == "quant"]*RR + log(100) )
    )

    for(i in seq(2,length(ll))){ ## starting at the second level (one below reference)
      Est[, i+1 ] <- exp(
        unique(msum$Intercept) + ## interaction intercept
          msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "categ"] +

          ## interaction slope
          (msum$Estimate[msum$Cov_Type == "quant"] +
             msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "q*c.int"])*RR
        + log(100))
    }

    Est %<>%
      tidyr::pivot_longer(-.data$Quant, names_to = "Fac", values_to = "Value")
  }

  if(quant_style == "discrete"){
    if(is.numeric(m_cov[,cov_str(...)[1]])) {
      levs <- sort(unique(m_cov[,cov_str(...)[2]]))
      ll <- paste(cov_str(...)[2], levs, sep = "")
    } else {
      levs <- sort(unique(m_cov[,cov_str(...)[1]]))
      ll <- paste(cov_str(...)[1], levs, sep = "")
    }

    low <- high <- numeric(length(ll))
    low[1] <- exp( unique(msum$Intercept) +
                  msum$Estimate[msum$Cov_Type == "quant"]*range[1] + log(100) )
    high[1] <- exp( unique(msum$Intercept) +
                   msum$Estimate[msum$Cov_Type == "quant"]*range[2] + log(100) )

    for(i in seq(2,length(ll))){ ## starting at the second level (one below reference)
      low[i] <- exp(
        unique(msum$Intercept) +
          msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "categ"] +
          (msum$Estimate[msum$Cov_Type == "quant"] +
             msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "q*c.int"]
          )*range[1] + log(100)
      )
    }

    for(i in seq(2,length(ll))){ ## starting at the second level (one below reference)
      high[i] <- exp(
        unique(msum$Intercept) +
          msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "categ"] +
          (msum$Estimate[msum$Cov_Type == "quant"] +
             msum$Estimate[stringr::str_detect(msum$Coef, ll[i]) & msum$Cov_Type == "q*c.int"]
          )*range[2] + log(100)
      )
    }

    Est <- data.frame(
      Est = c(low,high),
      Fac = factor(levs),
      HL = rep(c("Low","high"), each = length(levs))
    )
  }

  Est
}

qqint_est <- function(msum, range, ...){

  LL <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,1] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1]*range[2,1] +
               log(100))

  LH <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,2] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,1]*range[2,2] +
               log(100))

  HL <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,1] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2]*range[2,1] +
               log(100))

  HH <- exp( unique(msum$Intercept) +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               !stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2] +
               msum$Estimate[!stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[2,2] +
               msum$Estimate[stringr::str_detect(msum$Coef, cov_str(...)[1]) &
                               stringr::str_detect(msum$Coef, cov_str(...)[2])]*range[1,2]*range[2,2] +
               log(100))

  data.frame(Taxa = msum$Taxa,
             Est = c(LL,LH,HL,HH),
             Fac1 = rep(c("Low","High"), each = 2),
             Fac2 = rep(c("Low","High"), times = 2),
             nb.avg = mean( c(LL,LH,HL,HH) + log(100) )
  )
}

## Function to reorder Taxa to make sure "Other" is at the bottom of bar plots
Taxa_ord <- function(msum){

  if("Other" %in% unique(msum$Taxa)){
    msum$Taxa <- factor(msum$Taxa,
                        levels = c(as.character(sort(unique(msum$Taxa)[unique(msum$Taxa) != "Other"])),
                                   "Other")
    )
  } else{
    msum$Taxa <- factor(msum$Taxa,
                        levels = as.character(sort(unique(msum$Taxa)))
    )
  }

  msum
}
