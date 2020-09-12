###############################
##
## Project: tidyMicro
##
## Purpose: Helper functions for bb_mods
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2019-10-10
##
## ---------------------------
## Notes:
##
##
## ---------------------------
## HEADER ----

## Function to fit VGAM:: beta binomial models and provide summary output
FE_bb <- function(set, f., trace, method., SS_type., ...){
  # printing taxa to help keep track
  message("Fitting ", as.character(unique(set$Taxa)))

  m <- try(VGAM::vglm(stats::as.formula(f.), family = "betabinomial", trace = trace, data = set),
           silent = T)

  if(!inherits(m, "try-error")){

    ## Stops if confint doesn't work. Some models converge but confint breaks
    CI <- try(VGAM::confintvglm(m, method = method.), silent = T)

    lrt <- try(VGAM::lrt.stat.vlm(m, all.out = T, omit1s = F), silent = T)

    if(!inherits(CI, "try-error") & !any(is.na(CI)) & !inherits(lrt, "try-error")){

      s <- summary(m)
      conv <- data.frame(
        Taxa = as.character(unique(set$Taxa)),                        ## Taxa Names
        Coef = names(m@coefficients),                  ## Coefficient Names
        OR = bb_coef_est(m, ...) %>%
          dplyr::pull(.data$Est) %>% exp %>% round(4),       ## Odds Ratios
        CI_95 = bb_ci_est(m, ...),                     ## Wald intervals
        LRT = lrt$Lrt.stat2 %>% round(4),              ## LRT statistic
        P_val = lrt$pvalues,                           ## P_value (not rounded for FDR calculation)
        Intercept = m@coefficients[1],                 ## Intercepts for calculations
        Beta = m@coefficients,                         ## Raw Coefs for later use
        Cov_Type = cov_type(names(m@coefficients),
                            ..., RE = NULL),  ## Covariate types for Bar Charts (no random effects)
        ## CI for individual Beta's (profile likelihood)
        CI = paste0("(", round(CI[,1],4), ", ", round(CI[,2],4), ")"),
        FE_Converged = TRUE
      ) %>% Anova_SS_bb(m, SS_type.) ## Anova for each covariate

      ## Making empty rows when taxa models break
    } else{
      conv <- data.frame(
        Taxa = as.character(unique(set$Taxa)),
        Coef = character(1), OR = NA, CI_95 = character(1), LRT = NA, P_val = NA, Intercept = NA,
        Beta = NA, Cov_Type = character(1), CI = character(1), FE_Converged = FALSE, LRT = NA)
    }
  } else {
    conv <- data.frame(
      Taxa = as.character(unique(set$Taxa)),
      Coef = character(1), OR = NA, CI_95 = character(1), LRT = NA, P_val = NA, Intercept = NA,
      Beta = NA, Cov_Type = character(1), CI = character(1), FE_Converged = FALSE, LRT = NA)
  }

}

## Can change betabinomial to betabinomial(lrho = "rhobitlink") for rho close to 0.
## Might add extra step to fit using this model lrho = "rhobitlink" if original fails
## lrho = "rhobitlink" makes it an S4 class
## Names of linear predictors: logitlink(mu), logitlink(rho): The first intercept is the logit link of the mean

bb_coef_est <- function(m,...){

  covs <- data.frame(Coef = names(m@coefficients),
                     Est = m@coefficients,
                     Cov_Type = m@coefficients %>%
                       names %>%
                       cov_type(...))

  for(i in seq(1,nrow(covs))){
    ## detecting interactions
    if( stringr::str_detect(covs$Cov_Type[i], ".int") ){
      ## interaction pieces
      int. <- covs$Coef[i] %>% as.character %>% stringr::str_split(":") %>% unlist

      ## Exact tring matches of the main effects
      ## (there will only be 2 for any 2x2 interaction piece)
      me1 <- covs$Est[stringr::str_detect(covs$Coef, paste0("^", int.[1], "$"))]
      me2 <- covs$Est[stringr::str_detect(covs$Coef, paste0("^", int.[2], "$"))]

      ## interaction detected from if statement
      covs$Est[i] <- covs$Est[i] + me1 + me2
    }
  }

  covs
}


## Estimates profile likelihood CIs using confint for main effects
## and Wald CIs for interactions after adjusting estiamtes using coef_est
bb_ci_est <- function(m, ...){

  # Updating estimates in case there are interactions
  covs <- bb_coef_est(m = m, ...)

  CI <- matrix(rep(0, nrow(covs)*2), nrow = nrow(covs), ncol = 2,
               dimnames = list(rownames(covs), c("2.5 %", "97.5 %")))

  # covariance matrix
  s <- VGAM::vcov.vlm(m)

  # Confidence intervals
  for(r in seq(1,nrow(covs))){
    ## detecting interactions
    if( stringr::str_detect(covs$Cov_Type[r], ".int") ){

      # interaction pieces
      int. <- c(covs$Coef[r] %>% as.character %>%
                  stringr::str_split(":") %>% unlist,
                covs$Coef[r] %>% as.character)

      # pullin out rows and cols with interaction and main effect covariances
      cl <- colnames(s) %in% int. ; ro <- rownames(s) %in% int. ; s_ <- s[ro,cl]

      # summing vars and covars
      v <- 0
      for(i in seq(1,nrow(s_))){
        for(j in seq(1,ncol(s_))){

          if(j >= i){
            p <- ifelse(i == j, s_[i,j], 2*s_[i,j]) # Var or 2*Cov

            v <- v + p
          }
        }
      }

      # Wald interval for interaction
      ci <- covs$Est[r] + c(-1,1)*1.96*sqrt(v)

    } else{

      # Wald interval for main effects
      ci <- m@coefficients[r] + c(-1,1)*1.96*sqrt(diag(s))[r]
    }

    CI[r,] <- ci
  }

  ## Exponentiate and round
  CI %<>% exp %>% round(4)

  ## Make neat for the table
  paste0("(", CI[,1], ", ", CI[,2], ")")

}

## Function to run multivariable LRT and attach to output df
Anova_SS_bb <- function(conv, bb, SS_type){

  ## Coefs and LRT p-value
  aa <- VGAM::anova.vglm(bb, type = SS_type)
  a <- data.frame(coef = rownames(aa), Anova = aa$`Pr(>Chi)`)

  ## Adding LRT p-values from anova.vglm
  conv$LRT <- NA
  for(i in seq(1,nrow(a))){
    conv$LRT[grepl(a$coef[i], conv$Coef)] <- a$Anova[i]
  }

  conv
}

## Makes appropriate betabinomial formula from input
formula_fun_bb <- function(...){
  covs <- rlang::quos(...) %>%           ## Making "..." a quosure
    rlang::splice() %>%            ## Splicing quosure
    unlist %>%            ## Unlisting
    stringr::str_flatten(collapse = " + ") %>%            ## Combining elements of list with "+"
    stringr::str_replace_all(pattern = "~", replacement = "")   ## Removing "~" (why is it there?)

  paste("cbind(cts, Total - cts) ~ ", covs) ## pasting with "cts ~ " for final formula
}

