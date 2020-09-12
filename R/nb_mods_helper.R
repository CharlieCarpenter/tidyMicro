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


## Function for summarizing convergent models
FE.con <- function(micro_set, f, SS_type, ...){
  ## Defining formula for glm.nb with  w/ or w/o and Offset. Offset is the default

  ## Summary Stats of convergent models
  conv <- micro_set %>%
    plyr::ddply( ~ .data$Taxa, .fun = function(set, .f, .SS_type, ...){
      m <- try(MASS::glm.nb(stats::as.formula(.f), data = set), silent = T)

      if(!inherits(m, "try-error")){
        ## Stops if confint doesn't work. Some models converge but confint breaks
        CI <- try(stats::confint(m), silent = T)
        AA <- try(car::Anova(m, SS_type, test = "LR"), silent = T)
        if(!inherits(CI, "try-error") & !any(is.na(CI)) & !inherits(AA, "try-error")){
          data.frame(
            Taxa = unique(set$Taxa),                        ## Taxa Names
            Coef = names(m$coefficients),                  ## Coefficient Names
            RR = coef_est(m, ..., RE = NULL) %>%
              dplyr::pull(.data$Est) %>% exp %>% round(4),      ## Risk Ratios
            CI_95 = ci_est(m, ..., RE = NULL),             ## Wald intervals
            Z = stats::coef(summary(m))[,3] %>% round(4),         ## Z statistic
            P_val = stats::coef(summary(m))[,4],                  ## P_value (not rounded for FDR calculation)
            Intercept = stats::coef(summary(m))[1,1],             ## Intercepts for calculations
            Beta = m$coefficients,                         ## Raw Coefs for later use
            Cov_Type = cov_type(names(m$coefficients),
                                ..., RE = NULL),  ## Covariate types for Bar Charts (no random effects)
            ## CI for individual Beta's (profile likelihood)
            CI = paste0("(", round(CI[,1],4), ", ", round(CI[,2],4), ")"),
            FE_Converged = TRUE
          ) %>% Anova_SS(m, SS_type, test = "LR") ## LRT for each covariate
        } else{
          data.frame(
            Taxa = unique(set$Taxa),
            Coef = character(1), RR = NA, CI_95 = character(1), Z = NA, P_val = NA, Intercept = NA,
            Beta = NA, Cov_Type = character(1), CI = character(1), FE_Converged = FALSE, LRT = NA)
        }
      } else {
        data.frame(
          Taxa = unique(set$Taxa),
          Coef = character(1), RR = NA, CI_95 = character(1), Z = NA, P_val = NA, Intercept = NA,
          Beta = NA, Cov_Type = character(1), CI = character(1), FE_Converged = FALSE, LRT = NA)
      }
    }, .f = f, .SS_type = SS_type, ...)

  if(nrow(conv) == 0) stop("No taxa models converged.")

  ## False Discovery Rate adjustment
  conv %>% dplyr::mutate(FDR_Pval = stats::p.adjust(.data$P_val, method = "BH") %>% round(4))
}

# RE.con <- function(micro_set, f, SS_type, ..., RE){
#
#   ## Summary Stats of convergent models
#   conv <- plyr::ldply(unique(micro_set$Taxa), function(x){
#     m <- try(lme4::glmer.nb(stats::as.formula(f), data = micro_set %>%  ## The model
#                               dplyr::filter(Taxa == as.character(x))), silent = T)
#
#     if(!inherits(m, "try-error")){
#       CI <- try(stats::confint(m), silent = T)
#       AA <- try(car::Anova(m), silent = T)
#       if(!inherits(CI, "try-error") & !any(is.na(CI)) & !inherits(AA, "try-error")){
#         data.frame(
#           Taxa = as.character(x),                             ## Taxa Names
#           Coef = coef(summary(m)) %>% rownames,               ## Coefficient Names
#           RR = coef_est(m, !!!rlang::quos(...), RE = RE) %>%
#             .$Est %>% exp %>% round(4),    ## Risk Ratios
#           CI_95 = ci_est(m, !!!rlang::quos(...), RE = RE),           ## Upper profile likelihood conidence limits
#           Z = coef(summary(m))[,3] %>% round(4),              ## Z statistic
#           P_val = coef(summary(m))[,4],                       ## P_value (not rounded for FDR calculation)
#           Intercept = m@beta[1],                  ## Intercepts for calculations
#           Beta = m@beta,                    ## Raw Coefs for later use
#           Cov_Type = cov_type(names(coef(summary(m))[,1]),
#                               !!!rlang::quos(...), RE = RE),  ## Covariate types for Bar Charts
#           ## CI for betas
#           CI = paste("(", CI[-1,1] %>% round(4),", ",
#                      CI[-1,2] %>% round(4), ")", sep = ""),
#           RE_Converged = TRUE
#         ) %>% Anova_SS(m, SS_type, test = "Chisq" ) ## Anova for each covariate
#       } else{
#         data.frame(
#           Taxa = as.character(x),
#           Coef = character(1), RR = NA, CI_95 = character(1), Z = NA, P_val = NA, Intercept = NA,
#           Beta = NA, Cov_Type = character(1), CI = character(1), FE_Converged = FALSE, LRT = NA)
#       }} else{
#         data.frame(
#           Taxa = as.character(x),
#           Coef = character(1), RR = NA, CI_95 = character(1), Z = NA, P_val = NA, Intercept = NA,
#           Beta = NA, Cov_Type = character(1), CI = character(1), FE_Converged = FALSE, LRT = NA)
#       }
#   })
#
#   if(nrow(conv) == 0) stop("No taxa models converged.")
#
#   ## False Discovery Rate adjustment
#   conv %>% dplyr::mutate(FDR_Pval = stats::p.adjust(P_val, method = "BH") %>% round(4))
# }

## Add in summary measures of diversity measures?
## It could just be the "obj"
## Function for summarizing Taxa for all Taxa or nonconvergent taxa from NB_summary
N.con <- function(x, obj, ...){

  obj <- rlang::enquo(obj)

  if(missing(...)){ ## Don't always need ... (the covariates)
    x %>%
      taxa_summarize(!!obj) %>%
      return()
  }
  else{
    x %>% dplyr::select(.data$Taxa,!!obj, cov_str(...)) %>%
      dplyr::group_by_if(is.factor) %>%  ## Groups by covariates when given
      taxa_summarize(!!obj) %>%
      return()
  }
}

## Function for making formula from "..." in functions
formula_fun <- function(...){
 covs <-  rlang::quos(...) %>%           ## Making "..." a quosure
    rlang::splice() %>%            ## Splicing quosure
    unlist %>%            ## Unlisting
    stringr::str_flatten(collapse = " + ") %>%            ## Combining elements of list with "+"
    stringr::str_replace_all(pattern = "~", replacement = "") ## Removing "~" (why is it there?)

    paste("cts ~ ",covs)  ## pasting with "cts ~ " for final formula
}

## Function for Taxa Summaries
taxa_summarize <- function(x,obj){
  obj <- rlang::enquo(obj)
  x %>%
    dplyr::summarise(n = dplyr::n(), Percent_0 = round(sum(!!obj == 0, na.rm=TRUE)/length(!!obj)*100,2),
                     Mean = mean(!!obj, na.rm=TRUE), SD = stats::sd(!!obj, na.rm=TRUE),
                     Median = stats::median(!!obj, na.rm=TRUE), IQR = stats::IQR(!!obj, na.rm=TRUE),
                     Percentile_5th  = round(stats::quantile(!!obj, probs = 0.05, na.rm=TRUE),2),
                     Percentile_10th = round(stats::quantile(!!obj, probs = 0.10, na.rm=TRUE),2),
                     Percentile_25th = round(stats::quantile(!!obj, probs = 0.25, na.rm=TRUE),2),
                     Percentile_75th = round(stats::quantile(!!obj, probs = 0.75, na.rm=TRUE),2),
                     Percentile_90th = round(stats::quantile(!!obj, probs = 0.90, na.rm=TRUE),2),
                     Percentile_95th = round(stats::quantile(!!obj, probs = 0.95, na.rm=TRUE),2))
}

## Adds main effect estimates to their intercept before exponentiating
coef_est <- function(m, ..., RE = NULL){

  if(is.null(RE)){
    covs <- data.frame(Coef = names(m$coefficients),
                       Est = m$coefficients,
                       Cov_Type = m$coefficients %>%
                         names %>%
                         cov_type(...))
  } else{
    covs <- data.frame(Coef = names(stats::coef(summary(m))[,1] ),
                       Est = m@beta,
                       Cov_Type = stats::coef(summary(m))[,1] %>%
                         names %>%
                         cov_type(...))
  }

    ## Need to move through each interaction in case there
    ## is factor with 3 or more levels

    for(i in seq(1,nrow(covs))){
      ## detecting interactions
      if( stringr::str_detect(covs$Cov_Type[i], ".int") ){
        ## interaction pieces
        int. <- covs$Coef[i] %>% as.character %>% stringr::str_split(":") %>% unlist

        ## Exact tring matches of the main effects
        ## (there will only be 2 for any interaction piece)
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
ci_est <- function(m, ..., RE = NULL){

  # Updating estimates in case there are interactions
  covs <- coef_est(m = m, ..., RE = RE)

  CI <- matrix(rep(0, nrow(covs)*2), nrow = nrow(covs), ncol = 2,
               dimnames = list(rownames(covs), c("2.5 %", "97.5 %")))
  # Confidence intervals
  for(r in seq(1,nrow(covs))){
    ## detecting interactions
    if( stringr::str_detect(covs$Cov_Type[r], ".int") ){
      # interaction pieces
      int. <- c(covs$Coef[r] %>% as.character %>%
                  stringr::str_split(":") %>% unlist,
                covs$Coef[r] %>% as.character)

      # covariance matrix
      s <- summary(m)$cov.scaled
      # pullin out rows and cols with interaction and main effect covariances
      cl <- colnames(s) %in% int. ; ro <- rownames(s) %in% int. ; s <- s[ro,cl]

      # summing vars and covars
      v <- 0
      for(i in seq(1,nrow(s))){
        for(j in seq(1,ncol(s))){

          if(j >= i){
            p <- ifelse(i == j, s[i,j], 2*s[i,j])

            v <- v + p
          }
        }
      }

      # Wald interval for interaction
      ci <- covs$Est[r] + c(-1,1)*1.96*sqrt(v)

    } else{

      # Coef summary table
      s <- summary(m)$coef

      # Wald interval for main effects
      ci <- s[r,1] + c(-1,1)*1.96*s[r,2]
    }
    CI[r,] <- ci
  }

  ## Exponentiate and round
  CI %<>% exp %>% round(4)

  ## Make neat for the table
  paste("(", CI[,1], ", ", CI[,2], ")", sep = "")

}

## Breaks whatever '...' is down to a character string of each covariate
## ie. 'Group, Sex*Age' -> c("Group","Sex","Age")
cov_str <- function(...){

  suppressWarnings(
    Cov <- rlang::quos(...) %>%           ## Making "..." a quosure
      rlang::splice() %>%                   ## Splicing quosure
      unlist %>%                   ## Unlisting
      stringr::str_split(pattern = "[+]", simplify = T) %>%
      stringr::str_split(pattern = "[*]", simplify = T) %>%
      stringr::str_replace_all(pattern = "~", replacement = "") %>%
      stringr::str_trim(side = "both") ## Making strings nice to read
  )
  Cov <- Cov[stringr::str_length(Cov) > 0] ## removing empty strings
  Cov
}

## Helper function in cov_type for interaction covs
cov_int <- function(x,c){
  ifelse(
    stringr::str_split(x, pattern = ":") %>%
      lapply(function(x,c) sum(x %in% c), c = c) == 1,"q*c.int",

    ifelse(
      stringr::str_split(x, pattern = ":") %>%
        lapply(function(x,c) sum(x %in% c), c = c) > 1,"q*q.int", "c*c.int"))
}

## Function to identify type of covariate (categ, quant, c*c.int, q*c.int, q*q.int)
cov_type <- function(coefs, ...){

    Cov <- cov_str(...)
    ## Labeling coef types
    ifelse(coefs %in% c("(Intercept)", "(Intercept):1"), "Intercept", ## Intercept
           ifelse(coefs == "(Intercept):2", "rho",
                  ifelse(coefs %in% Cov, "quant", ## Exact matches will be quants
                         ifelse(stringr::str_detect(coefs, ":"), ## Interactions
                                cov_int(coefs,Cov) %>% suppressWarnings,
                                "categ"))))
}

## Calculates TypeII_SS for each variable for covergent summary
Anova_SS <- function(cov_sum, nb, SS_type, test){
  A <- car::Anova(nb, type = SS_type, test.statistic = test)

  P <- sapply(rownames(A),
              grep, cov_sum$Coef %>% as.character)

  TT <- NULL

  if(is.matrix(P)){
    TT[P] <-  A[,3]
  } else{
    for(i in seq(1,length(P))){
      TT[ P[[i]] ] <- A[i,3]
    }
  }

  cov_sum$LRT <- TT

  cov_sum
}


## A version of cov_type that will pick out RE's that are also main effects
## Don't need this right now

# cov_type <- function(coefs, ..., RE){
#
#   ## Identifies random effects separately from other covariates
#   if(!is.null(RE)){
#     RE %<>% stringr::str_replace_all("\\s+", "") ## removes any white space from a string
#     RE %<>% substr(4L, stringr::str_length(RE) - 1) ## extracts the RE term / cov
#
#     Cov <- cov_str(...)
#
#     ## Detects if a different covariates = RE AND start with RE characters
#     ## This will mess up the labeling below
#     if(sum( c(RE == Cov,  any(startsWith(Cov, RE) & !endsWith(Cov, RE))) ) > 1 ){
#       stop("Random effect is not distinguishable from other covariates.\n Please keep covariate names as distinct as possible.")
#     }
#
#     ## Labeling coef types
#     ifelse(coefs == "(Intercept)", "Intercept", ## Intercept
#            ## Exact matches that aren't RE will be quants
#            ifelse((coefs %in% Cov) & (coefs != RE),"quant",
#                   ## Should only match quant REs exactly
#                   ifelse(coefs == RE, "RE_quant",
#                          ifelse(startsWith(coefs, RE), "RE_categ",
#                                 ## Interactions
#                                 ifelse(stringr::str_detect(coefs, ":"),
#                                        cov_int(coefs,Cov) %>% suppressWarnings,
#                                        ## Only categoricals are left
#                                        "categ")))))
#   } else{
#     Cov <- cov_str(...)
#     ## Labeling coef types
#     ifelse(coefs == "(Intercept)", "Intercept", ## Intercept
#            ifelse(coefs %in% Cov, "quant", ## Exact matches will be quants
#                   ifelse(stringr::str_detect(coefs, ":"), ## Interactions
#                          cov_int(coefs,Cov) %>% suppressWarnings,
#                          "categ")))
#   }
#
# }

