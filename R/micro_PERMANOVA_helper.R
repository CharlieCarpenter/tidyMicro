###############################
##
## Project: tidyMicro
##
## Purpose: micro_PERMANOVA helper functions
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

## Function for making adonis formula from "..." in functions
adonis_formula <- function(...){
  rlang::quos(...) %>%           ## Making "..." a quosure
    rlang::splice() %>%            ## Splicing quosure
    unlist %>%            ## Unlisting
    stringr::str_flatten(collapse = " + ") %>%            ## Combining elements of list with "+"
    stringr::str_replace_all(pattern = "~", replacement = "") %>%     ## Removing "~" (why is it there?)
    #  paste("beta_div", ., sep = "~") %>%
    # as.formula %>%
    return()
}
