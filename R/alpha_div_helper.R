###############################
##
## Project: tidyMicro
##
## Purpose: Helper functions for alpha_div
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


## Sobs ----
sobs <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  ## Takes one bootstrap sample and adds up unique taxa
  sobs_fun <- function(data_set, min_lib.){
    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    sample(Seq, size = min_lib., replace = TRUE) %>%
      unique %>%
      length
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply(~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) sobs_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>% dplyr::rename(Sobs = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}

## Bias corrected Chao1 ----
chao1 <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  ## Takes one bootstrap sample and adds up unique taxa
  chao1_fun <- function(data_set, min_lib.){
    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    ## rarefied bootstrap sample
    boot.seq <- sample(Seq, size = min_lib., replace = TRUE)

    ## Number of species appearing exactly once
    sing <- sum(table(boot.seq) == 1)

    ## Number of species appearing exactly twice
    doub <- sum(table(boot.seq) == 2)

    ## Bias corrected Chao1
    length(unique(boot.seq)) +
      ((unique(data_set$Total) - 1)/unique(data_set$Total))*(sing*(sing - 1)/(2*(doub + 1)))
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply(~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) chao1_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>%
    dplyr::rename(Chao1 = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}

## Goods ----
goods <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  goods_fun <- function(data_set, min_lib.){
    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    ## rarefied bootstrap sample
    boot.seq <- sample(Seq, size = min_lib., replace = TRUE)

    ## Number of species appearing exactly once
    sing <- sum(table(boot.seq) == 1)

    ## Coverage
    (1 - (sing / length(Seq)))*100
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply( ~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) goods_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>%
    dplyr::rename(Goods = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}

## ShannonE ----

# Uses log2
shannonE <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  shannonE_fun <- function(data_set, min_lib.){

    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    ## rarefied bootstrap sample
    boot.seq <- sample(Seq, size = min_lib., replace = TRUE)

    p_i <- table(boot.seq) / min_lib.
    p_i <- p_i[p_i > 0]

    -sum(p_i * log2(p_i)) / log2(length(unique(boot.seq)))
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply(~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) shannonE_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>%
    dplyr::rename(ShannonE = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}

## ShannonH ----
## Uses log2
shannonH <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  shannonH_fun <- function(data_set, min_lib.){

    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    ## rarefied bootstrap sample
    boot.seq <- sample(Seq, size = min_lib., replace = TRUE)

    p_i <- table(boot.seq) / min_lib.
    p_i <- p_i[p_i > 0]

    -sum(p_i * log2(p_i))
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply(~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) shannonH_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>%
    dplyr::rename(ShannonH = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}

## SimpsonD ----
simpsonD <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  simpsonD_fun <- function(data_set, min_lib.){

    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    ## rarefied bootstrap sample
    boot.seq <- sample(Seq, size = min_lib., replace = TRUE)

    p_i <- table(boot.seq) / min_lib.

    1 / sum(p_i^2)
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply(~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) simpsonD_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>%
    dplyr::rename(SimpsonD = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}

## SimpsonE ----
simpsonE <- function(data_set, iter){

  min_lib <- min(data_set$Total) ## Rarefaction size of bootstrap samples

  simpsonE_fun <- function(data_set, min_lib.){

    ## The "raw" sequence
    Seq <- rep(data_set$Taxa, data_set$cts)

    ## rarefied bootstrap sample
    boot.seq <- sample(Seq, size = min_lib., replace = TRUE)

    p_i <- table(boot.seq) / min_lib.

    1 / (sum(p_i^2) * length(unique(boot.seq)))
  }

  ## Takes 'iter' bootstrap samples w/in Lib and takes their mean
  data_set %>%
    plyr::ddply(~ .data$Table + .data$Lib,
                .fun = function(data_set, min_lib. = min_lib, iter. = iter){
                  sapply(seq(1,iter.), function(x) simpsonE_fun(data_set, min_lib)) %>%
                    mean
                }
    ) %>%
    dplyr::rename(SimpsonE = .data$V1, Table = .data$`.data$Table`, Lib = .data$`.data$Lib`)
}
