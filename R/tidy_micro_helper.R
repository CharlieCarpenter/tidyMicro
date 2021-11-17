###############################
##
## Project: tidyMicro
##
## Purpose: tidy_micro helper functions
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

# Helpful functions
`%nin%` <- Negate(`%in%`)

## Naming global variables
  ## These are all taxa information pieces that we enquote (e.g. rlang::enquo(ra))
utils::globalVariables(c("cts", "bin", "clr", "ra"))

my_clr <- function(otu,total){

  x <- sweep(otu,1,1/total, "+") ## Adds 1/total to everything to remove 0s

  ## Log(x) and the Geometric mean
  l_x <- log(x)
  gm <- apply(l_x,1, function(x) exp(mean(x)) )

  clr <- sweep(l_x,1,gm) %>%
    as.data.frame

  clr
}

my_ra <- function(otu,total){
  ## Dividing everything by the first col (the root)
  ra <- sweep(otu, 1, total, "/") %>%
    ## Multiplying by 100 to get percentages
    sweep(1,100,"*") %>%
    as.data.frame

  ra
}

my_bin <- function(otu,total){
  bin <- ifelse(otu > 0,1,0) %>% as.data.frame

  bin
}

## Function for shortening Taxa names. Removes "Bacteria/"
## Keeps abbreviated phylum and Rank name.
## i.e. Class level OTUs will have names like Firm:Clostrdia
# taxa_names <- function(string){
#   sub("Bacteria/", "", string, ignore.case = FALSE, fixed = TRUE) %>% ## Removes Bacteria
#     sapply( function(x) strsplit(x, "/")) %>% ## Splits
#     lapply( function(x){
#       if(length(x) == 2L){
#         stringr::str_sub(x[1L],1L,5L) %>%
#           paste(x[2],sep = ":")
#       }
#        else {
#          if(length(x) > 2L){ ## This will only select Ranks lower than Class (longer names)
#           stringr::str_sub(x[1L], 1L, 5L) %>% ## Abreviates phylum
#              paste(stringr::str_sub(x[2L], 1L, 5L), sep = ":") %>% ## Abreviates Class
#              paste(x[length(x)],sep = ":")} ## Pastes with Rank level OTU name
#       else x ## Just prints phylum if Rank is phylum
#     }}) %>%
#     unlist %>%
#     as.character()
# }

mul_otu_long <- function(in_OTU, .meta){
  if(in_OTU %>% colnames %>% anyDuplicated){
    stop("OTU table contains repeated suject/library names.")
  }

  Lib <- colnames(in_OTU)[-1] ## Library Names
  if(!any(Lib %in% .meta$Lib)){
    stop("No library names match the column names of provided OTU table")
  }

  otu <- t(in_OTU)[-1,] %>% ## removing OTU Names
    apply(2, function(x) as.numeric(x)) %>%  ## Making cols numeric
    as.data.frame

  names(otu) <- in_OTU[,1] ## Replacing OTU Names
  total <- rowSums(otu) ## the total counts / seq depth

  ## Calculating transformations and reattaching on Lib
  ra <- my_ra(otu, total) %>% dplyr::mutate(Lib = Lib)
  clr <- my_clr(otu, total) %>% dplyr::mutate(Lib = Lib)
  bin <- my_bin(otu, total) %>% dplyr::mutate(Lib = Lib)

  ## Reattaching on Lib
  otu %<>% dplyr::mutate(Lib = Lib)

  ## pivot long data
  m.ra <- ra %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "ra")
  m.clr <- clr %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "clr")
  m.bin <- bin %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "bin")
  m.cts <- otu %>% tidyr::pivot_longer(-Lib, names_to = "Taxa", values_to = "cts")

  tots <- rep(total,
              each = length(unique(m.ra$Taxa)))

  long_OTU <- dplyr::left_join(m.bin,m.cts,by=c("Lib","Taxa")) %>%
    dplyr::left_join(m.clr,by=c("Lib","Taxa")) %>%
    dplyr::left_join(m.ra,by=c("Lib","Taxa")) %>%
    dplyr::mutate(Total= tots) %>% ## Creating total for each taxa
    dplyr::select(Lib, Taxa, dplyr::everything())

  suppressWarnings(
    long_OTU %>%
      dplyr::select(Lib, Taxa, dplyr::everything()) %>%
      dplyr::full_join(.meta, by="Lib")
  ) %>%
    return()
}
