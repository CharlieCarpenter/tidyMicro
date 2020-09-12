###############################
##
## Project: tidyMicro
##
## Purpose: rocky_mtn helper functions
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

## Function to reorder Taxa to make sure "Other" is the last level
phyl_ord <- function(cor_set){

  if("Other" %in% unique(cor_set$phyl)){
    cor_set$phyl <- factor(cor_set$phyl,
                           levels = c(as.character(sort(unique(cor_set$phyl)[unique(cor_set$phyl) != "Other"])),
                                      "Other")
    )
  } else{
    cor_set$phyl <- factor(cor_set$phyl,
                           levels = as.character(sort(unique(cor_set$phyl)))
    )
  }

  cor_set
}


## Function to match taxa with their full phylum name for better legend
phy_fun <- function(x){
  phy <- c("Acidobacteria","Actinobacteria","Aquificae","Armatimonadetes","Bacteroidetes",
           "Bacteria", "Caldiserica", "Chlamydiae", "Chlorobi","Chloroflexi","Chrysiogenetes",
           "Cyanobacteria","Deferribacteres","Deinococcus-Thermus","Dictyoglomi",
           "Elusimicrobia","Fibrobacteres","Firmicutes","Fusobacteria","Gemmatimonadetes",
           "Lentisphaerae","Nitrospirae","Planctomycetes","Proteobacteria","Spirochaetes",
           "Synergistetes","Tenericutes","Thermodesulfobacteria","Thermomicrobia",
           "Thermotogae","Verrucomicrobia")

  phyl <- ifelse(sum(stringr::str_detect(phy, stringr::str_sub(x, 1L, 8L))) == 1,
                 phy[stringr::str_detect(phy, stringr::str_sub(x, 1L, 7L))],
                 as.character(x)
  )
  phyl
}
