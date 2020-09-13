###############################
##
## Project: tidyMicro
##
## Purpose: filter_funs helper function
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


## Filter functions to be applied within tidy_mibi
mul_filter <- function(otu, prev_cutoff, ra_cutoff, total, ex, rr){

  if(prev_cutoff > 100 | prev_cutoff < 0) stop("Prevalence cutoff must be between 0 and 100")
  if(ra_cutoff > 100 | ra_cutoff < 0) stop("Relative abundance cutoff must be between 0 and 100")
  if(!is.null(ex) & !is.character(ex)) stop("'exclue_taxa' should be a character or left as NULL")

  ## Prints out curent filtering Table if any filtering is done
  if( any( c(!is.null(ex), prev_cutoff > 0, ra_cutoff > 0) ) ){
    message("Filter for ", rr, " counts\n\n")
  }
  ## Different steps
  if(!is.null(ex)) otu %<>% ex_filter(ex)
  if(prev_cutoff > 0) otu %<>% prev_filter(prev_cutoff)
  if(ra_cutoff > 0) otu %<>% ra_filter(ra_cutoff, total)
  otu
}

prev_filter <- function(otu, prev_cutoff){
  starting_cts <- sum(otu)
  prevs <- apply(otu, 2, function(x) 100*sum(x > 0)/length(x))
  ## calculating prevalence and converting to a percent

  message("Prevalence cutoff: ", prev_cutoff, "% (i.e., at least ", prev_cutoff,
      "% of libaries must be represented to keep OTU)\n")
  message("Found ", ncol(otu), " OTUs.\n")

  if(sum(prevs < prev_cutoff) == ncol(otu)) stop("All OTUs would be excluded with this cutoff. Halting operation")
  if(sum(prevs < prev_cutoff) > 0) {
    include <- otu[,prevs >= prev_cutoff]
    exclude <- otu[,prevs < prev_cutoff]

    if(is.data.frame(exclude) | is.matrix(exclude)) other <- rowSums(exclude)
    else if(is.numeric(exclude)) other <- exclude

    if ("Other" %in% names(include)) {
      include[,"Other"] <- include[,"Other"] + other
      message("Found 'Other' category in input data.\n")
    } else {
      include$Other <- other
      message("Created new 'Other' Taxa in OTU table for ", ncol(otu)+1, " counts.\n")
    }

    if("Other" %in% names(exclude)){
      message("Collapsed ",sum(prevs < prev_cutoff)-1," OTUs into 'Other' in OTU table.\n")
    }

    if("Other" %nin% names(exclude)){
      message("Collapsed ",sum(prevs < prev_cutoff)," OTUs into 'Other' in OTU table.\n")
    }

    message("Converted ",sum(exclude)," counts to 'Other' in otu category.\n")
    message("Remaining OTUs: ",ncol(include)," (Including 'Other').\n\n")

    ending_cts <- sum(include)
    if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")

    return(include)
  }
  else {
    message("No OTUs with prevalence < ", prev_cutoff, "\n\n")
    return(otu)
  }
}

ra_filter <- function(otu, ra_cutoff, total) {
  starting_cts <- sum(otu) ## Starting counts

  message("Relative abundance cutoff: ",ra_cutoff,
      "% (i.e., at least one library must have RA > ",ra_cutoff,
      "% to keep OTU).\n")
  message("Found ", ncol(otu), " OTUs.\n")

  ra <- my_ra(otu,total) ## Relative abundance
  ras <- apply(ra, 2, max) ## Max RA
  below_cutoff <- sum(ras < ra_cutoff) #store TRUE/FALSE for whether an OTU has a max abundance < cutoff
  if(below_cutoff == ncol(ra)) stop("All OTUs would be excluded with this RA cutoff. Halting operation!")

  if(below_cutoff > 0) {
    include <- otu[,ras >= ra_cutoff]  ## OTUs above ra_cutoff
    exclude <- otu[,ras < ra_cutoff]   ## OTUs below ra_cutoff

    if(is.data.frame(exclude) | is.matrix(exclude)) other <- rowSums(exclude)
    else if(is.numeric(exclude)) other <- exclude

    if("Other" %in% names(otu)){
      message("Found 'Other' category in input data.\n")
    } else message("Created new 'Other' Taxa in OTU table for ", ncol(otu)+1, " counts.\n")

    if ("Other" %in% names(include)) {
      include[,"Other"] <- include[,"Other"] + other
    } else include$Other <- other

    ending_cts <- sum(include)
    if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")

    if("Other" %in% names(exclude)){
      message("Collapsed ",below_cutoff - 1," OTUs into 'Other' in OTU table.\n")
    }

    if("Other" %nin% names(exclude)){
      message("Collapsed ",below_cutoff," OTUs into 'Other' in OTU table.\n")
    }

    message("Converted ",sum(exclude)," counts to 'Other' otu category.\n")
    message("Remaining OTUs: ",ncol(include)," (Including 'Other').\n\n")

    return(include)
  }
  else {
    message("No OTUs with max relative abundance < ", ra_cutoff, "%.\n\n")
    return(otu)
  }
}

ex_filter <- function(otu, ex){
  starting_cts <- sum(otu)  ## Starting counts
  starting_nms <- names(otu) ## Starting names

  if(!is.character(ex)) stop("exclude_taxa must be a character string")

  if("Other" %in% ex){
    warning("Can not remove 'Other' category from OTU tables. 'Other' will be the aggregated counts of filtered OTUs")
    ex <- ex[ex != "Other"] ## Removing "Other" from ex
  }

  if(any(ex %in% colnames(otu))){

    ## Printing the unwanted taxa that were found
    unc_sum <- 0
    for(i in seq(1,length(ex))) {

      if(ex[i] %in% names(otu)){
        message("Found '", ex[i], "' category in input data.\n")

        uncounts <- otu[,ex[i]] ## Pulling out the unwanted taxa
        unc_sum <- unc_sum + uncounts ## Running counter for counts collapsed

        if("Other" %in% colnames(otu)) {  ## Combining unclassified with other if Other exists
          if(is.data.frame(uncounts) | is.matrix(uncounts)){
            otu[ ,"Other"] <- otu[ ,"Other"] + rowSums(uncounts) ## rowSums if many unclassified
          } else otu[ ,"Other"] <- otu[ ,"Other"] + uncounts ## for only a vector of unclassified
          # message("Found 'Other' category in input data.\n")
        } else {
          if(is.data.frame(uncounts) | is.matrix(uncounts)){
            otu$Other <- rowSums(uncounts) ## rowSums if many unclassified
          } else otu$Other <- uncounts  ## Making Other column if one doesn't exist
          message("Created new 'Other' category.\n")
        }
      } else{
        message("'", ex[i], "' not found in input data. ","'", ex[i], "' might have been misspelled.\n")
      }
    }

    otu <- otu[,names(otu) %nin% ex] ## Taking out ex so that it isn't there twice

    message("Found ", length(starting_nms), " OTUs.\n")
    message("Collapsed ",sum(ex %in% starting_nms)," OTUs into 'Other' in OTU table.\n")
    message("Converted ",sum(unc_sum)," counts to 'Other' otu category.\n")
    message("Remaining OTUs: ", ncol(otu)," (Including 'Other').\n\n")

  } else {
    message("No taxa found in 'exlude_taxa' found: OTU counts unchanged.\n\n")
  }

  ending_cts <- sum(otu)

  if (ending_cts != starting_cts) stop("Ending counts not equal to starting counts!\n")
  otu
}

sum_fun <- function(Exp, obj, r){

  if(!is.numeric(Exp[,obj])) stop("obj from tidy_mibi must be numeric")
  if(!is.character(r)) stop("Table must be a character")

  Exp <- Exp %>% dplyr::filter(.data$Table==r)
  n <- tapply(Exp[,obj], Exp$Taxa, stats::median) %>% names()
  m <- tapply(Exp[,obj], Exp$Taxa, stats::median) %>% as.data.frame()
  i <- tapply(Exp[,obj], Exp$Taxa, stats::IQR) %>% as.data.frame()
  Q2 <- tapply(Exp[,obj], Exp$Taxa, stats::quantile,0.025) %>% as.data.frame()
  Q9 <- tapply(Exp[,obj], Exp$Taxa, stats::quantile,0.975) %>% as.data.frame()

  if(r != "phylum") message("Note: Phylum level names are unclassified taxa.")
  data.frame(
    n = n[!is.na(m)],
    Median = m[!is.na(m)],
    IQR = i[!is.na(m)],
    Quantile_2.5 = Q2[!is.na(Q2)],
    Quantile_97.5 = Q9[!is.na(Q9)]
  ) %>% return()
}
