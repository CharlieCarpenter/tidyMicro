% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_divs.R
\name{beta_div}
\alias{beta_div}
\title{Beta Diversity Calculations for tidy_micro}
\usage{
beta_div(micro_set, table, method = "bray")
}
\arguments{
\item{micro_set}{A tidy_micro data set}

\item{table}{Table you'd like to use when calculating alpha diversity. Your lowest level is recommended}

\item{method}{A dissimilarity method compatible with \code{\link[vegan]{vegdist}}}
}
\value{
A symmetrix distance matrix
}
\description{
Calculate beta diversities of your tidy_micro set. This function reformats the data into the original OTU table and then feeds that into the vegdist function
}
\examples{
data(bpd_phy); data(bpd_cla); data(bpd_ord); data(bpd_fam); data(bpd_clin)

otu_tabs <- list(Phylum = bpd_phy, Class = bpd_cla,
Order = bpd_ord, Family = bpd_fam)
set <- tidy_micro(otu_tabs = otu_tabs, clinical = bpd_clin) \%>\%
filter(day == 7) ## Only including the first week

## Bray-Curtis beta diversity
bray <- set \%>\% beta_div(table = "Family")

## Morisita-Horn beta diversity
horn <- set \%>\% beta_div(table = "Family", method = "horn")
}
\references{
\code{\link[vegan]{vegdist}}
}
