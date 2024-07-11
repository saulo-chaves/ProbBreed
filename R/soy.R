##' Soybean real data set
##'
##' This dataset belongs to the USDA Northern Region Uniform Soybean Tests,
##' and it is a subset of the data used by Krause et al. (2023). It contains the
##' empirical best linear unbiased estimates of genotypic means of the seed yield
##' from 39 experimental genotypes evaluated in 14 locations. The original data, available at the package
##' `SoyURT`, has 4,257 experimental genotypes evaluated at 63 locations and
##' 31 years resulting in 591 location-year combinations (environments) with
##' 39,006 yield values.
##'
##' @format ## `soy`
##'  A data frame with 823 rows and 3 columns:
##'  \describe{
##'    \item{Loc}{14 locations}
##'    \item{Gen}{39 experimental genotypes}
##'    \item{Y}{435 EBLUEs (phenotypes)}
##'  }
##'
##' @source
##'  \describe{
##'  Krause MD, Dias KOG, Singh AK, Beavis WD. 2023. Using soybean
##'  historical field trial data to study genotype by environment
##'  variation and identify mega-environments with the integration
##'  of genetic and non-genetic factors. bioRxiv : the preprint server for biology.
##'  doi: https://doi.org/10.1101/2022.04.11.487885
##'  }
##'
"soy"
