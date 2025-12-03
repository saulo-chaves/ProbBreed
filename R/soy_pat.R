#' Soybean PAT real data set
#'
#' This dataset belongs to the Soybean Pan-African Trials (PAT) which
#' evaluate 65 soybean genotypes across 19 environments (Araújo et al. 2025).
#' It contains the empirical best linear unbiased estimates of genotypic means of
#' grain yield (GY), plant height (PH) and number of days to maturity (NDM)
#' from 65 experimental genotypes evaluated in 19 locations.
#'
#' @format ## `soy_pat`
#'  A data frame with 540 rows and 5 columns:
#' \describe{
#'   \item{Env}{19 environments}
#'   \item{Gen}{65 experimental genotypes}
#'   \item{Plant_Height}{395 EBLUEs (phenotypes) - Plant height measurements}
#'   \item{Grain_Yield}{525 EBLUEs (phenotypes) - Grain yield measurements}
#'   \item{Days_to_Maturity}{312 EBLUEs (phenotypes) - Number of days to maturity}
#' }
#'
#' @source
#' \describe{
#' Araújo, Mauricio S., Saulo F. Chaves, Gérson N. C. Ferreira, Godfree Chigeza,
#' Erica P. Leles, Michelle F. Santos, Brian W. Diers, Peter Goldsmith, and José B. Pinheiro.
#' 2025. "High-Resolution Soybean Trial Data Supporting the Expansion of
#' Agriculture in Africa." Scientific Data.
#' doi: https://doi.org/10.1038/s41597-025-06190-3.}
#'
"soy_pat"
