#' Soybean Pan-African Trials data set
#'
#' This data set belongs to the Soybean Pan-African Trials (PAT). This subset has
#' the best linear unbiased estimates of grain yield (GY), plant height (PH) and
#' number of days to maturity (NDM) of 65 soybean genotypes evaluated over 19 environments.
#' The complete data set is available at Araújo et al. (2025) (check references).
#' It contains the empirical best linear unbiased estimates of
#' grain yield (GY), plant height (PH) and number of days to maturity (NDM)
#' from 65 experimental genotypes evaluated in 19 locations.
#'
#' @format ## `soy_pat`
#'  A data frame with 540 rows and 5 columns:
#' \describe{
#'   \item{Env}{19 environments}
#'   \item{Gen}{65 experimental genotypes}
#'   \item{Plant_Height}{395 BLUEs - Plant height measurements}
#'   \item{Grain_Yield}{525 BLUEs - Grain yield measurements}
#'   \item{Days_to_Maturity}{312 BLUEs - Number of days to maturity}
#' }
#'
#' @source
#' \describe{
#' Araújo, M. S., Chaves, S., Ferreira, G. N. C., Chigeza, G.,
#' Leles, E. P., Santos, M. F. S., Diers, B. W., Goldsmith, P.,
#' and Pinheiro, J. B. (2025). High-resolution soybean trial data supporting the
#' expansion of agriculture in Africa. \emph{Scientific Data}, 12:1908.
#' \doi{10.1038/s41597-025-06190-3}
#' }
#'
"soy_pat"
