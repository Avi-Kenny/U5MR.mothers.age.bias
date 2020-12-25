#' Convert dates (year+month) to CMC
#' @param year Calendar year
#' @param month Calendar month
#' @return Date in CMC format

dates_to_cmc <- function(year, month) {
  12*(year-1900)+month
}



#' Convert CMC to dates (year+month)
#' @param cmc CMC-formatted date
#' @return Date in CMC format
#' @return A list containing: \cr
#'     * `year`: calendar year  \cr
#'     * `month`: calendar month

cmc_to_dates <- function(cmc) {
  month <- ifelse(mod(cmc, 12)!=0, mod(cmc, 12), 12)
  return(list(
    year = 1900 + (cmc-month)/12,
    month = month
  ))
}

