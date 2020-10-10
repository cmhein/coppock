#' Coppock curve
#' @description
#' A momentum indicator, originally a "Very Long Term (VLT) Momentum Index",
#' defined as the 10-month weighted moving average (WMA) of the sum of a series'
#' 11-month and 14-month rates of change (ROC).
#' @details
#' Using the default parameters, the Coppock curve is defined as the
#' \code{10-WMA of (ROC[11] + ROC[14])}.
#' @param x A vector
#' @param n1 Rate of change period #1
#' @param n2 Rate of change period #2
#' @param n3 Weighted moving average period
#' @return A vector
#' @export
coppock_curve <- function(x, n1 = 14L, n2 = 11L, n3 = 10L) {
  x1 <- dplyr::lag(x, n1)
  x2 <- dplyr::lag(x, n2)
  r1 <- (x - x1) / x1
  r2 <- (x - x2) / x2
  TTR::WMA(r1 + r2, n3)
}

#' Solve Coppock curve for each \code{x}
#' @description
#' Solve for each value of \code{x} that gives Coppock curve value \code{y}.
#' @details
#' The calculation of the Coppock curve is such that changing one value of
#' \code{x} would affect the Coppock curve at other points. Therefore, the
#' solution for each \code{x} is independent and is found by holding constant
#' all other values of \code{x}.
#' @inheritParams coppock_curve
#' @param y Coppock curve value
#' @return A vector
#' @export
coppock_curve_solve_each <- function(x, n1 = 14L, n2 = 11L, n3 = 10L, y = 0) {
  # This is the first portion of the Coppock curve calculation, lacking only
  # the WMA component
  x1 <- dplyr::lag(x, n1)
  x2 <- dplyr::lag(x, n2)
  r1 <- (x - x1) / x1
  r2 <- (x - x2) / x2

  # The approach here is to compute the WMA using a period one less than the
  # desired. We'll then find the x-value that, when added to the WMA calcuation
  # gives us the y-value we seek. Note that in order to line-up the vector
  # elements for computation below, we shift the result of this calculation to
  # the right one element.
  n4 <- n3 - 1L
  wma <- TTR::WMA(r1 + r2, n4)
  wma <- c(NA, wma[-length(wma)])

  # Convert from weighted moving average to weighted moving sum by multiplying
  # by the divisor that was used to compute the average. Essentially, we're
  # backing-up a step so that we can extend the WMA period by one and add to
  # it the datapoint for which we're trying to solve.
  d <- n4 * (n4 + 1L) / 2L # sum of 1:n4
  wms <- wma * d

  # Solve for the x-value that would give the y-value provided by argument. This
  # is the Coppock curve equation rearranged (set equal to x). The wms parameter
  # is the weighted moving sum of the elements (sans the x for which we're
  # solving) used to compute the weighted moving avarege. We add the WMA period
  # (n3) to the WMA divisor (d) because the divisor was computed above using one
  # less than the desired WMA period.
  (((y * (d + n3) - wms) * x1 * x2 / n3) + (2 * x1 * x2)) / (x1 + x2)
}

#' Solve Coppock curve for next \code{x}
#' @description
#' Solve for next value of \code{x} that would give Coppock curve value
#' \code{y}.
#' @inheritParams coppock_curve_solve_each
#' @return A numeric
#' @export
coppock_curve_solve_next <- function(x, n1 = 14L, n2 = 11L, n3 = 10L, y = 0) {
  # Subset x so that we only compute lags and returns on data we need
  x_len <- length(x)
  x <- x[(x_len - max(n1, n2)):x_len]

  # Append a zero to the x vector as a simple means of solving for the next
  # x-value that would give the y-value provided by argument. Note that the zero
  # is arbitrary and does not affect the calculation.
  v <- coppock_curve_solve(c(x, 0), n1, n2, n3, y)
  v[length(v)]
}
