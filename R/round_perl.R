#' Round a number using Perl-style rounding
#'
#' This function rounds a numeric value to the nearest integer using
#' Perl-style rounding (i.e., adding 0.5 and taking the floor).
#'
#' @param number A numeric value to be rounded.
#' @return An integer rounded from `number`.
#' @examples
#' round_perl(2.3)
#' round_perl(2.7)
#' @export
round_perl <- function(number) {
  floor(number + 0.5)
}