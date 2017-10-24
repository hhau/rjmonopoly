# utilities

#' Return all samples in list, of length N
#'
#' @param sample_list list containing vectors of various length
#' @param N Length of vectors you wish to extract
#'
#' @return Matrix with N columns and an unknown number of rows
returnSamplesOfLengthN <- function(sample_list, N) {

  N_samples <- 0
  N_possible <- length(sample_list)
  samples_of_length_N <- matrix(NA, nrow = N_possible, ncol = N)

  for (ii in 1:N_possible) {
    if (length(sample_list[[ii]]) == N) {
      N_samples <- N_samples + 1
      samples_of_length_N[N_samples,] <- sample_list[[ii]]
    }
  }

  res <- samples_of_length_N[1:N_samples,]
  return(res)
}

#' Generate full set of Initial values
#'
#' Generates the full set of initial values for Gamma, so as to use them as
#' means for our independent proposals for the extra dimension when we jump.
#' This is done to avoid breaking the reversibility of the markov chain.
#'
#' @param y vector of (scaled) y values
#' @param x vector of (scaled) x values
#' @inheritParams calcPropRatio
#' @param a optional lower bound for monotonicity check, but given data should
#' scaled, shouldn't need to be adjusted
#' @param b optional upper bound for monotonicity check, but given data should
#' scaled, shouldn't need to be adjusted
#'
#' @return vector of initial gamma values, using the ideal fit from MonoPoly
.getFullInitial  <- function(y, x, d_max, a = 0, b = 1) {
  full_initial <- coef(MonoPoly::monpol(y ~ x, degree = d_max, a = a, b = b))
  return(full_initial)
}

# Given a series of lists as arguments, merge the lists from left to right (so
# that the right-most values override the left-most)
# thanks github/mbertolacci
.extend_list <- function(...) {
  lists <- list(...)
  output <- lists[[1]]
  for (value in lists[2 : length(lists)]) {
    for (name in names(value)) {
      output[[name]] <- value[[name]]
    }
  }
  return(output)
}

# unscale a scaled vector, with the list from the output object
.unscale <- function(vec_scl, scl_list, ...) {
  res <- (vec_scl * scl_list$scl_range) + scl_list$min
  res
}