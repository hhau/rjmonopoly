#mh helpers

#' Generate a proposal for our coefficients in the orthonormal space.
#'
#' @param gamma_old Vector containing the current value for gamma in the Markov
#' chain
#' @param innov_sigma Innovation variance associated with the proposal
#' distribution for sigma
#'
#' @return A proposed realisation for gamma
#'
genSingleGammaProp <- function(gamma_old, innov_sigma) {
  gamma_proposal <- rnorm(length(gamma_old), mean = gamma_old, sd = innov_sigma)
  return(gamma_proposal)
}

#' Generates a proposal for our variance parameter
#'
#' @param var_prev Current value of variance
#' @param innov_sd_var Variance associated with the innovation distribution for
#' our variance parameter
#'
#' @return A proposed realisation for the variance
#'
genVarProp <- function(var_prev, innov_sd_var) {
  var_prop <- rlnorm(1, meanlog = log(var_prev), sdlog = innov_sd_var)
  return(var_prop)
}



#' Calculate the section of the acceptance probabilty pertaining to the
#' posterior distribution of our parameters
#'
#' @param gamma_prop Vector containing proposal for Gamma
#' @param gamma_curr Vector containing the current value for Gamma
#' @param var_prop Proposed value for the variance
#' @param var_curr Current value for the variance
#' @param Q_full Full size (n_obs * (d_max + 1)) matrix for Q
#' @param d_prop Proposed polynomial degree
#' @param d_curr Current polynomial degree
#' @param y_vec Vector of y observations
#' @param dim_prior_vec Vector containing prior values for our dimension
#' parameter
#' @param R_inv_full Full size (d_max + 1) * (d_max + 1) basis inverse matrix
#'
#' @return A value for the posterior acceptance ratio
calcPostRatio <- function(gamma_prop, gamma_curr, var_prop, var_curr, Q_full,
                          d_prop, d_curr, y_vec, dim_prior_vec, R_inv_full) {

  mu_prop <- Q_full[, 1:(d_prop + 1)] %*% gamma_prop
  mu_curr <- Q_full[, 1:(d_curr + 1)] %*% gamma_curr

  beta_prop <- gammaToBeta(gamma_prop,
                           R_inv = R_inv_full[1:(length(gamma_prop)), 1:(length(gamma_prop))])

  mono_modifier <- 1
  if (!MonoPoly::ismonotone(beta_prop, a = 0, b = 1)) {
    mono_modifier <- 0

  }
  # flat prior on the variance, not always a good idea, here is where one
  # would change the prior on the variance.
  var_prior_ratio <-  (1 / var_prop) / (1 / var_curr)

  # our own prior on the dimension
  dim_prior_ratio <- dim_prior_vec[d_prop] / dim_prior_vec[d_curr]

  # likelihood term, evaled on the logscale and co-ordinatewise, for speed and
  # precision
  ratio <- sum(dnorm(x = (y_vec), mean = (mu_prop), sd = sqrt(var_prop), log = TRUE)) -
    sum(dnorm(x = (y_vec), mean = (mu_curr), sd = sqrt(var_curr), log = TRUE))

  ratio <- exp(ratio)
  return(ratio * var_prior_ratio * dim_prior_ratio * mono_modifier)
}

#' Calculate the part of the acceptance probability pertaining to the proposal
#' distributions. This function will make absolutely no sense unless looked at
#' after one has consulted the maths.
#'
#' @param gamma_k_prime_prop proposed value of gamma from the additional
#' dimension
#' @param k_prime_init mean value for the independent proposal we are using
#' for the additional regression parameter
#' @param innov_sigma innovation variance that pertains to the gamma proposal
#' @param d_prop proposed polynomial degree
#' @param d_curr current polynomial degree
#' @param d_min minimum allowable polynomial degree
#' @param d_max maximum allowable polynomial degree
#' @param var_prop proposed variance
#' @param var_curr current variance
#' @param innov_sd_var innovation variance pertaining to the variance proposal
#'
#' @return A value for the proposal acceptance ratio
calcPropRatio <- function(gamma_k_prime_prop, k_prime_init, innov_sigma,
                          d_prop, d_curr, d_min, d_max, var_prop, var_curr,
                          innov_sd_var) {

  #27/7/17 edit
  # Variance proposal is no longer symmetric random walk, so need that term.
  # I'M INCLUDING THIS TERM IN HERE, WHEN IT TECHINCALLY SHOULD BE IN THE
  # PropRatio FUNCTION JUST FOR CONVIENCED, IT ALSO MEANS ARE PRIORS ARE NOW
  # ON THE LOG SCALE FOR THE VARIANCE PARAMETER.

  var_proposal_ratio <- dlnorm(x = var_curr, meanlog = log(var_prop), sdlog = innov_sd_var) /
    dlnorm(x = var_prop, meanlog = log(var_curr), sdlog = innov_sd_var)

  if (d_prop == d_curr) {
    # if we don't propose a dimension change, then we don't have the dimension
    # jumping factor, nor the generation associated with the extra coefficient
    # so the proposal ratio is just 1 (everything else is random walk)

    return(1 * var_proposal_ratio)
  }

  q_gamma_k_prime <- dnorm(x = gamma_k_prime_prop, mean = k_prime_init,
                           sd = innov_sigma)

  dim_prop_ratio <- 1

  if ((d_prop == d_max) & (d_curr == (d_max - 1))) {
    there <- 1/3
    back <- 1/2
    dim_prop_ratio <- back / there

  } else if ((d_prop == (d_max - 1)) & (d_curr == d_max )) {
    there <- 1/2
    back <- 1/3
    dim_prop_ratio <- back / there

  } else if ((d_prop == d_min) & (d_curr == (d_min  + 1))) {
    there <- 1/3
    back <- 1/2
    dim_prop_ratio <- back / there

  } else if ((d_prop == (d_min + 1)) & (d_curr == d_min)) {
    there <- 1/2
    back <- 1/3
    dim_prop_ratio <- back / there

  }

  # I have the math to prove this needs to be in here twice
  res <- (1 / (q_gamma_k_prime ^ 2)) * dim_prop_ratio  * var_proposal_ratio

  return(res)

}


#' generates a discrete distribution over the allowable degree space, uses
#' the binomial distribution to adjust the distribution
#'
#' @param p  Probability passed to \link{dbinom}
#' @inheritParams calcPropRatio
#'
#' @return vector containing prior values
#'
genDimPrior <- function(d_min, d_max, p) {
  temp <- dbinom(x = 1:d_max, size = d_max, prob = p)
  temp[1:d_min] <- 0
  res <- temp / sum(temp)

  return(res)

}
