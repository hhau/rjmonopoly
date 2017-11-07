# rev jump helpers

#' Proposes new polynomial degree / dimension for parameter space
#'
#' @param d_current Current polynomial degree
#' @param d_min minimum allowable polynomial degree
#' @param d_max maximum allowable polynomial degree
#'
#' @return Proposed polynomial degree
dimProposer <- function(d_current, d_min, d_max) {

  if (d_current < d_min | d_current > d_max) {
    stop("current dimension is outside of allowed bounds")
  }

  if (d_current == d_min) {
    res <- sample(x = c(d_current, d_current + 1), size = 1)

  } else if (d_current == d_max) {
    res <- sample(x = c(d_current, d_current - 1), size = 1)

  } else {
    res <- sample(x = c(d_current - 1, d_current, d_current + 1), size = 1)

  }

  return(res)
}

#' Generates a proposal for Gamma, whilst considering proposed dimension
#'
#' @param current_gamma Vector containing the current value for Gamma (really
#' should rename this to gamma_curr like I have everywhere else)
#' @inheritParams calcPropRatio
#' @param R_inv_full full size (d_max + 1) * (d_max + 1) inverse matrix
#' @param full_inits full set of initial gamma values from the MonoPoly::monpol
#' fit
#'
#' @return Proposal for gamma
genGammaProp <- function(current_gamma, d_prop, d_curr,
                         innov_sigma, full_inits) {
  # make sure the full R_inv gets passed to this function, so we can pass
  # the correct dimension down

  if (d_prop == d_curr) {
    gamma_prop <- genSingleGammaProp(gamma_old = current_gamma,
                              innov_sigma = innov_sigma)

  } else if (d_prop ==  (1 + d_curr)) {
    prop_mean <- c(current_gamma, full_inits[d_prop + 1])
    gamma_prop <- genSingleGammaProp(gamma_old = prop_mean, innov_sigma = innov_sigma)

  } else if (d_prop == (d_curr - 1)) {
    prop_mean <- current_gamma[1:(d_prop + 1)]
    gamma_prop <- genSingleGammaProp(gamma_old = prop_mean, innov_sigma = innov_sigma)

  }

  # print("Double checking what gamma_proposal is:")
  # print(gamma_prop)


  return(gamma_prop)

}

#' Calculate the acceptance probability for a preposal, whilst enforcing
#' monotoncity.
#' @inheritParams calcPropRatio
#' @inheritParams calcPostRatio
#' @inheritParams genGammaProp
#'
#' @return An acceptance probability
calcAcceptProb <- function(gamma_prop, gamma_curr, d_prop, d_curr, d_min, d_max,
                           var_prop, var_curr, Q_full, full_inits, y_vec,
                           innov_sigma, dim_prior_vec, innov_sd_var = innov_sd_var,
                           R_inv_full) {

  accept_prob <- NULL

  if (d_prop == d_curr) {
    post <- calcPostRatio(gamma_prop = gamma_prop, gamma_curr = gamma_curr,
                          var_prop = var_prop, var_curr = var_curr,
                          Q_full = Q_full, d_prop = d_prop, d_curr = d_curr,
                          y_vec = y_vec, dim_prior_vec = dim_prior_vec,
                          R_inv_full = R_inv_full)
    prop <- calcPropRatio(gamma_k_prime_prop = gamma_prop[d_prop + 1],
                          k_prime_init = full_inits[d_prop + 1],
                          innov_sigma = innov_sigma, d_prop = d_prop, d_curr = d_curr,
                          d_min = d_min, d_max = d_max, var_prop = var_prop,
                          var_curr = var_curr, innov_sd_var = innov_sd_var)

    accept_prob <- min(1, post * prop)

  } else if (d_prop == (d_curr + 1)) {

    post <- calcPostRatio(gamma_prop = gamma_prop, gamma_curr = gamma_curr,
                          var_prop = var_prop, var_curr = var_curr,
                          Q_full = Q_full, d_prop = d_prop, d_curr = d_curr,
                          y_vec = y_vec, dim_prior_vec = dim_prior_vec,
                          R_inv_full = R_inv_full)

    prop <- calcPropRatio(gamma_k_prime_prop = gamma_prop[d_prop + 1],
                          k_prime_init = full_inits[d_prop + 1],
                          innov_sigma = innov_sigma, d_prop = d_prop, d_curr = d_curr,
                          d_min = d_min, d_max = d_max, var_prop = var_prop,
                          var_curr = var_curr, innov_sd_var = innov_sd_var)
    accept_prob <- min(1, post * prop)

  } else if (d_prop == (d_curr - 1)) {
    # different here, have to ~flip it n reverse it~
    # basically have to pretend proposal is current, and current is proposal
    # then invert after calculating accept_prob

    post <- calcPostRatio(gamma_prop = gamma_curr,
                          gamma_curr = gamma_prop,
                          var_prop = var_curr,
                          var_curr = var_prop,
                          Q_full = Q_full,
                          d_prop = d_curr,
                          d_curr = d_prop,
                          y_vec = y_vec, dim_prior_vec = dim_prior_vec,
                          R_inv_full = R_inv_full)

    prop <- calcPropRatio(gamma_k_prime_prop = gamma_curr[d_curr + 1],
                          k_prime_init = full_inits[d_curr + 1],
                          innov_sigma = innov_sigma,
                          d_prop = d_curr,
                          d_curr = d_prop,
                          d_min = d_min, d_max = d_max, var_prop = var_curr,
                          var_curr = var_prop, innov_sd_var = innov_sd_var)
    temp <- (post * prop)^(-1)

    # Here, we have to check after we computed the post ratio to see if the
    # proposal was monotonic, otherwise it can accept non-monotonic proposals
    # (as it is only checking gamma_curr in calcPostRatio technically.)

    beta_prop <- gammaToBeta(gamma = gamma_prop,
                             R_inv = R_inv_full[1:(length(gamma_prop)), 1:(length(gamma_prop))])
    if (!MonoPoly::ismonotone(beta_prop, a = 0, b = 1)) {
      temp <- 0
    }


    accept_prob <- min(1, temp)

  }


  return(accept_prob)

}