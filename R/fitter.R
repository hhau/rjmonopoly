# fit function

#' Reversible jump sampler for monotonic polynomials
#'
#' Other notes
#'
#' Note that the samples for Gamma and Sigma^2 are on the orthonormal scale.
#'
#' @param x vector of x values
#' @param y vector of y values
#' @param d_min minimum allowable polynomial degree
#' @param d_max maximum allowable polynomial degree
#' @param iter number of iterations
#' @param prior_prob probability for binomial prior
#' @param starting_var_val starting value for variance
#' @param control (optional) list for innovation variances
#' @param progress turn the progress bar on / off, TRUE is default, and shows
#' progress bar.
#' @param prior_option if set to string "flat", then will use a flat prior over
#' the allowable parameter space.
#'
#' @return object of type rjmonopol_fit containing lists of samples for all
#' regression parameters
#' @export
#'
#' @examples
#' library(fda)
#' x <- onechild$day
#' y <- onechild$height
#' fit <- rjmonopoly(x, y)
rjmonopoly <- function(
    x,
    y,
    d_min = 2,
    d_max = 10,
    iter = 50000,
    prior_option = NA,
    prior_prob = 0.5,
    starting_var_val = 0.001,
    control = list(
      innov_sd_beta = 0.025,
      innov_sd_var = 0.025
    ),
    progress = TRUE
  ) {

  if (!missing(control)) {
    control <- .extend_list(formals(rjmonopoly)$control, control)
  }

  # rescale inputs
  x_rescl <- (x - min(x)) / (diff(range(x)))
  y_rescl <- (y - min(y)) / (diff(range(y)))

  temp_df <- data.frame(x_rescl = x_rescl, y_rescl = y_rescl)
  temp_df <- temp_df[order(temp_df$x_rescl),]

  x_rescl <- temp_df$x_rescl
  y_rescl <- temp_df$y_rescl

  # define the prior distribution
  prior_vec <- genDimPrior(d_min, d_max, prior_prob)

  if (!missing(prior_option)) {
    if (prior_option == "flat") {
      # you need to pad this vector with zeros at the front, so that it has length
      # d_max. This is because elsewhere you use d / d_prop to get the dimension
      # prior probability.
      prior_vec <- c(rep(0, d_min - 1),
                     rep(1/length(d_min:d_max), length(d_min:d_max)))

    } else {
        stop("Unknown prior option, stopping")

      }
  }

  d_init <- round(median(d_min:d_max))
  beta_init_full <- coef(MonoPoly::monpol(y_rescl ~ x_rescl, degree = d_max, a = 0, b = 1))

  qr_mats <- genAllMatrices(x_rescl, d_max)
  Q_full <- qr_mats[[1]]
  R_inv_full <- qr_mats[[2]]

  gamma_init_full <- betaToGamma(beta_init_full, R_inv_full)

  gamma_samples <- list()
  gamma_samples[[1]] <- gamma_init_full[1:(d_init + 1)]
  innov_sd_beta <- control$innov_sd_beta

  var_samples <- matrix(NA, nrow = iter + 1, ncol = 1)
  var_samples[1] <- starting_var_val
  innov_sd_var <- control$innov_sd_var

  d_samples <- matrix(NA, nrow = iter + 1, ncol = 1)
  d_samples[1] <- d_init

  n_accept <- 0

  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }
  for (ii in 2:(iter + 1)) {
    d_prop <- dimProposer(d_current = d_samples[ii - 1],
                                       d_min = d_min,
                                       d_max = d_max)

    var_prop <- genVarProp(var_prev = var_samples[ii - 1],
                                        innov_sd_var = innov_sd_var)

    gamma_prop <- genGammaProp(current_gamma = gamma_samples[[ii - 1]],
                                             d_prop = d_prop,
                                             d_curr = d_samples[ii - 1],
                                             innov_sigma = innov_sd_beta,
                                             full_inits = gamma_init_full)

    accept_probability <- calcAcceptProb(gamma_prop = gamma_prop,
                                                      gamma_curr = gamma_samples[[ii - 1]],
                                                      d_prop = d_prop,
                                                      d_curr = d_samples[ii - 1],
                                                      d_min = d_min,
                                                      d_max = d_max,
                                                      var_prop = var_prop,
                                                      var_curr = var_samples[ii - 1],
                                                      Q_full = Q_full,
                                                      full_inits = gamma_init_full,
                                                      y_vec = y_rescl,
                                                      innov_sigma = innov_sd_beta,
                                                      dim_prior_vec = prior_vec,
                                                      innov_sd_var = innov_sd_var,
                                                      R_inv_full = R_inv_full)

    if (runif(1) < accept_probability) {
      gamma_samples[[ii]] <- gamma_prop
      d_samples[ii] <- d_prop
      var_samples[ii] <- var_prop
      n_accept <- n_accept + 1

    } else {
      gamma_samples[[ii]] <- gamma_samples[[ii - 1]]
      d_samples[ii] <- d_samples[ii - 1]
      var_samples[ii] <- var_samples[ii - 1]

    }
    if (progress) {
      utils::setTxtProgressBar(pb = pb, value = ii)
    }
  }

  res <- list(gamma_samples = gamma_samples,
              var_samples = var_samples,
              d_samples = d_samples,
              n_accepted = n_accept,
              scaling_factors = list(
                y = list(min = min(y),
                         max = max(y),
                         scl_range = diff(range(y))),
                x = list(min = min(x),
                         max = max(x),
                         scl_range = diff(range(x)))
              ),
              x = x_rescl,
              y = y_rescl,
              Q_full = Q_full,
              R_inv_full = R_inv_full
          )

  class(res) <- "rjmonopolyfit"
  return(res)

}