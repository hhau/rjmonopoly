#' Thin an rjmonopolyfit
#'
#' Thins an object of type "rjmonopolyfit" by *thin_period*. All quantites are
#' thinned by the same period.
#'
#' @param obj object of type "rjmonopolyfit"
#' @param thin_period period of thinned samples
#' @param ... other
#' @param warm_up number of samples that will be removed from the beggining of
#' the original object
#'
#' @return an object of type rjmonopolyfit, with n_accepted changed to
#' n_accepted_orig(inally)
#' @export
#'
thin <- function(obj, thin_period, warm_up = 1, ...) {

  if (class(obj) != "rjmonopolyfit") {
    stop("Object not of type rjmonopolyfit")
  }

  n_samples <- length(obj$d_samples)
  thin_vector <- seq(from = warm_up, to = n_samples, by = thin_period)

  res <- list(gamma_samples = obj$gamma_samples[thin_vector],
              var_samples = obj$var_samples[thin_vector],
              d_samples = obj$d_samples[thin_vector],
              n_accept_orig = obj$n_accept,
              scaling_factors = obj$scaling_factors,
              x = obj$x,
              y = obj$y,
              Q_full = obj$Q_full,
              R_inv_full = obj$R_inv_full)
  class(res) <- "rjmonopolyfit"
  return(res)
}
