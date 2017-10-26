# some plotting functions

#' Title
#'
#' @param obj object of type "rjmonopolyfit" for plotting
#' @param ... other optional parameters for controlling plot style
#' @param only if specified, plot only the barplot
#'
#' @return obj
#' @export
#'
plotDegreePost <- function(obj, only = "barplot", ...) {
  if (!class(obj) == "rjmonopolyfit") {
    print("Incorrect object type")
    return(NULL)

  }

  samples <- obj$d_samples

  if (missing(only)) {
    par(mfrow = c(2,1))
    barplot(table(samples) / length(samples), xlab = "Polynomial degree", ylab = "Posterior probability")
    plot(samples, xlab = "Sample index", ylab = "Polynomial degree", type = "l")
    # dev.off()
    par(mfrow = c(1,1))

  } else if (only == "barplot") {
    barplot(table(samples) / length(samples), xlab = "Polynomial degree", ylab = "Posterior probability")
  }

  invisible(obj)

}

#' Plot of fitted curve, with credible interval. Either uses the posterior
#' mode degree, or degree closest to posetrior mean.
#'
#' @param obj object of type
#' @param mode which degree to plot, "mean" for rounded mean, "mode" for
#' posterior mode
#' @param degree (optional) Specify a specific degree to plot
#' @param ...
#'
#' @return obj
#' @export
plotFit <- function(obj, mode = "mean", degree = NA, ...) {
  if (!class(obj) == "rjmonopolyfit") {
    print("Incorrect object type")
    return(NULL)
  }

  plot_deg <- NA

  if (mode == "mean") {
    plot_deg <- round(mean(obj$d_samples))

  } else if (mode == "mode") {
    degs <- unique(obj$d_samples)
    plot_deg <- degs[which.max(tabulate(match(obj$d_samples, degs)))]
  }

  if (!missing(degree)) {
    # this could be a lot more robust.
    plot_deg <- degree
  }

  gamma_samples <- obj$gamma_samples
  gamma_samples <- rjmonopoly:::returnSamplesOfLengthN(gamma_samples, (plot_deg + 1))

  gamma_mean <- apply(gamma_samples, 2, mean)

  fitted_post <- obj$Q_full[, 1:(plot_deg + 1)] %*% t(gamma_samples)
  fitted_lower <- apply(fitted_post, 1, function(x){quantile(x, 0.025)})
  fitted_upper <- apply(fitted_post, 1, function(x){quantile(x, 0.975)})
  fitted_mean <- obj$Q_full[, 1:(plot_deg + 1)] %*% gamma_mean

  # how to plot using fitted values, but with original scale? I solved this once
  # somewhere else.

  plot_x <- rjmonopoly:::.unscale(obj$x, obj$scaling_factors$x)

  plot_y <- rjmonopoly:::.unscale(obj$y, obj$scaling_factors$y)
  plot_mean <- rjmonopoly:::.unscale(fitted_mean, obj$scaling_factors$y)
  plot_upper <- rjmonopoly:::.unscale(fitted_upper, obj$scaling_factors$y)
  plot_lower <- rjmonopoly:::.unscale(fitted_lower, obj$scaling_factors$y)

  plot_df <- data.frame(x = plot_x, y = plot_y, fitted = plot_mean,
                        lower = plot_lower, upper = plot_upper)

  plot_obj <- ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = fitted), colour = "red") +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(x = x, ymin = lower, ymax = upper), alpha = 0.25) +
    ggplot2::theme_light()

  print(plot_obj)
  invisible(plot_obj)
}