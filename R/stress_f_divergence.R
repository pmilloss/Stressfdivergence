#' stress_mean_div
#'
#' Provides weights on simulated scenarios from a baseline stochastic model, such that a stressed model component (random variable) fulfill a moment constraint. Scenario weights are selected by constrained minimisation of a selected divergence to baseline model.
#'
#' @param f A function that, applied to \code{x}, constitute the moment constraints.
#'
#' @param k A vector indicating which columns of \code{x} the function \code{f} operates on. By default, the first columnn is selected and the identity function is applied.
#'
#' @param m Numeric, the stressed moments of \code{f(x)}. Must be in the
#' range of \code{f(x)}.
#'
#' @param div Character. One of "Chi2", "KL", "Hellinger", "Alpha", "Triangular", "Jeffrey" or "user". When a user specified divergence is chosen, the additional parameters "inv.div" (inverse of the divergence function) and "d.div" (derivative of the divergence function) must be passed. For the "Alpha" divergence, the numeric parameter "alpha" must be provided (when alpha is 1 the KL divergence is used)
#'
#' @param inv.div A function specifying the inverse of the divergence function, when a user divergence function is used.
#'
#' @param d.div A function specifying the derivative of the divergence function, when a user divergence function is used.
#'
#' @param d.inv An optional function specifying the derivative of the inverse of the divergence function, when a user divergence function is used.
#'
#' @param p a set of optional weights specifying the baseline model.
#'
#' @param normalise Logical. If true, values of \code{f(x)} are linearly scaled to the unit interval.
#'
#' @param show Logical. If true, print the result of the call to \code{\link[nleqslv]{nleqslv}}.
#'
#' @param start A numeric vector with two elements, the starting values for the+++++++
#'
#' @param sumRN Logical. If true, the++++++++++
#'
#' @param ... Additional arguments to be passed to \code{\link[nleqslv]{nleqslv}}.
#'
#' @seealso See \code{\link{stress_moment}} for a more flexible function allowing to perform multiple joint stresses under the KL divergnce.
#'
#' @return A named list
#' @export

stress_mean_div <- function(x, f = function(x)x, k = 1, m, div = c("Chi2", "KL", "Hellinger", "Alpha", "Triangular", "Jeffrey", "user"), inv.div = NULL, d.div = NULL, d.inv = NULL, p = rep(1 / length(x), length(x)), alpha = NULL, normalise = TRUE, show = FALSE, start = NULL, suMRN = FALSE, ...){

  if (is.SWIM(x)) x_data <- SWIM::get_data(x) else x_data <- as.matrix(x)
  if (anyNA(x_data)) warning("x contains NA")
  if (!is.function(f)) stop("f must be a function")
  if (!is.numeric(k)) stop("k must be a numeric vector")
  if (!is.numeric(m)) stop("m must be numeric")

  if (!is.numeric(p) | any(p < 0) | anyNA(p)) stop("p must be a vector of nonnegative weights") else p <- p / sum(p)

  z <- f(x_data[, k])
  min.fz <- min(z)
  max.fz <- max(z)
  if (m < min.fz || m > max.fz) stop("Values in m must be in the range of f(x)")

  if (normalise == TRUE) {
    z <- apply(z, 2, .scale)
    m <- (m - min.fz) / (max.fz - min.fz)
    }

  if (div == "user" & (is.null(inv.div) | is.null(d.inv))) stop("For a user defined divergence, the arguments 'inv.div' and 'd.div' must be provided")

  if (div == "Alpha" & (is.null(alpha) | !is.numeric(alpha))) stop("For the Alpha divergence, the numeric argument 'alpha' must be provided")

  if (div == "Chi2" | (div == "Alpha" & alpha == 2)) {
    inv.div <- function(x)0.5 * x
    d.div <- function(x)2 * x
    d.div0 <- d.div(0)
    d.inv <- function(x)0.5
  } else if (div == "KL") {
    inv.div <- function(x)exp(x - 1)
    d.div <- function(x)1 + log(x)
    d.div0 <- -Inf
    d.inv <- function(x)inv.div(x)
  } else if (div == "Hellinger") {
    inv.div <- function(x)1 / (1 - x) ^ 2
    d.div <- function(x)(1 - 1 / sqrt(x))
    d.div0 <- -Inf
    d.inv <- function(x)2 / (1 - x) ^ 3
  } else if (div == "Alpha") {
    inv.div <- function(x)(-2 * x * (1 + alpha)) ^ (-2 / (1 + alpha))
    d.div <- function(x)(-0.5 * (x ^ (-0.5 * (1 + alpha))) / (1 + alpha))
    d.div0 <- d.div(0)
    d.inv <- function(x)4 * (-2 * x * (1 + alpha)) ^ (- (3 + alpha) / (1 + alpha))
  } else if (div == "Triangular") {
    inv.div <- function(x)-1 + 2 / sqrt(ifelse(x < 1, 1 - x, 0))
      # ifelse(x < 1, -1 + 2 / sqrt(1 - x), Inf)
    d.div <- function(x)(x - 1) * (x + 3) / (x + 1) ^ 2
    d.div0 <- -3
    # d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
    d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
  } else if (div == "Jeffrey") {
  } else if (div == "user") {
    d.div0 <- d.div(0)
  } else stop("The argument 'div' must be one of 'Chi2', 'KL', 'Hellinger', 'Alpha', 'Triangular', 'Jeffrey' or 'user' ")

  constraints <- function(L){
    RN <- g(pmax(d.div0, L[1] + L[2] * x))
    F1 <- sum(p * RN) - 1
    F2 <- sum(p * RN * x) - new_mean
    return(c(F1, F2))
  }

  if (is.null(start)) start <- c(L1 = d.div(1), L2 = 0)

  if (!is.null(d.inv)) {
    J <- function(L) {
      RN.der <- d.inv(L[1] + L[2] * x)
      ind <- (L[1] + L[2] * x > d.div0)
      J11 <- sum(p * RN.der * ind)
      J12 <- sum(p * RN.der * x * ind)
      J22 <- sum(p * RN.der * x * x * ind)
      mat <- matrix(c(J11, J12, J12, J22), nrow = 2, byrow = TRUE)
      return(mat)
    }
    sol <- nleqslv::nleqslv(x = start, fn = constraints, jac = J, ...)
  } else {
    sol <- nleqslv::nleqslv(x = start, fn = constraints, ...)
    }

  if (sol$termcd != 1) warning(paste("nleqslv terminated with code ", sol$termcd))

  constr_moment <- list("k" = k, "m" = m, "f" = f)
  constr <- list(constr_moment)

  # Name stresses
  if (is.null(names)) {
    temp <- paste("stress", 1)
  } else {
    temp <- names
  }
  if (length(temp) != 1) stop("length of names are not the same as the number of models")
  names(constr) <- temp


  if (is.null(colnames(x_data))) colnames(x_data) <-  paste("X", 1:ncol(x_data), sep = "")

  w <- inv.div(pmax(d.div0, sol$x[1] + sol$x[2] * z))
  if(sumRN) w <- w / mean(w)
  new_weights = list()
  new_weights[[temp]] <- w

  type <- list("moment") #++++++++++++++++++++
  my_list <- SWIM("x" = x_data, "new_weights" = new_weights, "type" = type, "specs" = constr)
  if (is.SWIM(x)) my_list <- merge(x, my_list)
  if (show == TRUE) print(sol)
  m.ac <- colMeans(z * as.vector(exp(z %*% sol$x)))[-1]
  if (normalise == TRUE){
    m <- min.fz + (max.fz - min.fz) * m
    m.ac <- min.fz + (max.fz - min.fz) * m.ac
  }
  err <- m - m.ac
  rel.err <- (err / m) * (m != 0)
  outcome <- data.frame(cols = as.character(k), required_moment = m, achieved_moment = m.ac, abs_error = err, rel_error = rel.err)
  print(outcome)

  if (log) {
    summary_weights(my_list)
  }

  return(my_list)

}

set.seed(0)
n <- 10 ^ 5
x <- rlnorm(n, meanlog = 1, sdlog = 1)
mean(x)
stress <- 1.2
t <- stress * mean(x)

res_KL <- stress_mean_div(x = x, new_mean = t, div = "KL", trace = T)
plot(x, res_KL)

res_Chi2 <- stress_mean_div(x = x, new_mean = t, div = "Chi2", trace = T)
plot(x, res_Chi2)

res_Hell <- stress_mean_div(x = x, new_mean = t, div = "Hellinger", trace = T)# fails after inspection
res_Hell1 <- stress_mean_div(x = x, new_mean = t, div = "Hellinger", trace = T, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 1000)) # works better, changing the scaling
plot(x, res_Hell1, xlim = c(0, 150))
plot(x, res_Hell1, xlim = c(0, 150), ylim = c(0, 10))
lines(x, res_Hell, xlim = c(0, 150), ylim = c(0, 10))# madness

res_alpha <- stress_mean_div(x = x, new_mean = t, div = "Alpha", alpha = -2, trace = T)
plot(x, res_alpha)

res_alpha2 <- stress_mean_div(x = x, new_mean = t, div = "Alpha", alpha = 0.5, trace = T) #ok
res_alpha3 <- stress_mean_div(x = x, new_mean = t, div = "Alpha", alpha = 1.5, trace = T) #fails
plot(x, res_alpha)
plot(x, res_alpha2)

res_trian <- stress_mean_div(x = x, new_mean = t, div = "Triangular", trace = T)# fails
res_trian <- stress_mean_div(x = x, new_mean = t, div = "Triangular", trace = T, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 1000))# this works better
