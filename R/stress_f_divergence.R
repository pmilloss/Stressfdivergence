#' stress_mean_div
#'
#' Provides weights on simulated scenarios from a baseline stochastic model, such that (i) a stressed model component fulfills a moment constraint and scenario weights are selected by minimisation of a selected divergence to the baseline model, or (ii) a stressed model fulfills a divergence constraint and scenario weights are selected by maximisation of the moment of a model component. Case (i) is obtained by specifying the stressed moment parameter 'm', case (ii) by specifying the divergence constraint 'theta'.
#'
#' @param x A vector, matrix or data frame containing realisations of random variables. Columns of \code{x} correspond to random variables; or a \code{SWIM} object, where \code{x} corresponds to the underlying data of the \code{SWIM} object.
#'
#' @param f A function that, applied to \code{x}, generates the valued of the variable that will satisfy a mean constraint (and the divergence will be minimized) or whose mean will be maximized (under a divergence constraint). By default it is the identity function.
#'
#' @param k A vector indicating which columns of \code{x} the function \code{f} operates on. By default/, the first columnn is selected.
#'
#' @param m Numeric, the stressed moments of \code{f(x)}. Must be in the
#' range of \code{f(x)}.
#'
#' @param theta Numeric, the stressed divergence of the model.  range of \code{f(x)}. +++++++++++++++++
#'
#' @param dvg Character. One of "Chi2", "KL", "Hellinger", "Alpha", "Triangular", "Jeffrey", "LeCam" or "user". When a user specified divergence is chosen, the additional list 'div' containing ++++++++parameters "inv" (inverse of the first derivative of the divergence function) and "d.div0" (derivative of the divergence function) must be passed. For the "Alpha" divergence, the numeric parameter "alpha" must be provided (when alpha is 1 or 2 the KL, respectively the Chi2, divergence is used). +++++++++++++++++++ADD EQUATIONS FOR THE DIVERGENCES+++++++++++++++++
#'
#' @param div.usr When a user divergence function is used, a named list containing at least the following objects: a function 'inv' giving the inverse of the first derivative of the divergence; a numeric 'div0' giving the the derivative at 0 of the divergence (possibly -Inf); when 'theta' is specified (divergence constraint), a function 'div' giving the divergence. Optionally, a function 'd.div' specifying the first derivative of the divergence.
#'
#' @param p Numeric. A set of optional nonnegative scenario weights specifying the baseline model.
#'
#' @param alpha Numeric. The 'alpha' parameter when the Alpha divergence is used.
#'
#' @param normalise Logical. If true, values of \code{f(x)} are linearly scaled to the unit interval. Not used when 'theta' is specified.
#'
#' @param show Logical. If true, print the result of the call to \code{\link[nleqslv]{nleqslv}}.
#'
#' @param names Character vector, the names of stressed models.
#'
#' @param start A numeric vector with two elements, the starting values for the coefficients lambda1 and lambda2. Defaults to div$d.div(1) and 0 respectively, guaranteeing that the initial set of scenario weights is constant.
#'
#' @param sumRN Logical. If true, the scenario weights are normalized so as to average to 1 exactly.
#'
#' @param log Logical, the option to print weights' statistics.
#'
#' @param ... Additional arguments to be passed to \code{\link[nleqslv]{nleqslv}}.
#'
#' @seealso See \code{\link{stress_moment}} for a more flexible function allowing to perform multiple joint stresses under the KL divergence.
#'
#' @author Pietro Millossovich
#'
#' @return A named list
#' @export

stress_mean_div <- function(x, f = function(x)x, k = 1, m = NULL, theta = NULL, dvg = c("Chi2", "KL", "Hellinger", "Alpha", "Triangular", "Jeffrey", "LeCam", "user"), div.usr = NULL, p = rep(1 / length(x), length(x)), alpha = NULL, normalise = TRUE, show = FALSE, names = NULL, start = NULL, sumRN = FALSE, log = FALSE, ...){

  min.d <- !is.null(m)
  max.l <- !is.null(theta)

  if (min.d + max.l != 1) stop("exactly one of m and theta must be provided")
  if (SWIM:::is.SWIM(x)) x_data <- SWIM::get_data(x) else x_data <- as.matrix(x)
  if (anyNA(x_data)) warning("x contains NA")
  if (!is.function(f)) stop("f must be a function")
  if (!is.numeric(k)) stop("k must be a numeric vector")
  if (min.d & !is.numeric(m)) stop("m must be numeric")
  if (max.l & !is.numeric(theta)) stop("theta must be numeric")

  if (!is.numeric(p) | any(p < 0) | anyNA(p)) stop("p must be a vector of nonnegative weights") else p <- p / sum(p)

  if (dvg == "user" & (!is.function(div.usr$inv) | !is.numeric(div.usr$div0))) stop("For a user defined divergence, the list 'div.usr' must contain at least the objects 'inv' and 'div0' must be provided")

  if (dvg == "user" & max.l & !is.function(div.usr$div)) stop("For a user defined divergence, when 'theta' is specified, the list 'div.usr' must also contain the object 'div'")

  if (dvg == "Alpha" & (is.null(alpha) | !is.numeric(alpha))) stop("For the Alpha divergence, the numeric argument 'alpha' must be provided")

  if (dvg == "Chi2" | (dvg == "Alpha" & identical(alpha, 2))) {
    div <- function(x)x ^ 2 - 1
    d.div <- function(x)2 * x
    inv <- function(x)0.5 * x
    div0 <- d.div(0)
    d.inv <- function(x)0.5
  } else if (dvg == "KL" | (dvg == "Alpha" & identical(alpha, 1))) {
    div <- function(x)ifelse(x > 0, x * log(x), 0)
    inv <- function(x)exp(x - 1)
    d.div <- function(x)1 + log(x)
    div0 <- -Inf
    d.inv <- function(x)inv(x)
  } else if (dvg == "Hellinger") {
    div <- function(x)(sqrt(x) - 1) ^ 2
    inv <- function(x)1 / (1 - x) ^ 2
    d.div <- function(x)1 - 1 / sqrt(x)
    div0 <- -Inf
    d.inv <- function(x)2 / (1 - x) ^ 3
  } else if (dvg == "Alpha") {
    div <- function(x)(x ^ alpha - alpha * (x - 1) - 1) / (alpha * (alpha - 1))
    inv <- function(x)(-2 * x * (1 + alpha)) ^ (-2 / (1 + alpha))
    d.div <- function(x)(-0.5 * (x ^ (-0.5 * (1 + alpha))) / (1 + alpha))
    div0 <- d.div(0)
    d.inv <- function(x)4 * (-2 * x * (1 + alpha)) ^ (- (3 + alpha) / (1 + alpha))
    dvg <- paste(dvg, " a=", alpha)
  } else if (dvg == "Triangular") {
    div <- function(x)(x - 1) ^ 2 / (x + 1)
    inv <- function(x)-1 + 2 / sqrt(ifelse(x < 1, 1 - x, 0))
      # ifelse(x < 1, -1 + 2 / sqrt(1 - x), Inf)
    d.div <- function(x)(x - 1) * (x + 3) / (x + 1) ^ 2
    div0 <- -3
    # d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
    d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
  } else if (dvg == "Jeffrey") {
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++

  } else if (dvg == "LeCam") {
    div <- function(x)0.5 * (1 - x) / (x + 1)
  } else if (dvg == "user") {
    div0 <- div.usr$div0
    inv <- div.usr$inv
    if (max.l) div <- div.usr$div
  } else stop("The argument 'div' must be one of 'Chi2', 'KL', 'Hellinger', 'Alpha', 'Triangular', 'Jeffrey', 'LeCam' or 'user' ")

  z <- f(x_data[, k])
  min.fz <- min(z)
  max.fz <- max(z)
  if (min.d & (m < min.fz | m > max.fz)) stop("m must be in the range of f(x)")
  # add a check if there are too few points above m as in stress VaR????

  pn <- p[which.max(z)]
  max.div <- div(1 / pn) * pn
  if (max.l & (theta < 0 | theta > max.div)) stop("theta must be in the range [0,max.div], where max.div=div(pn)/pn and pn is the weight corresponding to the largest observation") # ++++++++++++++++

  if (min.d & normalise == TRUE) {
    z <- (z - min.fz) / (max.fz - min.fz)
    m <- (m - min.fz) / (max.fz - min.fz)
  }

  if (min.d) constr <- function(L){
    RN <- inv(pmax(div0, L[1] + L[2] * z))
    C1 <- sum(p * RN) - 1
    C2 <- sum(p * RN * z) - m
    return(c(C1, C2))
  } else constr <- function(L){
    RN <- inv(pmax(div0, L[1] + L[2] * z))
    C1 <- sum(p * RN) - 1
    C2 <- sum(p * div(RN)) - theta
    return(c(C1, C2))
  }

  if (is.null(start) & dvg != "user") start <- c(L1 = d.div(1), L2 = 0)


  if (exists(d.inv) & min.d) {
    J <- function(L) {
      d.RN <- d.inv(L[1] + L[2] * z)
      ind <- (L[1] + L[2] * z > div0)
      J11 <- sum(p * d.RN * ind)
      J12 <- sum(p * d.RN * z * ind)
      J22 <- sum(p * d.RN * z * z * ind)
      jac <- matrix(c(J11, J12, J12, J22), nrow = 2)
      return(jac)
    } else if (exists(d.inv) & max.l) {
      J <- function(L) {
        d.RN <- d.inv(L[1] + L[2] * z)
        ind <- (L[1] + L[2] * z > div0)
        J11 <- sum(p * d.RN * ind)
        J12 <- sum(p * d.RN * z * ind)
        J22 <- sum(p * d.RN * z * z * ind)
        J21 <- L[1] * J11 + L[2] * J12
        J22 <- L[1] * J12 + L[2] * J22
        jac <- matrix(c(J11, J21, J12, J22), nrow = 2)
        return(jac)
    }
    sol <- nleqslv::nleqslv(x = start, fn = constr, jac = J, ...)
  } else {
    sol <- nleqslv::nleqslv(x = start, fn = constr, ...)
    }

  if (sol$termcd != 1) warning(paste("nleqslv terminated with code ", sol$termcd))

  w <- inv(pmax(div0, sol$x[1] + sol$x[2] * z))
  if(sumRN) w <- w / mean(w)

  if (min.d) {
    m.ac <- sum(p * w * z)
    if(normalise == TRUE) {
      m <- min.fz + (max.fz - min.fz) * m
      m.ac <- min.fz + (max.fz - min.fz) * m.ac
    }
    err <- m - m.ac
    rel.err <- (err / m) * (m != 0)
    outcome <- c(required_moment = m, achieved_moment = m.ac, abs_error = err, rel_error = rel.err)
    constr_moment <- list("k" = k, "m" = m, "f" = f, "div" = dvg)
  } else {
    theta.ac <- sum(p * div(w))
    err <- theta - theta.ac
    rel.err <- (err / theta) * (theta != 0)
    outcome <- c(required_divergence = theta, achieved_divergence = theta.ac, abs_error = err, rel_error = rel.err)
    constr_moment <- list("k" = k, "theta" = theta, "f" = f, "div" = dvg)
  }

  print(outcome)
  constr <- list(constr_moment)

  # Name stresses +++++++++++++++++++++++++
  if (is.null(names)) {
    temp <- paste("stress", 1)
  } else {
    temp <- names
  }
  if (length(temp) != 1) stop("length of names are not the same as the number of models")
  names(constr) <- temp

  if (is.null(colnames(x_data))) colnames(x_data) <-  paste("X", 1:ncol(x_data), sep = "")

  new_weights = list()
  new_weights[[temp]] <- w

  type <- list("moment")
  my_list <- SWIM:::SWIM("x" = x_data, "new_weights" = new_weights, "type" = type, "specs" = constr)
  if (SWIM:::is.SWIM(x)) my_list <- merge(x, my_list)

  if (show == TRUE) print(sol)

  if (log) {
    SWIM::summary_weights(my_list)
  }

  return(my_list)

}
