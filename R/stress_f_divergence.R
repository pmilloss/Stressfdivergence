#' stress_mean_div
#'
#' Provides weights on simulated scenarios from a baseline stochastic model, such that (i) a stressed model component fulfills a moment constraint and scenario weights minimize a selected divergence to the baseline model, or (ii) a stressed model fulfills a divergence constraint and scenario weights maximize the moment of a model component. Case (i) is obtained by specifying the moment parameter '\code{m}', case (ii) by specifying the divergence constraint '\code{theta}'.
#'
#' @param x A vector, matrix or data frame containing realizations of random variables. Columns of \code{x} correspond to random variables; OR \cr
#' A \code{SWIM} object, where \code{x} corresponds to the underlying data of the \code{SWIM} object.
#' @param f A function that, applied to columns of \code{x}, generates the values of the variable that will satisfy a mean constraint (and the divergence will be minimized) or whose mean will be maximized (under a divergence constraint). By default it is the identity function.
#' @param k A vector indicating which columns of \code{x} the function \code{f} operates on. By default, the first column is selected.
#' @param m Numeric, the stressed moments of \code{f(x)}. Must be in the range of \code{f(x)}.
#' @param theta Numeric, the stressed divergence of the model. The range of possible divergence values goes from 0 to \code{max.div=div(0) * (1 - 1 / n) + div(n) / n}, where \code{div} is the divergence function and \code{n} is the number of simulations (rows of \code{x}).
#' @param dvg Character. One of "\code{Chi2}", "\code{KL}", "\code{revKL}", "\code{Hellinger}", "\code{Alpha}", "\code{Triangular}", "\code{LeCam}" or "\code{user}". For a user specified divergence, see the additional list 'div.usr' that must be passed. For the "\code{Alpha}" divergence, the numeric parameter "\code{alpha}" must be provided (when \code{alpha=1} the KL divergence is used; when alpha is 2, the Chi2 divergence is used). +++++++++++++++++++ADD EQUATIONS FOR THE DIVERGENCES+++++++++++++++++
#' @param div.usr When a user divergence function is chose, a named list containing at least the following objects: a function 'inv' giving the inverse of the first derivative of the divergence; a numeric 'div0' giving the the derivative at 0 of the divergence (possibly -Inf); when 'theta' is specified (divergence constraint), a function 'div' giving the divergence. Optionally, a function 'd.div' specifying the first derivative of the divergence. The divergence function must be one for which ###########
#' @param alpha Numeric. The 'alpha' parameter when the Alpha divergence is used.
#' @param normalise Logical. If true, values of \code{f(x)} are linearly scaled to the unit interval. Not used when 'theta' is specified.
#' @param show Logical. If true, print the result of the call to \code{\link[nleqslv]{nleqslv}}.
#' @param names Character vector, the names of stressed models.
#' @param start A numeric vector with two elements, the starting values for the coefficients lambda1 and lambda2. Defaults to div$d.div(1) and 0 respectively, guaranteeing that the initial set of scenario weights is constant.
#' @param sumRN Logical. If true, the scenario weights are normalized so as to average to 1.
#' @param log Logical, the option to print weights' statistics.
#' @param use.jac Logical, should the Jacobian matrix be used in the call to \code{\link[nleqslv]{nleqslv}}. Set to TRUE when a user divergence is used and the argument div.usr$d.div is provided.
#' @param ... Additional arguments to be passed to \code{\link[nleqslv]{nleqslv}}.
#'
#' @seealso See \code{\link{stress_moment}} for a more flexible function allowing to perform multiple joint stresses under the KL divergence.
#'
#' @details The moment constraints are given by \code{E^Q( f(x) ) = m},where \code{E^Q} denotes the expectation under the stressed model. \code{stress_moment} solves the subsequent set of equations with respect to theta, using \code{\link[nleqslv]{nleqslv}} from package \code{\link[nleqslv]{nleqslv}}:
#'
#'  \deqn{E^Q( f(x) ) = E( f(x) * exp(theta * f(x)) ) = m.}
#'
#' There is no guarantee that the set of equations has a solution, or that the solution is unique. \code{SWIM} will
#'     return a warning if the termination code provided by \code{nleqslv} is
#'     different from 1 (convergence has been achieved). It is recommended to
#'     check the result of the call to \code{nleqslv} using the "show" argument. The
#'     user is referred to the \code{\link[nleqslv]{nleqslv}} documentation for
#'     further details.
#'
#'
#' @author Pietro Millossovich
#'
#' @examples
#' # example code
#'
#'
#' @return A named list
#'
#' @export

stress_mean_div <- function(x, f = function(x)x, k = 1, m = NULL, theta = NULL, dvg = c("Chi2", "KL", "revKL", "Hellinger", "Alpha", "Triangular", "Jeffrey", "LeCam", "user"), div.usr = NULL, alpha = NULL, normalise = TRUE, show = FALSE, names = NULL, start = NULL, sumRN = FALSE, log = FALSE, use.jac = FALSE, ...){

  min.d <- !is.null(m)
  max.l <- !is.null(theta)

  if (min.d + max.l != 1) stop("exactly one of m and theta must be provided")
  if (SWIM:::is.SWIM(x)) x_data <- SWIM::get_data(x) else x_data <- as.matrix(x)
  if (anyNA(x_data)) warning("x contains NA")
  if (!is.function(f)) stop("f must be a function")
  if (!is.numeric(k)) stop("k must be a numeric vector")
  if (min.d & !is.numeric(m)) stop("m must be numeric")
  if (max.l & !is.numeric(theta)) stop("theta must be numeric")

  n <- nrow(x_data)

  if (dvg == "user" & (!is.function(div.usr$inv) | !is.numeric(div.usr$div0))) stop("For a user defined divergence, the list 'div.usr' must contain at least the objects 'inv' and 'div0'")

  if (dvg == "user" & max.l & !is.function(div.usr$div)) stop("For a user defined divergence, when 'theta' is specified, the list 'div.usr' must also contain the object 'div'")

  if (dvg == "Alpha" & (is.null(alpha) | !is.numeric(alpha))) stop("For the Alpha divergence, the numeric argument 'alpha' must be provided")

  if (dvg == "Chi2" | (dvg == "Alpha" & identical(alpha, 2))) {
    div <- function(x)ifelse(x >= 0, x ^ 2 - 1, Inf)
    d.div <- function(x)2 * x
    inv <- function(x)0.5 * x
    div0 <- d.div(0)
    d.inv <- function(x)0.5
  } else if (dvg == "KL" | (dvg == "Alpha" & identical(alpha, 1))) {
    div <- function(x)ifelse(x > 0, x * log(x), Inf)
    inv <- function(x)exp(x - 1)
    d.div <- function(x)ifelse(x > 0, 1 + log(x), -Inf)
    div0 <- -Inf
    d.inv <- function(x)inv(x)
  } else if (dvg == "revKL" | (dvg == "Alpha" & identical(alpha, 0))) {
    div <- function(x)ifelse(x > 0, -log(x), Inf)
    inv <- function(x)ifelse(x < 0, -1 / x, Inf)
    d.div <- function(x)ifelse(x > 0, -1 / x, -Inf)
    div0 <- -Inf
    d.inv <- function(x)ifelse(x < 0, 1 / x ^ 2, Inf)
  } else if (dvg == "Hellinger") {
    div <- function(x)ifelse(x >= 0, (sqrt(x) - 1) ^ 2, Inf)
    inv <- function(x)ifelse(x < 1, 1 / (1 - x) ^ 2, Inf)
    d.div <- function(x)ifelse(x > 0, 1 - 1 / sqrt(x), -Inf)
    div0 <- -Inf
    d.inv <- function(x)ifelse(x < 1, 2 / (1 - x) ^ 3, Inf)
  } else if (dvg == "Alpha" & alpha > 1) {
    div <- function(x)ifelse(x >= 0, (x ^ alpha - alpha * (x - 1) - 1) / (alpha * (alpha - 1)), Inf)
    div0 <- -1 / (alpha - 1)
    inv <- function(x)ifelse(x > div0, (1 + x * (alpha - 1)) ^ (1 / (alpha - 1)), 0)
    d.div <- function(x)ifelse(x > 0, (x ^ (alpha - 1) - 1) / (alpha -1), div0)
    d.inv <- function(x)ifelse(x > 0, (1 + x * (alpha - 1)) ^ (1 / (alpha - 1) - 1), )
    dvg <- paste(dvg, " a=", alpha)
  } else if (dvg == "Alpha" & alpha > 2) {
    div <- function(x)(x ^ alpha - alpha * (x - 1) - 1) / (alpha * (alpha - 1))
    inv <- function(x)(-2 * x * (1 + alpha)) ^ (-2 / (1 + alpha))
    d.div <- function(x)(-0.5 * (x ^ (-0.5 * (1 + alpha))) / (1 + alpha))
    div0 <- d.div(0)
    d.inv <- function(x)4 * (-2 * x * (1 + alpha)) ^ (- (3 + alpha) / (1 + alpha))
    dvg <- paste(dvg, " a=", alpha)
  } else if (dvg == "Triangular") {
    div <- function(x)ifelse(x >= 1, (x - 1) ^ 2 / (x + 1), Inf)
    inv <- function(x)-1 + 2 / sqrt(ifelse(x < 1, 1 - x, 0))
      # ifelse(x < 1, -1 + 2 / sqrt(1 - x), Inf)
    d.div <- function(x)(x - 1) * (x + 3) / (x + 1) ^ 2
    div0 <- -3
    # d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
    d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
  } else if (dvg == "LeCam") {
    div <- function(x)0.5 * (1 - x) / (x + 1)
  } else if (dvg == "user") {
    div0 <- div.usr$div0
    inv <- div.usr$inv
    if (max.l) div <- div.usr$div
    use.jac <- TRUE
  } else stop("The argument 'div' must be one of 'Chi2', 'KL', 'Hellinger', 'Alpha', 'Triangular', 'LeCam' or 'user' ") # check+++++++++++++

  z <- apply(X = x_data[, k, drop = FALSE], MARGIN = 1, FUN = f)
  min.z <- min(z)
  max.z <- max(z)
  if (min.d) {
    if (m < min.z | m > max.z) stop("m must be in the range of f(x)")
    }

  if (max.l) {
    max.div <- div(0) * (1 - 1 / n) + div(n) / n
    if (theta < 0 | theta > max.div) stop(paste("theta must be in the range [0,max.div], where max.div=", max.div,".")) # ++++++++++++++++
    }

  if (min.d & normalise == TRUE) {
    z <- (z - min.z) / (max.z - min.z)
    m <- (m - min.z) / (max.z - min.z)
  }

  if (min.d) constr <- function(L){
    q <- inv(pmax(div0, L[1] + L[2] * z)) / n
    C1 <- sum(q) - 1
    C2 <- sum(q * z) - m
    return(c(C1, C2))
  } else constr <- function(L){
    w <- inv(pmax(div0, L[1] + L[2] * z))
    C1 <- sum(w / n) - 1
    C2 <- sum(div(w) / n) - theta
    return(c(C1, C2))
  }

  if (is.null(start) & dvg != "user") start <- c(L1 = d.div(1), L2 = 0) ######

  if (exists("d.inv") & use.jac == TRUE) {
    if (min.d) {
    J <- function(L) {
      ind <- (L[1] + L[2] * z > div0)
      q1 <- d.inv(L[1] + L[2] * z) * ind / n
      J11 <- sum(q1)
      J12 <- sum(q1 * z)
      J22 <- sum(q1 * z * z)
      jac <- matrix(c(J11, J12, J12, J22), nrow = 2)
      return(jac)
      }
    } else {
      J <- function(L) {
        ind <- (L[1] + L[2] * fx > div0)
        q1 <- d.inv(L[1] + L[2] * z) * ind / n
        J11 <- sum(q1)
        J12 <- sum(q1 * z)
        J22 <- sum(q1 * z * z)
        J21 <- L[1] * J11 + L[2] * J12
        J22 <- L[1] * J12 + L[2] * J22
        jac <- matrix(c(J11, J21, J12, J22), nrow = 2)
        return(jac) ################## double check
      }
    }
    sol <- nleqslv::nleqslv(x = start, fn = constr, jac = J, ...)
  } else sol <- nleqslv::nleqslv(x = start, fn = constr, ...)

  if (sol$termcd != 1) warning(paste("nleqslv terminated with code ", sol$termcd))

  w <- inv(pmax(div0, sol$x[1] + sol$x[2] * z))########
  # if(sumRN) w <- w / mean(w) ################

  if (min.d) {
    m.ac <- sum(w * z / n) ############
    if(normalise == TRUE) {
      m <- min.z + (max.z - min.z) * m
      m.ac <- min.z + (max.z - min.z) * m.ac
    }
    err <- m - m.ac
    rel.err <- (err / m) * (m != 0)
    outcome <- c(required_moment = m, achieved_moment = m.ac, abs_error = err, rel_error = rel.err)
    constr_moment <- list("k" = k, "m" = m, "f" = f, "div" = dvg)
  } else {
    theta.ac <- sum(div(w) / n) ##############
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
