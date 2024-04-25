#' stress_mean_div
#'
#' Provides weights on simulated scenarios from a baseline stochastic model, such that (i) a stressed model component fulfills a moment constraint and scenario weights minimize a selected divergence to the baseline model, or (ii) a stressed model fulfills a divergence constraint and scenario weights maximize the moment of a model component. Case (i) is obtained by specifying the moment parameter '\code{m}', case (ii) by specifying the divergence constraint '\code{theta}'.
#'
#' @param x A vector, matrix or data frame containing realizations of random variables. Columns of \code{x} correspond to random variables; OR \cr
#' A \code{SWIM} object, where \code{x} corresponds to the underlying data of the \code{SWIM} object.
#' @param f A function to be applied to the columns of \code{x}.
#' @param k A vector indicating which columns of \code{x} the function \code{f} operates on. By default, the first column is selected.
#' @param m Numeric, the stressed moments of \code{f(x)}.
#' @param theta Numeric, the stressed divergence of the model.
#' @param dvg Character. One of "\code{Chi2}", "\code{KL}", "\code{revKL}", "\code{Hellinger}", "\code{Alpha}" or "\code{user}". See 'Details'
#' @param div.usr a list specifying the user divergence, see 'Details'.
#' @param alpha Numeric. The 'alpha' parameter when the Alpha divergence is used.
#' @param normalise Logical. If true, values of \code{f(x)} are linearly scaled to the unit interval. Not used when 'theta' is specified.
#' @param show Logical. If true, print the result of the call to \code{\link[nleqslv]{nleqslv}}.
#' @param names Character vector, the names of stressed models.
#' @param start A numeric vector, see 'Details'.
#' @param sumRN Logical. If true, the scenario weights are normalized so as to average to 1.
#' @param log Logical, the option to print weights' statistics.
#' @param use.jac Logical, should the Jacobian matrix be used in the call to \code{\link[nleqslv]{nleqslv}}. Set to TRUE when a user divergence is used and the argument div.usr$d.div is provided.
#' @param ... Additional arguments to be passed to \code{\link[nleqslv]{nleqslv}}.
#'
#' @seealso See \code{\link{stress_moment}} for a more flexible function allowing to perform multiple joint stresses under the KL divergence.
#'
#' @details When the argument \code{m} is specified, \code{stress_mean_div} solves the problem
#'
#' \deqn{\min_W E(h(W))}
#' under the constraints
#' \deqn{W\geq 0,\, E(W)=1,\,E(WZ)=m}
#' where \code{W=dQ/dP} is the Radon-Nikodym derivative and \code{Z} is a random variable obtained from the columns of \code{x} (see later). The moment constraint above can be restated as \code{E^{Q}(Z)=m}, where \code{E^{Q}} denotes the expectation under the stressed model. \code{stress_mean_div} solves the subsequent set of equations with respect to \code{\theta}, using \code{\link[nleqslv]{nleqslv}} from package \code{\link[nleqslv]{nleqslv}}:
#'
#' \deqn{E^Q( f(x) ) = E( f(x) * exp(theta * f(x)) ) = m.}
#'
#' When the argument \code{theta} is specified, \code{stress_mean_div} solves the problem
#'
#' \deqn{\max_W E(WZ)}
#' under the constraints
#' \deqn{W\geq 0,\, E(W)=1,\,E(h(W))\leq \theta}
#'
#' Here \code{E(WZ)=E^{Q}(Z)} where \code{E^{Q}} is the expectation operator under the stressed model  where \code{W=dQ/dP}
#'
#' There is no guarantee that the set of equations has a solution, or that the solution is unique. \code{stress_mean_div} will return a warning if the termination code provided by \code{nleqslv} is different from 1 (convergence has been achieved). It is recommended to check the result of the call to \code{nleqslv} using the \code{show} argument. The user is referred to the \code{\link[nleqslv]{nleqslv}} documentation for further details.
#'
#' The stressed moment \code{m} must be in the range of \code{f(x)}. The divergence value \code{theta} must be between 0 and \code{max.div=h(0) * (1-1/n)+h(n)/n}, where \code{h} is the divergence function and \code{n} is the number of simulations (number of rows of \code{x}).
#'
#' The function \code{f} operates on the columns of \code{x} indexed by the vector \code{k} to determine the values of the variable \code{Z}. By default \code{Z} is the first column of \code{x}.
#'
#' For a user specified divergence, the additional argument \code{div.usr} must be passed. This is a named list containing at least the following items: a function \code{inv} giving the inverse of the first derivative of the divergence function \code{h}; a numeric \code{div0} giving the the derivative at 0 of\code{h} (possibly \code{-Inf}); when \code{theta'} is specified (divergence constraint), \code{div.usr} must also contain the argumen \code{div} giving the divergence function \code{h}. Optionally, a function \code{d.div}, specifying the first derivative of the divergence, which will be used for the Jacobian of the +++++++. The divergence function must be one for which ###########

#'
#' For the "\code{Alpha}" divergence, the numeric parameter "\code{alpha}" must be provided (when \code{alpha=1} the KL divergence is used; when \code{alpha=2}, the \code{Chi2} divergence is used; when \code{alpha=0} the \code{revKL} divergence is used). The solution is numerically unstable for values of \code{alpha} +++++++++++
#'
#' The following divergences are considered:
#'
#' Chi2
#'
#' \deqn{h(x)=x^2-1}
#'
#' Kullback-Leibler
#'
#' \deqn{h(x)=x * log(x)}
#'
#' reverse Kullback-Leibler
#'
#' \deqn{h(x)=-log(x)}
#'
#' Hellinger
#'
#' \deqn{h(x)=(1 - x) ^ {-2}}
#'
#' Alpha
#'
#' \deqn{h(x)=\frac{x ^ \alpha - \alpha * (x - 1) - 1}{\alpha * (\alpha - 1)}}
#'
#'
#'
#'
#' with two elements, the starting values for the coefficients lambda1 and lambda2. Defaults to div$d.div(1) and 0 respectively, guaranteeing that the initial set of scenario weights is constant.
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

stress_mean_div <- function(x, f = function(x)x, k = 1, m = NULL, theta = NULL, dvg = c("Chi2", "KL", "revKL", "Hellinger", "Alpha", "user"), div.usr = NULL, alpha = NULL, normalise = TRUE, show = FALSE, names = NULL, start = NULL, sumRN = FALSE, log = FALSE, use.jac = FALSE, ...){

  min.d <- !is.null(m)
  max.l <- !is.null(theta)

  if (min.d + max.l != 1) stop("exactly one of m and theta must be provided")
  if (SWIM::is.SWIM(x)) x_data <- SWIM::get_data(x) else x_data <- as.matrix(x)
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
  }
  # else if (dvg == "Triangular") {
  #   div <- function(x)ifelse(x >= 1, (x - 1) ^ 2 / (x + 1), Inf)
  #   inv <- function(x)-1 + 2 / sqrt(ifelse(x < 1, 1 - x, 0))
  #     # ifelse(x < 1, -1 + 2 / sqrt(1 - x), Inf)
  #   d.div <- function(x)(x - 1) * (x + 3) / (x + 1) ^ 2
  #   div0 <- -3
  #   # d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
  #   d.inv <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
  # } else if (dvg == "LeCam") {
  #   div <- function(x)0.5 * (1 - x) / (x + 1)
  # } else if (dvg == "user") {
  #   div0 <- div.usr$div0
  #   inv <- div.usr$inv
  #   if (max.l) div <- div.usr$div
  #   use.jac <- TRUE
  # }
  else stop("The argument 'div' must be one of 'Chi2', 'KL', 'revKL', Hellinger', 'Alpha' or 'user' ")

  z <- apply(X = x_data[, k, drop = FALSE], MARGIN = 1, FUN = f)
  min.z <- min(z)
  max.z <- max(z)
  if (min.d) {
    if (m < min.z | m > max.z) stop("m must be in the range of f(x)")
    }

  if (max.l) {
    max.div <- div(0) * (1 - 1 / n) + div(n) / n
    if (theta < 0 | theta > max.div) stop(paste("theta must be in the range [0,max.div], where max.div=", max.div,"."))
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
  my_list <- SWIM::SWIM("x" = x_data, "new_weights" = new_weights, "type" = type, "specs" = constr)
  if (SWIM::is.SWIM(x)) my_list <- merge(x, my_list)

  if (show == TRUE) print(sol)

  if (log) {
    SWIM::summary_weights(my_list)
  }

  return(my_list)

}
