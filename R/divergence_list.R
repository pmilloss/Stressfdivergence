#' divergence
#'
#' @param dvg Character. One of "\code{Chi2}", "\code{KL}", "\code{revKL}", "\code{Hellinger}", "\code{Alpha}" or "\code{user}". See 'Details'
#'
#' @param alpha the value of the alpha parameter for the "Alpha" divergence
#'
#' @author Pietro Millossovich
#'
#' @details The following divergence functions are considered:
#'
#' Chi2 (\code{dvg=Chi2}):
#'
#' \deqn{h(x)=x^2-1}
#'
#' Kullback-Leibler (\code{dvg=KL}):
#'
#' \deqn{h(x)=x \dot \log(x)}
#'
#' reverse Kullback-Leibler (\code{dvg=revKL}):
#'
#' \deqn{h(x)=-\log(x)}
#'
#' Hellinger (\code{dvg=Hellinger}):
#'
#' \deqn{h(x)=(1 - x) ^ {-2}}.
#'
#' Alpha (\code{dvg=Alpha}):
#'
#' \deqn{h(x)=\frac{x ^ \alpha - \alpha * (x - 1) - 1}{\alpha * (\alpha - 1)}}
#'
#'
#' @examples
#' # example code
#'
#'
#'
#'

divergence <- function(dvg = c("Chi2", "KL", "revKL", "Hellinger", "Alpha"), alpha = NULL)
{

  if (dvg == "Chi2" | (dvg == "Alpha" & identical(alpha, 2))) {
    div <- function(x)ifelse(x >= 0, x ^ 2 - 1, Inf)
    d.div <- function(x)2 * x
    inv <- function(x)0.5 * x
    d.inv <- function(x)0.5
  } else if (dvg == "KL" | (dvg == "Alpha" & identical(alpha, 1))) {
    div <- function(x)ifelse(x > 0, x * log(x), Inf)
    inv <- function(x)exp(x - 1)
    d.div <- function(x)ifelse(x > 0, 1 + log(x), -Inf)
    d.inv <- function(x)inv(x)
  } else if (dvg == "revKL" | (dvg == "Alpha" & identical(alpha, 0))) {
    div <- function(x)ifelse(x > 0, -log(x), Inf)
    inv <- function(x)ifelse(x < 0, -1 / x, Inf)
    d.div <- function(x)ifelse(x > 0, -1 / x, -Inf)
    d.inv <- function(x)ifelse(x < 0, 1 / x ^ 2, Inf)
  } else if (dvg == "Hellinger") {
    div <- function(x)ifelse(x >= 0, (sqrt(x) - 1) ^ 2, Inf)
    inv <- function(x)ifelse(x < 1, 1 / (1 - x) ^ 2, Inf)
    d.div <- function(x)ifelse(x > 0, 1 - 1 / sqrt(x), -Inf)
    d.inv <- function(x)ifelse(x < 1, 2 / (1 - x) ^ 3, Inf)
  } else if (dvg == "Alpha"){
    h <- function(x)(x ^ alpha - alpha * (x - 1) - 1) / (alpha * (alpha - 1))
    h1 <- function(x)(x ^ (alpha - 1) - 1) / (alpha -1)
    g <- function(x)(1 + x * (alpha - 1)) ^ (1 / (alpha - 1))
    g1 <- function(x)(1 + x * (alpha - 1)) ^ ((2 - alpha) / (alpha - 1))
    d0 <- -1 / (alpha - 1)
    dvg <- paste(dvg, " a=", alpha)
    if (alpha > 1){
      div <- function(x)ifelse(x >= 0, h(x), Inf)
      inv <- function(x)ifelse(x > d0, g(x), 0)
      d.div <- function(x)ifelse(x >= 0, h1(x), d0)
      if (alpha > 2) d.inv <- function(x)ifelse(x > d0, g1(x), Inf) else d.inv <- function(x)ifelse(x > d0, g1(x), 0)
    } else if (alpha < 1) {
      div <- function(x)ifelse(x > 0, h(x), Inf)
      inv <- function(x)ifelse(x < d0, g(x), Inf)
      d.div <- function(x)ifelse(x > 0, h1(x), -Inf)
      d.inv <- function(x)ifelse(x < d0, g1(x), Inf)
    }
    dvg <- paste(dvg, " alpha=", alpha)
  }

  diver = list(div = div, inv = inv, d.div = d.div, d.inv = d.inv)
  return(diver)

  }
