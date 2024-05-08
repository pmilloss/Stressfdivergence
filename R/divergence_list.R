#'
#'
#'
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
  } else if (dvg == "Alpha"){
    if (alpha > 1){
      div <- function(x)ifelse(x >= 0, (x ^ alpha - alpha * (x - 1) - 1) / (alpha * (alpha - 1)), Inf)
      div0 <- -1 / (alpha - 1)
      inv <- function(x)ifelse(x > div0, (1 + x * (alpha - 1)) ^ (1 / (alpha - 1)), 0)
      d.div <- function(x)ifelse(x > 0, (x ^ (alpha - 1) - 1) / (alpha -1), div0)
      d.inv <- function(x)ifelse(x > div0, (1 + x * (alpha - 1)) ^ (1 / (alpha - 1) - 1), )
      dvg <- paste(dvg, " a=", alpha)
    } else if (alpha > 1) {
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


}

  diver = list(div = div, inv = inv, div0 = div0, d.div = d.div, d.inv = d.inv)
  return(diver)

  }
