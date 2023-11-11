require(nleqslv)

stress_mean_div <- function(x, new_mean, div = c("Chi2", "KL", "Hellinger", "Alpha", "Triangular", "Jeffrey", "user"), g = NULL, f.der = NULL, g.der = NULL, p = rep(1 / length(x), length(x)), alpha = NULL, trace = FALSE,...){
  
  # f.der: first derivative of the divergence function
  # g: inverse function of f.der
  # g.der: first derivative of g (optional)
  
  if (div == "user" & (is.null(g) | is.null(f.der))) stop("For a user defined divergence, the arguments 'g' and 'f.der' must be provided")
  
  if (div == "Alpha" & is.null(alpha)) stop("For the Alpha divergence, the numeric argument 'alpha' must be provided")
  
  if (div == "Chi2") {
    g <- function(x)0.5 * x
    f.der <- function(x)2 * x
    f.der0 <- f.der(0)
    g.der <- function(x)0.5
  } else if (div == "KL") {
    g <- function(x)exp(x - 1)
    f.der <- function(x)1 + log(x)
    f.der0 <- -Inf
    g.der <- function(x)g(x)
  } else if (div == "Hellinger") {
    g <- function(x)1/(1 - x) ^ 2
    f.der <- function(x)(1 - 1 / sqrt(x))
    f.der0 <- -Inf
    g.der <- function(x)2 / (1 - x) ^ 3
  } else if (div == "Alpha") {
    g <- function(x)(-2 * x * (1 + alpha)) ^ (-2 / (1 + alpha))
    f.der <- function(x)(-0.5 * (x ^ (-0.5 * (1 + alpha))) / (1 + alpha))
    f.der0 <- f.der(0)
    g.der <- function(x)4 * (-2 * x * (1 + alpha)) ^ (- (3 + alpha) / (1 + alpha))
    #print(f.der0)
  } else if (div == "Triangular") {
    g <- function(x)-1 + 2 / sqrt(ifelse(x < 1, 1 - x, 0))
      # ifelse(x < 1, -1 + 2 / sqrt(1 - x), Inf)
    f.der <- function(x)(x - 1) * (x + 3) / (x + 1) ^ 2
    f.der0 <- -3
    # g.der <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
    g.der <- function(x)ifelse(x < 1, (1 - x) ^ -1.5, Inf)
  } else if (div == "Jeffrey") {
  } else if (div == "user") {
    f.der0 <- f.der(0)
  } else stop("The argument 'div' must be one of 'Chi2', 'KL', 'Hellinger', 'Alpha', 'Triangular', 'Jeffrey' or 'user' ")
  
  constraints <- function(L){
    RN <- g(pmax(f.der0, L[1] + L[2] * x))
    F1 <- sum(p * RN) - 1
    F2 <- sum(p * RN * x) - new_mean
    return(c(F1, F2))
  }
  
  starting <- c(L1 = f.der(1), L2 = 0)
  
  if (!is.null(g.der)){
    
    J <- function(L){
      RN.der <- g.der(L[1] + L[2] * x)
      ind <- (L[1] + L[2] * x > f.der0)
      J11 <- sum(p * RN.der * ind)
      J12 <- sum(p * RN.der * x * ind)
      J22 <- sum(p * RN.der * x * x * ind)
      mat <- matrix(c(J11, J12, J12, J22), nrow = 2, byrow = TRUE)
      return(mat)
    }
    
    res <- nleqslv(x = starting, fn = constraints, jac = J, ...)
  } else {
    res <- nleqslv(x = starting, fn = constraints, ...)
  }
  # print(c(res$x[1], res$x[2])) 
  # print(f.der0)
  #print(pmax(f.der0, res$x[1] + res$x[2] * x))
  RN <- g(pmax(f.der0, res$x[1] + res$x[2] * x))
  #print(RN)
  RN <- RN / mean(RN)
  if (trace)print(res)
  return(RN)
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




