library(SWIM)

set.seed(0)
x <- as.data.frame(cbind(
  "normal" <- rnorm(1000),
  "gamma" <- rgamma(1000, shape = 2)))

# stress entropy using SWIM
res1 <- stress(type = "VaR", x = x, alpha = 0.9, q_ratio = 1.05)

# stress entropy using stress_mean_div
newVaR <- quantile(x[, 1], probs = 0.9, type = 1) * 1.05
res2 <- stress_mean_div(res1, k = 1, f = function(x)1 * (x > newVaR), m = 0.1, dvg = "KL")
all.equal(get_weights(res2)[, 1], get_weights(res2)[, 2])

# stress mean with KL and chi2 divergence
newMean <- mean(x[, 1]) + 1
res1.1 <- stress_mean_div(x, k = 1, m = newMean, dvg = "KL")
res1.2 <- stress_mean_div(res1.1, k = 1, m = newMean, dvg = "Chi2")
res1.3 <- stress_mean_div(res1.2, k = 1, m = newMean, dvg = "Alpha", alpha = 1.5)
res1.4 <- stress_mean_div(res1.3, k = 1, m = newMean, dvg = "Alpha", alpha = 2.5)

summary(res1.4, xCol = 1)

get_specs(res1.4)
plot_weights(res1.4)

res1.5 <- stress_mean_div(x, k = 1, m = newMean, dvg = "revKL")
res1.51 <- stress_mean_div(x, k = 1, m = newMean * 1.1, dvg = "revKL")
res1.52 <- stress_mean_div(x, k = 1, m = newMean * 2, dvg = "revKL")




mean(x[, 1])
res1.6 <- stress_mean_div(x, k = 1, m = newMean, dvg = "Hellinger")
res1.61 <- stress_mean_div(x, k = 1, m = newMean * 1.1, dvg = "Hellinger")
res1.61 <- stress_mean_div(x, k = 1, m = newMean * 2, dvg = "Hellinger")



# stress the sum of the mean of columns 1, 2
res2.1 <- stress_mean_div(x, f = sum, k = 1:2, m = sum(colMeans(x)) * 1.1, dvg = "Chi2")
get_specs(res2.1)

# stress divergence
res3.1 <- stress_mean_div(x, theta = 10, dvg = "Chi2")
debug(stress_mean_div)
colMeans(x)
mean_stressed(res3.1, xCol = 1)
res3.2 <- stress_mean_div(res3.1, theta = 20, dvg = "Chi2")
mean_stressed(res3.2, xCol = 1)


quantile_stressed(res1, probs = 0.9, type = "(i-1)/(n-1)")
quantile_stressed(res1, probs = 0.9, type = "i/(n+1)")
quantile_stressed(res1, probs = 0.9, type = "i/n")
quantile_stressed(res2, probs = 0.9, wCol = 2, type = "i/n")

plot_weights(res2, wCol = 1)
plot_weights(res2, wCol = 2)

sum(get_weights(res2, wCol = 2))

## calling stress_VaR directly
## stressing "gamma"
res2 <- stress_VaR(x = x, alpha = 0.9,
  q_ratio = c(1.03, 1.05), k = 2)
get_specs(res2)
summary(res2)


set.seed(0)
n <- 10 ^ 4
x <- rlnorm(n, meanlog = 1, sdlog = 1)
mean(x)
stress <- 1.2
t <- stress * mean(x)
range(x)

# res_KL <- stress_mean_div(x = x, theta = 0.1, dvg = "KL", show = T, use.jac = TRUE) # not run - does not work as Jacobian is singular
res_KL <- stress_mean_div(x = x, theta = 0.1, dvg = "KL", show = T, use.jac = TRUE, control = list(allowSingular = TRUE))
res_KL1 <- stress_mean_div(x = x, theta = 0.1, dvg = "KL", show = T) # do not use Jacobian


all.equal(get_weights(res_KL), get_weights(res_KL1))

summ_KL <- summary(res_KL, base = TRUE)
mean_KL <- summ_KL$`stress 1`$X1[1]

# now stress the mean
res_KL2 <- stress_mean_div(x, m = mean_KL, dvg = "KL", show = T)
summary(res_KL2)
summary(res_KL)

all.equal(get_weights(res_KL), get_weights(res_KL2))






res_KL <- stress_mean_div(x = x, m = t, div = "KL", show = T)
res_KL.1 <- stress_mean_div(x = x, m = t, div = "KL", normalise = FALSE, show = T)
SWIM::plot_weights(res_KL, xCol = 1)

all.equal(res_KL$new_weights, res_KL.1$new_weights)

res_KL.2 <- stress_mean_div(res_KL = x, m = t, div = "Alpha", show = T, alpha = 1)
all.equal(res_KL$new_weights, res_KL.2$new_weights)

res_KL <- stress_mean_div(x = x, m = t, div = "KL", show = T, control = list(maxit = 1000))
plot(x, res_KL)


# reverse KL
# res_KL <- stress_mean_div(x = x, theta = 0.1, dvg = "KL", show = T, use.jac = TRUE) # not run - does not work as Jacobian is singular
res_KL <- stress_mean_div(x = x, theta = 0.1, dvg = "KL", show = T, use.jac = TRUE, control = list(allowSingular = TRUE))
res_KL1 <- stress_mean_div(x = x, theta = 0.1, dvg = "KL", show = T) # do not use Jacobian


all.equal(get_weights(res_KL), get_weights(res_KL1))

summ_KL <- summary(res_KL, base = TRUE)
mean_KL <- summ_KL$`stress 1`$X1[1]

# now stress the mean
res_KL2 <- stress_mean_div(x, m = mean_KL, dvg = "KL", show = T)
summary(res_KL2)
summary(res_KL)

all.equal(get_weights(res_KL), get_weights(res_KL2))


res_Chi2 <- stress_mean_div(x = x, m = mean(x) * 1.2, dvg = "Chi2", show = T, use.jac = TRUE)
w_Chi2 <- res_Chi2$new_weights$`stress 1`
div_Chi2 <- mean((w_Chi2) ^ 2 - 1)

res_Chi2.1 <- stress_mean_div(x = x, theta = div_Chi2, dvg = "Chi2", show = T, use.jac = TRUE, control = list(allowSingular = TRUE))
all.equal(get_weights(res_Chi2), get_weights(res_Chi2.1))






plot_weights(res_Chi2, xCol = 1)


# Hellinger

t <- mean(x[, 2]) * 1.2
res_Hell0 <- stress_mean_div(k = 2, x = x, m = t, dvg = "Hellinger", normalise = FALSE, show = T)
plot_weights(res_Hell0, xCol = 2)


, y_limits = c(0, 20))
res_Hell <- stress_mean_div(x = res_Hell0, m = t, dvg = "Hellinger", show = T)# normalisation
res_Hell1 <- stress_mean_div(x = x, m = t, dvg = "Hellinger", show = T)# normalisation

res_Hell2 <- stress_mean_div(x = x, m = t, dvg = "Hellinger", show = T, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 1000)) # works better, changing the scaling
plot_weights(res_Hell2, y_limits = c(0, 20))

plot(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 150))
plot(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 150), ylim = c(0, 10))
lines(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 150), ylim = c(0, 10))# madness
plot(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 10), ylim = c(0, 2))

all(diff(res_Hell1$new_weights$`stress 1`[order(x)]) >= 0)


res_Hell4 <- stress_mean_div(x = x, theta = 0.5, dvg = "Hellinger", show = T)






alpha <- seq(-1.5, 3, by = 0.5)
res_alpha <- stress_mean_div(x = x, m = t, div = "Alpha", alpha = -2, show = FALSE)
for (i in seq_along(alpha))res_alpha <- stress_mean_div(x = res_alpha, m = t, div = "Alpha", alpha = alpha[i], show = FALSE)


plot(x, res_alpha)
res_alpha2 <- stress_mean_div(x = x, m = t, div = "Alpha", alpha = -1.5, show = TRUE)
res_alpha2 <- stress_mean_div(x = x, new_mean = t, div = "Alpha", alpha = 0.5, trace = T) #ok
res_alpha3 <- stress_mean_div(x = x, new_mean = t, div = "Alpha", alpha = 1.5, trace = T) #fails
plot(x, res_alpha)
plot(x, res_alpha2)

res_alpha2 <- stress_mean_div(x = x, m = mean(x) * 1.1, div = "Alpha", alpha = 2.5, show = TRUE, control = list(ftol = 10e-5, xtol = 10e-5, maxit = 10000), method = "Broyden", global = "cline") # fails


res_alpha3 <- stress_mean_div(x = x, m = mean(x) * 1.1, div = "Alpha", alpha = 4, show = TRUE) # fails




# triangular

res_trian0 <- stress_mean_div(k = 2, x = x, m = t, dvg = "Triangular", show = TRUE, normalise = FALSE)# fails


res_trian0 <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "Triangular", show = TRUE, normalise = FALSE)# fails
res_trian1 <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "Triangular", show = TRUE) # works after normalisation
res_trian2 <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "Triangular", show = TRUE, normalise = FALSE, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 10000)) # works after "ad-hoc" normalisation

w_trian1 <- get_weights(res_trian1)
div_trian1 <- mean((w_trian1 - 1) ^ 2 / (w_trian1 + 1))






# reverse KL

res_revKL <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "user", div.usr = list(div = function(x)-log(x), inv = function(x)-1 / x, div0 = -Inf, d.div = function(x)-1 / x), show = TRUE, normalise = FALSE)# fails

debug(stress_mean_div)
