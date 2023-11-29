library(SWIM)

set.seed(0)
n <- 10 ^ 5
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






res_KL <- stress_mean_div(x = x, m = t, div = "KL", show = T)
res_KL.1 <- stress_mean_div(x = x, m = t, div = "KL", normalise = FALSE, show = T)
SWIM::plot_weights(res_KL, xCol = 1)

all.equal(res_KL$new_weights, res_KL.1$new_weights)

res_KL.2 <- stress_mean_div(res_KL = x, m = t, div = "Alpha", show = T, alpha = 1)
all.equal(res_KL$new_weights, res_KL.2$new_weights)

res_KL <- stress_mean_div(x = x, m = t, div = "KL", show = T, control = list(maxit = 1000))
plot(x, res_KL)



res_Chi2 <- stress_mean_div(x = x, m = mean(x) * 1.2, dvg = "Chi2", show = T, use.jac = TRUE)
w_Chi2 <- res_Chi2$new_weights$`stress 1`
div_Chi2 <- mean((w_Chi2) ^ 2 - 1)

res_Chi2.1 <- stress_mean_div(x = x, theta = div_Chi2, dvg = "Chi2", show = T, use.jac = TRUE, control = list(allowSingular = TRUE))
all.equal(get_weights(res_Chi2), get_weights(res_Chi2.1))






plot_weights(res_Chi2, xCol = 1)


# Hellinger

res_Hell0 <- stress_mean_div(x = x, m = t, dvg = "Hellinger", normalise = FALSE, show = T)# fails after inspection
res_Hell <- stress_mean_div(x = x, m = t, div = "Hellinger", show = T)# fails after inspection
res_Hell1 <- stress_mean_div(x = x, m = t, div = "Hellinger", show = T, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 1000)) # works better, changing the scaling
plot(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 150))
plot(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 150), ylim = c(0, 10))
lines(x, res_Hell1$new_weights$`stress 1`, xlim = c(0, 150), ylim = c(0, 10))# madness


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

res_trian0 <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "Triangular", show = TRUE, normalise = FALSE)# fails
res_trian1 <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "Triangular", show = TRUE) # works after normalisation
res_trian2 <- stress_mean_div(x = x, m = mean(x) * 1.5, dvg = "Triangular", show = TRUE, normalise = FALSE, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 10000)) # works after "ad-hoc" normalisation

w_trian1 <- get_weights(res_trian1)
div_trian1 <- mean((w_trian1 - 1) ^ 2 / (w_trian1 + 1))






# reverse KL



