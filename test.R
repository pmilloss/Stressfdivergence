set.seed(0)
n <- 10 ^ 5
x <- rlnorm(n, meanlog = 1, sdlog = 1)
mean(x)
stress <- 1.2
t <- stress * mean(x)
range(x)

debug(stress_mean_div)
res_KL <- stress_mean_div(x = x, theta = 3, div = "KL", show = T)

res_KL <- stress_mean_div(x = x, m = t, div = "KL", show = T)
res_KL.1 <- stress_mean_div(x = x, m = t, div = "KL", normalise = FALSE, show = T)
SWIM::plot_weights(res_KL, xCol = 1)

all.equal(res_KL$new_weights, res_KL.1$new_weights)

res_KL.2 <- stress_mean_div(x = x, m = t, div = "Alpha", show = T, alpha = 1)
all.equal(res_KL$new_weights, res_KL.2$new_weights)

res_KL <- stress_mean_div(x = x, m = t, div = "KL", show = T, control = list(maxit = 1000))
plot(x, res_KL)



res_Chi2 <- stress_mean_div(x = x, m = mean(x) * 1.1, div = "Chi2", show = T)
SWIM::plot_weights(res_Chi2, xCol = 1)


res_Hell0 <- stress_mean_div(x = x, m = t, div = "Hellinger", normalise = FALSE, show = T)# fails after inspection
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

res_trian0 <- stress_mean_div(x = x, m = t, div = "Triangular", show = TRUE, normalise = FALSE)# fails
res_trian1 <- stress_mean_div(x = x, m = t, div = "Triangular", show = TRUE) # works after normalisation
res_trian2 <- stress_mean_div(x = x, m = t, div = "Triangular", show = TRUE, normalise = FALSE, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 1000)) # works after "ad-hoc" normalisation

res_trian <- stress_mean_div(x = x, new_mean = t, div = "Triangular", trace = T, xscalm = "auto", control = list(ftol = 10e-15, xtol = 10e-15, maxit = 1000))# this works better




