# logistic interpretation

# forecast derived effect
x <- seq(0, 1, length.out = 100)
# s_k based on a sqrt shape

s_k <- 6 * (sqrt(x) - mean(sqrt(x)))

plot(x, s_k, type = "l")

# logit function
x_2 <- seq(-5, 5, length.out = 100)
plot(x_2, plogis(x_2), type = "l")
plot(x_2, plogis(x_2) * (1 - plogis(x_2)), type = "l")
plot(x_2, plogis(x_2) * (1 - plogis(x_2)) * (1 - 2 * plogis(x_2)), type = "l")


plogis_der <- function(x) plogis(x) * (1 - plogis(x))
plogis_der2 <- function(x) plogis(x) * (1 - plogis(x)) * (1 - 2 * plogis(x))
plogis_tailor2 <- function(h, x_0 = 0) {
  plogis(x_0) + plogis_der(x_0) * h + 1 / 2 * plogis_der2(x_0) * h^2
}


other_terms <- 0.1
intercept <- 0.2
x_0 <- other_terms + intercept
h <- seq(-1, 1, length.out = 100)
plot(x_0 + h, plogis(x_0 + h), type = "l")
plogis_der()
lines(
  x_0 + h,
  plogis_tailor2(h, x_0),
  col = "darkred"
)
plogis(intercept + other_terms)
s_k <- 6 * (sqrt(x) - mean(sqrt(x)))
plot(s_k, plogis(intercept + other_terms + s_k), type = "l")
lines(
  s_k,
  plogis_tailor2(
    h = s_k,
    x_0 = intercept + other_terms
  ),
  col = "darkred"
)


other_terms <- 0
intercept <- 0
s_k <- 6 * (sqrt(x) - mean(sqrt(x)))
plot(x, plogis(intercept + other_terms + s_k), type = "l")


plot(x, s_k, type = "l")
s_k_flat <- 6 * (sqrt(x) - mean(sqrt(x))) + 1 * (x - 0.5)
plot(x, s_k_flat, type = "l")

other_terms <- 0
intercept <- 0
s_k_flat <- 6 * (sqrt(x) - mean(sqrt(x))) + 1 * (x - 0.5)
plot(x, plogis(intercept + other_terms + s_k_flat), type = "l")


s_k <- 6 * (sqrt(x) - mean(sqrt(x)))
plot(x, s_k, type = "l")
lines(x, qlogis(x), col = "darkred")
