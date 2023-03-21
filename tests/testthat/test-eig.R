# analytical expression
# equation 5 in https://projecteuclid.org/journals/statistical-science/volume-10/issue-3/Bayesian-Experimental-Design-A-Review/10.1214/ss/1177009939.full
ig <- function(X){
  # assume know sigma = 1, noise of linear model
  # assume R = 1, prior variance of parameter
  # assume k = 1, one parameter, no intercept
  sqrt(t(X) %*% X + 1)
}
ig_min <- optimize(ig, c(-10, 10))$minimum

# test model
priors <- c(
  brms::prior(normal(0, 1), class = "b"),
  brms::prior(constant(1), class = "sigma")
)
model <- brms::brm(
  test | weights(weight) ~ 0 + x,
  data = data.frame(x = 0, test = 0, weight = 0),
  prior = priors,
  sample_prior = "only",
  inits = 0
)

test_that("batch_nmc_eig works", {
  xs <- -10:10
  designs <- purrr::map(xs, ~ data.frame(x = .x))

  eigs <- batch_nmc_eig(designs, model)

  expect_length(eigs, length(xs)) # one eig for each design

  eig_min <- xs[which.min(eigs)]
  expect_equal(eig_min, ig_min) # same minimum value
})

test_that("inner_monte_carlo works", {
  N <- 100
  M <- 200
  mat <- matrix(0, nrow = M, ncol = N)
  lq_N <- inner_monte_carlo(mat, N = N, P = 1)

  expect_equal(length(lq_N), N)
  expect(all(lq_N == 0), "All lq_N != 0")

  mat <- matrix(1, nrow = M, ncol = N)
  lq_N <- inner_monte_carlo(mat, N = N, P = 1)
  expect(all(lq_N == 1), "All lq_N != log(1)")

  # this is kinda hacky but it is a nice analytical comparison
  # we can use this function to estimate entropy of a std normal
  N <- 1
  M <- 100000
  # we need to take the log of the log likelihood since inner_monte_carlo first exponentiates
  # but we cannot take log of log-likelihood because they are all negative values
  # in order to take log we must first multiply by the most negative number in set
  mat <- matrix(dnorm(rnorm(N*M), log = TRUE), nrow = M, ncol = N)
  min_mat <- min(mat)
  mat <- log(mat * min_mat)
  # we recover actual value by dividing min number back out at end.
  mc_entropy <- -exp(inner_monte_carlo(mat, N = N, P = 1)) / min_mat
  true_entropy <- 0.5 * log(2 * pi) + 0.5
  expect_equal(mc_entropy, true_entropy, tolerance = 0.01)
})

test_that("update_model works", {
  N <- 20
  M <- 20

  # stan kicks out lots of fitting diagnostic warnings so we suppress those.
  up_model <- suppressWarnings(update_model(model, N, M, cores = 1))
  fit <- as.matrix(up_model$fit)
  expect_equal(dim(fit)[1], N + M) # expect correct number of draws
  expect_equal(up_model$algorithm, "sampling")

  up_model <- suppressWarnings(update_model(model, N, M, cores = 2))
  fit <- as.matrix(up_model$fit)
  expect_equal(dim(fit)[1], N + M) # expect correct number of draws
  expect_equal(up_model$algorithm, "sampling")
})

test_that("predict_design works", {
  expect_error(
    predict_design(data.frame(y_=1), model, 1),
    regexp = "Column `y_` already exists and cannot overwrite."
  )

  N <- 7
  newdata <- data.frame(x = 1, dummy = 0)
  newpreds <- predict_design(newdata, model, N)
  expect_equal(nrow(newpreds), N)
  expect_setequal(colnames(newpreds), c("x", "test", "dummy"))
})
