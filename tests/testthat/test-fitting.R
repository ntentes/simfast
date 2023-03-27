## Generate Data for Tests

set.seed(1337)

d <- 2
preds <- matrix(stats::runif(d*500, min = -4, max = 4), 500, d)

# Choose true alpha in R^2 with magnitude 1

alpha <- c(-2, 0.5)/sqrt(4.25)

# Set true ridge function (piecewise)

fn <- function(x){
  #  RANGE              #  VALUES
  (x< -3)             * (x+5)               +
  (x< -1)*(x>= -3)    * 2                   +
  (x< 1.525)*(x>= -1) * (0.1*(x + 1)^3 + 2) +
  (x>=1.525)          * (0.4*x + 3)
}

# Generate Gaussian response values with the same number of rows

indexvals <- preds %*% alpha
mu  <- fn(indexvals)
y <- stats::rnorm(500, mean = mu)

# Set weights

weights <- rep(c(5, 2, 3, 5, 2, 2, 1, 5, 3, 1), 50)


################################################################

set.seed(1337) # fitting method has a stochastic component
sf_matrix <- simfast_m(x = preds, y = y, weights = weights)

df <- data.frame(preds, y, wts = weights)
set.seed(1337)
sf_formula <- simfast(formula = y ~ X1 + X2, weights = wts, data = df)

testthat::test_that("simfast_m provides same coefficients as simfast with formula", {
  testthat::expect_equal(coefficients(sf_matrix), coefficients(sf_formula))
})

testthat::test_that("simfast_m provides same fitted values as simfast with formula -- stochastic fitting procedure", {
  testthat::expect_equal(fitted(sf_matrix), fitted(sf_formula))
})

# exact method fits
sf_matrix_exact <- simfast_m(x = preds, y = y, weights = weights, method = "exact")
sf_formula_exact <- simfast(formula = y ~ X1 + X2, weights = wts, data = df, method = "exact")

testthat::test_that("simfast_m provides same fitted values as simfast with formula -- exact fitting procedure", {
  testthat::expect_equal(fitted(sf_matrix_exact), fitted(sf_formula_exact))
})

#################################################################

## test offset in formula

set.seed(1337)
y_pois <- round(stats::rpois(500, lambda = mu))
y_pois_offset <- y_pois + 3

set.seed(1337) # fitting method has a stochastic component
sf_matrix_offset <- simfast_m(x = preds, y = y_pois, weights = weights, offset = y_pois_offset, family = poisson)

df_offset <- cbind(data.frame(preds), y_pois, wts = weights, os_vec = y_pois_offset)
set.seed(1337)
sf_formula_offset <- simfast(formula = y_pois ~ X1 + X2 + offset(os_vec), weights = wts,
                             data = df_offset, family = poisson)

testthat::test_that("simfast_m provides same fitted values as simfast with formula -- stochastic fitting with offset", {
  testthat::expect_equal(fitted(sf_matrix_offset), fitted(sf_formula_offset))
})
