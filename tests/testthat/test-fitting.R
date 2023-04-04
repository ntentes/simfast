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

# need to implement coefficients function
# testthat::test_that("simfast_m provides same coefficients as simfast with formula", {
#   testthat::expect_equal(coefficients(sf_matrix), coefficients(sf_formula))
# })

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


#####################################################################


## test perfect separation data in binomial case

set.seed(1337)
ps_data <- data.frame(y = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                  x1 = c(0, 0, 0, 0, 1, 10, 25, 1, 1, 1, 10, 10, 10, 25, 25, 25),
                  x2 = c(0, 1, 10, 100, 0, 0, 0, 1, 10, 100, 1, 10, 100, 1, 10, 100))


## error when logLik "Delta/check" value is NaN
expect_error(simfast(y ~ x1 + x2, weights = rep(3, 16), family = binomial, data = ps_data))



#####################################################################


## expect warnings when values are missing either in data, offset or weights

set.seed(1337)

df_missing <- cbind(as.data.frame(lapply(df[, -3], function(x) replace(x, sample(NROW(df), 5), NA))), y = df[, 3])
df_offset_missing <- cbind(as.data.frame(lapply(df_offset[, -3], function(x) replace(x, sample(NROW(df_offset), 5), NA))), y_pois = df_offset[, 3])
df_response_missing <- df
df_response_missing$y[sample(NROW(df), 10)] <- NA


## fit stochastic missing predictor values
expect_warning(simfast_m(x = as.matrix(df_missing[, 1:2]), y = df_missing[, "y"]))
expect_warning(simfast(y ~ X1 + X2, data = df_missing))



## same fit whether missing predictor/weight values are left in or removed
set.seed(1337)
xmat_missing_index <- complete.cases(df_missing[, 1:3])
fitted0_missing <- fitted(simfast_m(x = as.matrix(df_missing[xmat_missing_index, 1:2]),
                                    y = df_missing[xmat_missing_index, "y"],
                                    weights = df_missing[xmat_missing_index, "wts"]))
set.seed(1337)
xmat_missing_index <- complete.cases(df_missing[, 1:2])
fitted1_missing <- fitted(simfast_m(x = as.matrix(df_missing[xmat_missing_index, 1:2]),
                                    y = df_missing[xmat_missing_index, "y"],
                                    weights = df_missing[xmat_missing_index, "wts"]))
set.seed(1337)
fitted2_missing <- fitted(simfast_m(x = as.matrix(df_missing[, 1:2]), y = df_missing[, "y"], weights = df_missing[, "wts"]))

expect_equal(fitted0_missing, fitted1_missing)
expect_equal(fitted1_missing, fitted2_missing)



## same fit whether missing offset values are left in or removed (formula vs matrix)
xmat_offset_missing_index <- complete.cases(df_offset_missing[, 1:4])
set.seed(1337)
fitted0_offset_missing <- fitted(simfast(formula = y_pois ~ X1 + X2 + offset(os_vec), weights = wts,
                                         data = df_offset_missing[xmat_offset_missing_index, ], family = poisson))
set.seed(1337)
fitted1_offset_missing <- fitted(simfast_m(x = as.matrix(df_offset_missing[, 1:2]),
                                           y = df_offset_missing$y_pois,
                                           weights = df_offset_missing$wts,
                                           offset = df_offset_missing$os_vec,
                                           family = poisson))

expect_equal(fitted0_offset_missing, fitted1_offset_missing)



## same fit whether missing response values are left in or removed (formula vs matrix)
response_missing_index <- complete.cases(df_response_missing[, "y"])

set.seed(1337)
fitted0_response_missing <- fitted(simfast(formula = y ~ X1 + X2, weights = wts,
                                         data = df_response_missing[response_missing_index, ]))
set.seed(1337)
fitted1_response_missing <- fitted(simfast_m(x = as.matrix(df_response_missing[, 1:2]),
                                             y = df_response_missing[, "y"],
                                             weights = df_response_missing[, "wts"]))

expect_equal(fitted0_response_missing, fitted1_response_missing)
