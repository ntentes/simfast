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

## does fit stay the same when data is shuffled?

df <- data.frame(preds, y, wts = weights)

set.seed(1337)
sf_formula <- simfast(formula = y ~ X1 + X2, weights = wts, data = df)

df_shuffled_index <- sample(1:NROW(df), NROW(df))
set.seed(1337)
sf_formula_shuffled <- simfast(formula = y ~ X1 + X2, weights = wts, data = df[df_shuffled_index, ])


testthat::test_that("simfast provides same fitted values when data is shuffled", {
  testthat::expect_equal(fitted(sf_formula)[df_shuffled_index], fitted(sf_formula_shuffled))
})

################################################################

## goodness of fit test requires data that was fit -- fit object should not inlcude missing data
set.seed(1337)
y_pois <- round(stats::rpois(500, lambda = mu))
y_pois_offset <- y_pois + 3
df_offset <- cbind(data.frame(preds), y_pois, wts = weights, os_vec = y_pois_offset)

df_offset_missing <- as.data.frame(lapply(df_offset, function(x) replace(x, sample(NROW(df_offset), 5), NA)))

sfm_os_missing <- simfast_m(x = as.matrix(df_offset_missing[, c("X1", "X2")]),
                            y = df_offset_missing[, "y_pois"],
                            weights = df_offset_missing[, "wts"],
                            offset = log(df_offset_missing[, "os_vec"]),
                            family = poisson)

sf_os_missing <- simfast(y_pois ~ X1 + X2 + offset(log(os_vec)),
                         weights = wts,
                         family = poisson,
                         data = df_offset_missing)

expect_equal(length(fitted(sfm_os_missing)), 475)
expect_equal(length(fitted(sf_os_missing)), 475)

