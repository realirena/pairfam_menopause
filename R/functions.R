#functions 

library(survival)
library(tidyverse)
library(haven)

# Generic comparison: population-marginal model cumulative hazard vs. a
# Turnbull/KM reference curve, with CI-band overlap diagnostics.
#
# surv_fn(t, draws): returns an (iterations x I) matrix of each individual's
#   survival probability at age t, for every posterior draw. This is the
#   only family-specific piece.
# draws: a list holding whatever surv_fn needs (varies by family, see below).
# age_grid: ages to evaluate at.
# turnbull_dat: your existing plot_dat (columns: time, cumhaz, cumhaz_lower, cumhaz_upper).

compare_model_to_turnbull <- function(surv_fn, draws, age_grid, turnbull_dat) {

  model_summary <- t(sapply(age_grid, function(t) {
    S_mat <- surv_fn(t, draws)
    S_pop <- rowMeans(S_mat)
    ch <- -log(S_pop)
    c(time = t,
      cumhaz = median(ch),
      cumhaz_lower = unname(quantile(ch, 0.025)),
      cumhaz_upper = unname(quantile(ch, 0.975)))
  }))
  model_dat <- as.data.frame(model_summary)

  turnbull_lower <- approx(turnbull_dat$time, turnbull_dat$cumhaz_lower, xout = age_grid)$y
  turnbull_upper <- approx(turnbull_dat$time, turnbull_dat$cumhaz_upper, xout = age_grid)$y

  overlap_check <- data.frame(
    age = age_grid,
    model_lower = model_dat$cumhaz_lower,
    model_upper = model_dat$cumhaz_upper,
    turnbull_lower = turnbull_lower,
    turnbull_upper = turnbull_upper
  )
  overlap_check$bands_overlap  <- with(overlap_check, model_lower <= turnbull_upper & turnbull_lower <= model_upper)
  overlap_check$model_too_high <- with(overlap_check, model_lower > turnbull_upper)
  overlap_check$model_too_low  <- with(overlap_check, turnbull_lower > model_upper)

  list(
    model_dat = model_dat,
    overlap_check = overlap_check,
    summary = list(
      n_grid_points   = sum(!is.na(overlap_check$bands_overlap)),
      age_range       = range(overlap_check$age[!is.na(overlap_check$bands_overlap)]),
      prop_non_overlap    = mean(!overlap_check$bands_overlap, na.rm = TRUE),
      prop_model_too_high = mean(overlap_check$model_too_high, na.rm = TRUE),
      prop_model_too_low  = mean(overlap_check$model_too_low, na.rm = TRUE)
    )
  )
}

weibull_surv_fn <- function(t, draws) {
  exp(-sweep(t / draws$hazard, 1, draws$c, "^"))
}
res_weibull <- compare_model_to_turnbull(
  weibull_surv_fn,
  draws = list(hazard = post_weibull$hazard, c = as.numeric(post_weibull$c)),
  age_grid, plot_dat
)

lognormal_surv_fn <- function(t, draws) {
  z <- sweep(log(t) - draws$mu, 1, draws$sigma, "/")
  pnorm(z, lower.tail = FALSE)
}
res_lognormal <- compare_model_to_turnbull(
  lognormal_surv_fn,
  draws = list(mu = post_lognormal$mu, sigma = as.numeric(post_lognormal$sigma_ln)),
  age_grid, plot_dat
)

loglogistic_surv_fn <- function(t, draws) {
  1 / (1 + sweep(t / draws$scale, 1, draws$shape, "^"))
}
res_loglogistic <- compare_model_to_turnbull(
  loglogistic_surv_fn,
  draws = list(scale = post_loglogistic$scale_ll, shape = as.numeric(post_loglogistic$shape_ll)),
  age_grid, plot_dat
)
