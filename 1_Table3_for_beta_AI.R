library(dplyr)
library(tidyr)
library(lubridate)
library(writexl)

## --------------------------------------------------------------------------
## 0. User settings
## --------------------------------------------------------------------------

sample_start <- as.Date("1965-01-01")
sample_end   <- as.Date("1989-12-31")

# The original paper does not report its minimum daily-observation rule.
# This is an explicit replication choice; later check 100/150/200 sensitivity.
min_daily_obs <- 100L

# Map these three names to columns in crsp_data.
# Typical alternatives are date/dlycaldt and ret/dlyret.
daily_date_col       <- "date"
daily_stock_ret_col  <- "ret"

# This must be a DAILY market-return series, preferably the CRSP
# value-weighted index return, repeated consistently by date in crsp_data.
# If crsp_data does not contain it, create a separate market_daily object
# as described immediately below and set use_market_from_crsp_data <- FALSE.
daily_market_ret_col <- "vwretd"
use_market_from_crsp_data <- TRUE

# When use_market_from_crsp_data is FALSE, create market_daily before running:
# market_daily <- YOUR_DAILY_MARKET_DATA %>%
#   transmute(
#     date = as.Date(YOUR_DATE_COLUMN),
#     market_ret = as.numeric(YOUR_CRSP_VW_RETURN_COLUMN)
#   )

## --------------------------------------------------------------------------
## 1. Input checks
## --------------------------------------------------------------------------

required_monthly <- c(
  "date", "permno", "mthret", "mthcap", "cum_ret_6", "st6"
)

missing_monthly <- setdiff(required_monthly, names(dat))

if (length(missing_monthly) > 0L) {
  stop(
    "dat is missing required columns: ",
    paste(missing_monthly, collapse = ", ")
  )
}

monthly_duplicate_keys <- dat %>%
  count(date, permno, name = "n") %>%
  filter(n != 1L)

if (nrow(monthly_duplicate_keys) > 0L) {
  stop("dat has duplicate date-permno keys. Fix them before forming portfolios.")
}

required_daily <- c("permno", daily_date_col, daily_stock_ret_col)

if (use_market_from_crsp_data) {
  required_daily <- c(required_daily, daily_market_ret_col)
}

missing_daily <- setdiff(required_daily, names(crsp_data))

if (length(missing_daily) > 0L) {
  stop(
    "crsp_data is missing required columns: ",
    paste(missing_daily, collapse = ", "),
    ". Change the daily column mappings in Section 0."
  )
}

## --------------------------------------------------------------------------
## 2. Standardize daily stock and market data
## --------------------------------------------------------------------------

daily_stock <- crsp_data %>%
  transmute(
    date = as.Date(.data[[daily_date_col]]),
    permno = permno,
    stock_ret = as.numeric(.data[[daily_stock_ret_col]])
  )

daily_stock_duplicate_keys <- daily_stock %>%
  count(date, permno, name = "n") %>%
  filter(n != 1L)

if (nrow(daily_stock_duplicate_keys) > 0L) {
  stop(
    "Daily stock data have duplicate date-permno keys. ",
    "Do not continue until the duplicate rows are resolved."
  )
}

if (use_market_from_crsp_data) {
  market_daily_check <- crsp_data %>%
    transmute(
      date = as.Date(.data[[daily_date_col]]),
      market_ret = as.numeric(.data[[daily_market_ret_col]])
    ) %>%
    filter(is.finite(market_ret)) %>%
    distinct() %>%
    count(date, name = "n") %>%
    filter(n != 1L)

  if (nrow(market_daily_check) > 0L) {
    stop(
      "crsp_data contains multiple daily market returns for the same date. ",
      "Construct a unique market_daily series first."
    )
  }

  market_daily <- crsp_data %>%
    transmute(
      date = as.Date(.data[[daily_date_col]]),
      market_ret = as.numeric(.data[[daily_market_ret_col]])
    ) %>%
    filter(is.finite(market_ret)) %>%
    distinct(date, market_ret)
} else {
  if (!exists("market_daily")) {
    stop(
      "Create market_daily with unique date and market_ret columns before running."
    )
  }

  market_daily <- market_daily %>%
    transmute(
      date = as.Date(date),
      market_ret = as.numeric(market_ret)
    )
}

market_daily_duplicate_keys <- market_daily %>%
  count(date, name = "n") %>%
  filter(n != 1L)

if (nrow(market_daily_duplicate_keys) > 0L) {
  stop("market_daily must have exactly one row per date.")
}

## --------------------------------------------------------------------------
## 3. Estimate prior-calendar-year Scholes-Williams stock betas
## --------------------------------------------------------------------------

safe_slope <- function(y, x, min_obs) {
  valid <- is.finite(y) & is.finite(x)

  if (sum(valid) < min_obs) {
    return(NA_real_)
  }

  x_valid <- x[valid]
  y_valid <- y[valid]
  x_var <- stats::var(x_valid)

  if (!is.finite(x_var) || x_var <= 0) {
    return(NA_real_)
  }

  stats::cov(y_valid, x_valid) / x_var
}

market_daily <- market_daily %>%
  filter(is.finite(market_ret)) %>%
  mutate(estimation_year = year(date)) %>%
  group_by(estimation_year) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    market_lag = lag(market_ret, 1L),
    market_lead = lead(market_ret, 1L),
    market_rho = stats::cor(
      market_ret,
      market_lag,
      use = "complete.obs"
    )
  ) %>%
  ungroup()

daily_beta_input <- daily_stock %>%
  mutate(estimation_year = year(date)) %>%
  left_join(
    market_daily,
    by = c("date", "estimation_year"),
    relationship = "many-to-one"
  )

stock_beta_year <- daily_beta_input %>%
  group_by(permno, estimation_year) %>%
  summarise(
    n_daily = sum(is.finite(stock_ret) & is.finite(market_ret)),
    beta_lag = safe_slope(stock_ret, market_lag, min_daily_obs),
    beta_0 = safe_slope(stock_ret, market_ret, min_daily_obs),
    beta_lead = safe_slope(stock_ret, market_lead, min_daily_obs),
    market_rho = first(market_rho),
    .groups = "drop"
  ) %>%
  mutate(
    sw_denominator = 1 + 2 * market_rho,
    beta_sw = if_else(
      is.finite(beta_lag) &
        is.finite(beta_0) &
        is.finite(beta_lead) &
        is.finite(sw_denominator) &
        abs(sw_denominator) > 1e-8,
      (beta_lag + beta_0 + beta_lead) / sw_denominator,
      NA_real_
    ),
    formation_year = estimation_year + 1L
  ) %>%
  select(
    permno,
    formation_year,
    beta_sw,
    n_daily,
    beta_lag,
    beta_0,
    beta_lead,
    market_rho
  )

rm(
  list = intersect(
    c(
      "daily_stock",
      "daily_beta_input",
      "market_daily_check",
      "daily_stock_duplicate_keys",
      "market_daily_duplicate_keys",
      "missing_daily",
      "missing_monthly",
      "monthly_duplicate_keys",
      "required_daily",
      "required_monthly",
      "safe_slope"
    ),
    ls()
  )
)

## --------------------------------------------------------------------------
## 4. Form All, size-tercile, and beta-tercile momentum portfolios
## --------------------------------------------------------------------------

# All-stock momentum deciles already constructed by the user's st6 variable.
formation_all <- dat %>%
  filter(!is.na(st6), is.finite(cum_ret_6)) %>%
  transmute(
    date,
    permno,
    subsample = "All",
    momr = as.integer(st6)
  )

# Size is measured at the end of the formation month because the signal uses
# returns through that month and holding starts in the following month.
# If the user's mthcap definition differs, replace mthcap here explicitly.
formation_size <- dat %>%
  filter(
    is.finite(cum_ret_6),
    is.finite(mthcap),
    mthcap > 0
  ) %>%
  group_by(date) %>%
  mutate(size_tile = ntile(mthcap, 3L)) %>%
  group_by(date, size_tile) %>%
  mutate(momr = ntile(cum_ret_6, 10L)) %>%
  ungroup() %>%
  transmute(
    date,
    permno,
    subsample = paste0("S", size_tile),
    momr = as.integer(momr)
  )

# Each formation-year beta is estimated exclusively from the preceding
# calendar year's daily returns.
formation_beta <- dat %>%
  filter(is.finite(cum_ret_6)) %>%
  transmute(
    date,
    permno,
    cum_ret_6,
    formation_year = year(date)
  ) %>%
  left_join(
    stock_beta_year %>%
      select(permno, formation_year, beta_sw, n_daily),
    by = c("permno", "formation_year"),
    relationship = "many-to-one"
  ) %>%
  filter(is.finite(beta_sw)) %>%
  group_by(date) %>%
  mutate(beta_tile = ntile(beta_sw, 3L)) %>%
  group_by(date, beta_tile) %>%
  mutate(momr = ntile(cum_ret_6, 10L)) %>%
  ungroup() %>%
  transmute(
    date,
    permno,
    subsample = paste0("B", beta_tile),
    momr = as.integer(momr)
  )

table3_formation <- bind_rows(
  formation_all,
  formation_size,
  formation_beta
)

formation_duplicate_keys <- table3_formation %>%
  count(date, permno, subsample, name = "n") %>%
  filter(n != 1L)

if (nrow(formation_duplicate_keys) > 0L) {
  stop("Table III formation data contain duplicate stock assignments.")
}

formation_sort_check <- table3_formation %>%
  group_by(date, subsample) %>%
  summarise(
    n_stocks = n(),
    n_portfolios = n_distinct(momr),
    min_portfolio = min(momr),
    max_portfolio = max(momr),
    .groups = "drop"
  )

bad_formation_sorts <- formation_sort_check %>%
  filter(
    n_portfolios != 10L |
      min_portfolio != 1L |
      max_portfolio != 10L
  )

if (nrow(bad_formation_sorts) > 0L) {
  warning(
    "Some formation-month subsamples do not contain all 10 momentum portfolios. ",
    "Inspect bad_formation_sorts."
  )
}

rm(formation_all, formation_size, formation_beta, formation_duplicate_keys)

## --------------------------------------------------------------------------
## 5. Expand each formation cohort to six holding months
## --------------------------------------------------------------------------

holding_return <- dat %>%
  transmute(
    holding_date = as.Date(date),
    permno,
    mthret
  )

holding_return_duplicate_keys <- holding_return %>%
  count(holding_date, permno, name = "n") %>%
  filter(n != 1L)

if (nrow(holding_return_duplicate_keys) > 0L) {
  stop("Holding-return data contain duplicate holding_date-permno keys.")
}

table3_cohort_dat <- table3_formation %>%
  mutate(K = 6L) %>%
  uncount(K, .remove = TRUE, .id = "holding_month") %>%
  mutate(
    holding_date = date %m+% months(holding_month)
  ) %>%
  left_join(
    holding_return,
    by = c("holding_date", "permno"),
    relationship = "many-to-one"
  )

table3_cohort <- table3_cohort_dat %>%
  group_by(subsample, holding_date, date, momr) %>%
  summarise(
    cohort_ret = if (
      all(is.na(mthret))
    ) {
      NA_real_
    } else {
      mean(mthret, na.rm = TRUE)
    },
    n_stocks = sum(is.finite(mthret)),
    .groups = "drop"
  )

table3_portfolio_ret <- table3_cohort %>%
  group_by(subsample, holding_date, momr) %>%
  summarise(
    portfolio_ret = if (
      all(is.na(cohort_ret))
    ) {
      NA_real_
    } else {
      mean(cohort_ret, na.rm = TRUE)
    },
    n_cohorts = sum(is.finite(cohort_ret)),
    .groups = "drop"
  ) %>%
  filter(
    holding_date >= sample_start,
    holding_date <= sample_end
  ) %>%
  mutate(
    subsample = factor(
      subsample,
      levels = c("All", "S1", "S2", "S3", "B1", "B2", "B3")
    )
  ) %>%
  arrange(subsample, holding_date, momr)

rm(
  holding_return,
  holding_return_duplicate_keys,
  table3_cohort_dat,
  table3_formation
)

## --------------------------------------------------------------------------
## 6. Validate the completed Table III return panel
## --------------------------------------------------------------------------

table3_panel_check <- table3_portfolio_ret %>%
  group_by(subsample) %>%
  summarise(
    rows = n(),
    months = n_distinct(holding_date),
    portfolios = n_distinct(momr),
    min_date = min(holding_date),
    max_date = max(holding_date),
    missing = sum(is.na(portfolio_ret)),
    nonfinite = sum(!is.finite(portfolio_ret)),
    .groups = "drop"
  )

table3_duplicate_keys <- table3_portfolio_ret %>%
  count(subsample, holding_date, momr, name = "n") %>%
  filter(n != 1L)

table3_bad_cohorts <- table3_portfolio_ret %>%
  filter(n_cohorts != 6L)

if (nrow(table3_duplicate_keys) > 0L) {
  stop("Completed Table III panel has duplicate subsample-date-portfolio keys.")
}

if (any(table3_panel_check$nonfinite > 0L)) {
  warning("Completed Table III panel contains nonfinite returns.")
}

if (nrow(table3_bad_cohorts) > 0L) {
  warning(
    "Some Table III portfolio-months have fewer than six active cohorts. ",
    "Inspect table3_bad_cohorts before interpreting the table."
  )
}

## --------------------------------------------------------------------------
## 7. Panel A: average monthly returns and t-statistics
## --------------------------------------------------------------------------

mean_t_result <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(
      tibble(
        estimate = NA_real_,
        t_stat = NA_real_,
        n_months = length(x)
      )
    )
  }

  test <- stats::t.test(x, mu = 0)

  tibble(
    estimate = unname(test$estimate),
    t_stat = unname(test$statistic),
    n_months = length(x)
  )
}

panel_a_portfolios <- table3_portfolio_ret %>%
  group_by(subsample, momr) %>%
  group_modify(~ mean_t_result(.x$portfolio_ret)) %>%
  ungroup() %>%
  mutate(portfolio = paste0("P", momr)) %>%
  select(subsample, portfolio, estimate, t_stat, n_months)

table3_spread <- table3_portfolio_ret %>%
  select(subsample, holding_date, momr, portfolio_ret) %>%
  pivot_wider(
    names_from = momr,
    values_from = portfolio_ret,
    names_prefix = "P"
  ) %>%
  mutate(portfolio_ret = P10 - P1)

panel_a_spread <- table3_spread %>%
  group_by(subsample) %>%
  group_modify(~ mean_t_result(.x$portfolio_ret)) %>%
  ungroup() %>%
  mutate(portfolio = "P10-P1") %>%
  select(subsample, portfolio, estimate, t_stat, n_months)

panel_a <- bind_rows(panel_a_portfolios, panel_a_spread) %>%
  mutate(
    portfolio = factor(
      portfolio,
      levels = c(paste0("P", 1:10), "P10-P1")
    )
  ) %>%
  arrange(subsample, portfolio)

# Joint test for Panel A:
# H0: E[P1] = E[P2] = ... = E[P10].
panel_a_joint_test <- function(data) {
  wide <- data %>%
    select(holding_date, momr, portfolio_ret) %>%
    pivot_wider(
      names_from = momr,
      values_from = portfolio_ret,
      names_prefix = "P"
    ) %>%
    arrange(holding_date)

  returns <- as.matrix(wide[, paste0("P", 1:10)])
  returns <- returns[complete.cases(returns), , drop = FALSE]

  time_count <- nrow(returns)
  restriction_count <- 9L
  restriction <- matrix(0, nrow = restriction_count, ncol = 10L)

  for (j in 1:restriction_count) {
    restriction[j, j] <- 1
    restriction[j, 10L] <- -1
  }

  restricted_returns <- returns %*% t(restriction)
  restricted_mean <- colMeans(restricted_returns)
  restricted_cov <- stats::cov(restricted_returns)

  t_squared <- tryCatch(
    time_count * drop(
      t(restricted_mean) %*%
        solve(restricted_cov) %*%
        restricted_mean
    ),
    error = function(e) NA_real_
  )

  f_stat <- (
    (time_count - restriction_count) /
      (restriction_count * (time_count - 1L))
  ) * t_squared

  tibble(
    f_stat = f_stat,
    df1 = restriction_count,
    df2 = time_count - restriction_count,
    p_value = stats::pf(
      f_stat,
      df1 = restriction_count,
      df2 = time_count - restriction_count,
      lower.tail = FALSE
    ),
    n_months = time_count
  )
}

panel_a_f_test <- table3_portfolio_ret %>%
  group_by(subsample) %>%
  group_modify(~ panel_a_joint_test(.x)) %>%
  ungroup()

## --------------------------------------------------------------------------
## 8. Panel B: market-model alphas and t-statistics
## --------------------------------------------------------------------------

fama_duplicate_keys <- fama %>%
  count(holding_date, name = "n") %>%
  filter(n != 1L)

if (nrow(fama_duplicate_keys) > 0L) {
  stop("fama has duplicate holding_date keys.")
}

table3_regression_data <- table3_portfolio_ret %>%
  select(subsample, holding_date, momr, portfolio_ret) %>%
  left_join(
    fama %>% select(holding_date, Mkt.RF, RF),
    by = "holding_date",
    relationship = "many-to-one"
  ) %>%
  mutate(portfolio_excess = portfolio_ret - RF)

market_model_result <- function(data) {
  fit <- stats::lm(
    portfolio_excess ~ Mkt.RF,
    data = data
  )

  coefficient_table <- summary(fit)$coefficients

  tibble(
    alpha = unname(coef(fit)[["(Intercept)"]]),
    alpha_t = unname(coefficient_table["(Intercept)", "t value"]),
    market_beta = unname(coef(fit)[["Mkt.RF"]]),
    n_months = stats::nobs(fit)
  )
}

panel_b_portfolios <- table3_regression_data %>%
  group_by(subsample, momr) %>%
  group_modify(~ market_model_result(.x)) %>%
  ungroup() %>%
  mutate(portfolio = paste0("P", momr)) %>%
  select(
    subsample,
    portfolio,
    alpha,
    alpha_t,
    market_beta,
    n_months
  )

panel_b_spread_data <- table3_spread %>%
  transmute(
    subsample,
    holding_date,
    portfolio_ret
  ) %>%
  left_join(
    fama %>% select(holding_date, Mkt.RF, RF),
    by = "holding_date",
    relationship = "many-to-one"
  ) %>%
  mutate(
    # RF cancels in P10 minus P1, but using this common interface is harmless.
    portfolio_excess = portfolio_ret
  )

panel_b_spread <- panel_b_spread_data %>%
  group_by(subsample) %>%
  group_modify(~ market_model_result(.x)) %>%
  ungroup() %>%
  mutate(portfolio = "P10-P1") %>%
  select(
    subsample,
    portfolio,
    alpha,
    alpha_t,
    market_beta,
    n_months
  )

panel_b <- bind_rows(panel_b_portfolios, panel_b_spread) %>%
  mutate(
    portfolio = factor(
      portfolio,
      levels = c(paste0("P", 1:10), "P10-P1")
    )
  ) %>%
  arrange(subsample, portfolio)

# Gibbons-Ross-Shanken joint test for Panel B:
# H0: alpha_P1 = ... = alpha_P10 = 0.
grs_test_one_subsample <- function(data) {
  wide <- data %>%
    select(holding_date, momr, portfolio_ret, Mkt.RF, RF) %>%
    pivot_wider(
      names_from = momr,
      values_from = portfolio_ret,
      names_prefix = "P"
    ) %>%
    arrange(holding_date)

  required_columns <- c(paste0("P", 1:10), "Mkt.RF", "RF")
  wide <- wide[complete.cases(wide[, required_columns]), ]

  raw_returns <- as.matrix(wide[, paste0("P", 1:10)])
  excess_returns <- sweep(raw_returns, 1L, wide$RF, FUN = "-")
  factors <- matrix(wide$Mkt.RF, ncol = 1L)

  time_count <- nrow(excess_returns)
  asset_count <- ncol(excess_returns)
  factor_count <- ncol(factors)

  design <- cbind(1, factors)
  coefficients <- solve(
    crossprod(design),
    crossprod(design, excess_returns)
  )

  alpha <- coefficients[1L, ]
  residuals <- excess_returns - design %*% coefficients
  residual_cov <- crossprod(residuals) /
    (time_count - factor_count - 1L)

  factor_mean <- colMeans(factors)
  factor_cov <- stats::cov(factors)

  numerator <- tryCatch(
    drop(t(alpha) %*% solve(residual_cov) %*% alpha),
    error = function(e) NA_real_
  )

  denominator_adjustment <- tryCatch(
    1 + drop(
      t(factor_mean) %*%
        solve(as.matrix(factor_cov)) %*%
        factor_mean
    ),
    error = function(e) NA_real_
  )

  df1 <- asset_count
  df2 <- time_count - asset_count - factor_count

  grs_stat <- (
    (time_count - asset_count - factor_count) / asset_count
  ) * numerator / denominator_adjustment

  tibble(
    f_stat = grs_stat,
    df1 = df1,
    df2 = df2,
    p_value = stats::pf(
      grs_stat,
      df1 = df1,
      df2 = df2,
      lower.tail = FALSE
    ),
    n_months = time_count
  )
}

panel_b_f_test <- table3_regression_data %>%
  group_by(subsample) %>%
  group_modify(~ grs_test_one_subsample(.x)) %>%
  ungroup()

## --------------------------------------------------------------------------
## 9. Final formatting and export
## --------------------------------------------------------------------------

panel_a_output <- panel_a %>%
  mutate(
    estimate = round(estimate, 4L),
    t_stat = round(t_stat, 2L)
  )

panel_b_output <- panel_b %>%
  mutate(
    alpha = round(alpha, 4L),
    alpha_t = round(alpha_t, 2L),
    market_beta = round(market_beta, 4L)
  )

panel_a_f_output <- panel_a_f_test %>%
  mutate(
    f_stat = round(f_stat, 4L),
    p_value = round(p_value, 6L)
  )

panel_b_f_output <- panel_b_f_test %>%
  mutate(
    f_stat = round(f_stat, 4L),
    p_value = round(p_value, 6L)
  )

write_xlsx(
  list(
    Panel_A = panel_a_output,
    Panel_A_F = panel_a_f_output,
    Panel_B = panel_b_output,
    Panel_B_F = panel_b_f_output,
    Panel_Check = table3_panel_check,
    Bad_Cohorts = table3_bad_cohorts,
    Formation_Check = formation_sort_check,
    SW_Beta_Diagnostics = stock_beta_year
  ),
  "results/table3.xlsx"
)

rm(
  list = intersect(
    c(
      "fama_duplicate_keys",
      "formation_sort_check",
      "market_model_result",
      "mean_t_result",
      "panel_a_joint_test",
      "grs_test_one_subsample",
      "table3_regression_data",
      "panel_b_spread_data",
      "panel_a_portfolios",
      "panel_a_spread",
      "panel_b_portfolios",
      "panel_b_spread"
    ),
    ls()
  )
)

# Core objects intentionally retained for inspection:
#   stock_beta_year
#   table3_cohort
#   table3_portfolio_ret
#   table3_panel_check
#   table3_bad_cohorts
#   panel_a_output
#   panel_a_f_output
#   panel_b_output
#   panel_b_f_output
