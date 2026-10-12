## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE,include=FALSE------------------------
library(postshock)


## ----load-cop-data------------------------------------------------------------
data("cop_oil_postshock", package = "postshock")

names(cop_oil_postshock)

## ----inspect-cop-data---------------------------------------------------------
cop_structure <- data.frame(
  Position = seq_along(cop_oil_postshock$Y),
  Series = names(cop_oil_postshock$Y),
  Covariates = names(cop_oil_postshock$X),
  Shock_time = cop_oil_postshock$shock_time_vec,
  Shock_length = cop_oil_postshock$shock_length_vec
)

knitr::kable(
  cop_structure,
  caption = "Alignment of the target and donor inputs in `cop_oil_postshock`."
)

## ----fit-cop-dbw, results='hide', message=FALSE, warning=FALSE----------------
cop_dbw <- dbw(
  X = cop_oil_postshock$X,
  dbw_indices = seq_len(ncol(cop_oil_postshock$X[[1]])),
  shock_time_vec = cop_oil_postshock$shock_time_vec,
  center = TRUE,
  scale = TRUE,
  sum_to_1 = TRUE,
  bounded_below_by = 0,
  bounded_above_by = 1,
  normchoice = "l2",
  penalty_normchoice = "l2",
  penalty_lambda = 1e-2
)

## ----inspect-cop-dbw----------------------------------------------------------
names(cop_dbw)

cop_dbw$convergence
cop_dbw$loss
sum(cop_dbw$opt_params)

## ----display-cop-dbw----------------------------------------------------------
cop_weight_table <- data.frame(
  Donor = sub(
    "^Y_",
    "",
    names(cop_oil_postshock$Y)[-1]
  ),
  Weight = as.numeric(cop_dbw$opt_params)
)

knitr::kable(
  cop_weight_table,
  digits = 3,
  caption = "Donor balancing weights for the ConocoPhillips example."
)

## ----fit-cop-synthprediction, results='hide', message=FALSE, warning=FALSE----
cop_fit <- SynthPrediction(
  Y_series_list = cop_oil_postshock$Y,
  covariates_series_list = cop_oil_postshock$X,
  shock_time_vec = cop_oil_postshock$shock_time_vec,
  shock_length_vec = cop_oil_postshock$shock_length_vec,
  k = 1L,
  dbw_indices = seq_len(ncol(cop_oil_postshock$X[[1]])),
  dbw_center = TRUE,
  dbw_scale = TRUE,
  use_dbw = TRUE,
  covariate_indices = NULL,
  seasonal = TRUE,
  penalty_lambda = 1e-2,
  penalty_normchoice = "l2"
)

## ----inspect-cop-fit----------------------------------------------------------
names(cop_fit)
names(cop_fit$predictions)
names(cop_fit$meta)

## ----cop-forecast-comparison--------------------------------------------------
forecast_index <- cop_oil_postshock$shock_time_vec[1] + 1L

observed_value <- as.numeric(
  cop_oil_postshock$Y[[1]][forecast_index]
)

forecast_values <- c(
  as.numeric(cop_fit$predictions$unadjusted[1]),
  as.numeric(cop_fit$predictions$arithmetic_mean[1]),
  as.numeric(cop_fit$predictions$adjusted[1])
)

cop_forecast_table <- data.frame(
  Method = c(
    "Observed value",
    "Unadjusted",
    "Arithmetic mean",
    "DBW adjusted"
  ),
  Value = c(
    observed_value,
    forecast_values
  ),
  Absolute_error = c(
    NA_real_,
    abs(forecast_values - observed_value)
  )
)

knitr::kable(
  cop_forecast_table,
  digits = 3,
  col.names = c("Method", "Value", "Absolute error"),
  caption = paste(
    "One-step forecast comparison for the March 2020",
    "ConocoPhillips shock episode."
  )
)

## ----inspect-cop-donor-results------------------------------------------------
cop_donor_summary <- data.frame(
  Donor = sub(
    "^Y_",
    "",
    names(cop_oil_postshock$Y)[-1]
  ),
  Weight = as.numeric(cop_fit$linear_combinations),
  Estimated_effect = as.numeric(cop_fit$meta$omega_vec)
)

knitr::kable(
  cop_donor_summary,
  digits = 3,
  col.names = c(
    "Donor",
    "DBW weight",
    "Estimated post-shock effect"
  ),
  caption = "Donor weights and estimated effects for the COP example."
)

## ----inspect-cop-adjustment---------------------------------------------------
data.frame(
  Unadjusted_forecast =
    as.numeric(cop_fit$predictions$unadjusted[1]),
  Combined_donor_effect =
    as.numeric(cop_fit$meta$combined_omega),
  Adjusted_forecast =
    as.numeric(cop_fit$predictions$adjusted[1])
) |>
  knitr::kable(
    digits = 3,
    col.names = c(
      "Unadjusted forecast",
      "Combined donor effect",
      "DBW-adjusted forecast"
    ),
    caption = "Construction of the DBW-adjusted COP forecast."
  )

## ----define-cop-pools---------------------------------------------------------
cop_pools <- list(
  supply = "Y_2014_11_28",
  demand = c(
    "Y_2008_03_17",
    "Y_2008_09_09",
    "Y_2008_09_15",
    "Y_2008_09_29"
  )
)

cop_pool_table <- data.frame(
  Pool = names(cop_pools),
  Donors = unname(
    vapply(
      cop_pools,
      function(x) {
        paste(
          gsub("_", "-", sub("^Y_", "", x)),
          collapse = ", "
        )
      },
      character(1)
    )
  )
)

knitr::kable(
  cop_pool_table,
  row.names = FALSE,
  col.names = c("Pool", "Donor episodes"),
  caption = "Structured donor pools for the ConocoPhillips example."
)

## ----fit-cop-multipool, results='hide', message=FALSE, warning=FALSE----------
cop_multipool <- SynthPrediction(
  Y_series_list = cop_oil_postshock$Y,
  covariates_series_list = cop_oil_postshock$X,
  shock_time_vec = cop_oil_postshock$shock_time_vec,
  shock_length_vec = cop_oil_postshock$shock_length_vec,
  k = 1L,
  covariate_indices = NULL,
  seasonal = TRUE,
  dbw_indices = seq_len(
    ncol(cop_oil_postshock$X[[1]])
  ),
  dbw_center = TRUE,
  dbw_scale = TRUE,
  penalty_lambda = 1e-2,
  penalty_normchoice = "l2",
  pools = cop_pools,
  pool_agg = "sum",
  base_use_dbw = FALSE,
  pool_use_dbw = TRUE
)

## ----inspect-cop-pool-effects-------------------------------------------------
cop_pool_effects <- data.frame(
  Pool = names(
    cop_multipool$meta$multi_pool$pool_effects
  ),
  Estimated_adjustment = as.numeric(
    cop_multipool$meta$multi_pool$pool_effects
  )
)

knitr::kable(
  cop_pool_effects,
  row.names = FALSE,
  digits = 3,
  col.names = c(
    "Pool",
    "Estimated pool adjustment"
  ),
  caption = "Pool-specific adjustments for the ConocoPhillips example."
)

## ----compare-cop-multipool----------------------------------------------------
cop_multipool_values <- c(
  as.numeric(
    cop_multipool$predictions$unadjusted[1]
  ),
  as.numeric(
    cop_fit$predictions$adjusted[1]
  ),
  as.numeric(
    cop_multipool$predictions$adjusted_multipool[1]
  )
)

cop_multipool_table <- data.frame(
  Method = c(
    "Observed value",
    "Unadjusted",
    "Single-pool DBW",
    "MultiPool (sum)"
  ),
  Value = c(
    observed_value,
    cop_multipool_values
  ),
  Absolute_error = c(
    NA_real_,
    abs(cop_multipool_values - observed_value)
  )
)

knitr::kable(
  cop_multipool_table,
  digits = 3,
  col.names = c("Method", "Value", "Absolute error"),
  caption = paste(
    "Comparison of single-pool and MultiPool forecasts",
    "for the ConocoPhillips example."
  )
)

## ----load-iyg-data------------------------------------------------------------
data("iyg_onestep_input", package = "postshock")

names(iyg_onestep_input)
names(iyg_onestep_input$inputs)

## ----inspect-iyg-input--------------------------------------------------------
iyg_restricted <- iyg_onestep_input$inputs$restricted

names(iyg_restricted)

## ----summarize-iyg-input------------------------------------------------------
iyg_restricted$target_event
iyg_restricted$donor_events

c(
  number_of_series =
    length(iyg_restricted$Y_series_list),
  number_of_donors =
    length(iyg_restricted$donor_events),
  donor_covariate_sets =
    length(iyg_restricted$covariates_series_list),
  shock_time_entries =
    length(iyg_restricted$shock_time_vec),
  shock_length_entries =
    length(iyg_restricted$shock_length_vec)
)

## ----fit-iyg-volatility, results='hide', message=FALSE, warning=FALSE---------
iyg_fit <- SynthVolForecast(
  Y_series_list =
    iyg_restricted$Y_series_list,
  covariates_series_list =
    iyg_restricted$covariates_series_list,
  target_covariates =
    iyg_restricted$target_covariates,
  shock_time_vec =
    iyg_restricted$shock_time_vec,
  shock_length_vec =
    iyg_restricted$shock_length_vec,
  k = 1L,
  dbw_indices = seq_len(
    ncol(iyg_restricted$target_covariates)
  ),
  dbw_center = TRUE,
  dbw_scale = TRUE,
  penalty_lambda = 1e-2,
  penalty_normchoice = "l2",
  covariate_indices = NULL,
  garch_order = NULL,
  max.p = 3L,
  max.q = 3L,
  vcov.type = "robust",
  return_fits = FALSE
)

## ----inspect-iyg-orders-------------------------------------------------------
iyg_selected_orders <- do.call(
  rbind,
  c(
    iyg_fit$meta$donor_orders,
    list(iyg_fit$meta$target_order)
  )
)

iyg_order_table <- data.frame(
  Model = c(
    paste("Donor:", iyg_restricted$donor_events),
    paste("Target:", iyg_restricted$target_event)
  ),
  GARCH_q = iyg_selected_orders[, 1],
  ARCH_p = iyg_selected_orders[, 2],
  Asymmetry_o = iyg_selected_orders[, 3]
)

knitr::kable(
  iyg_order_table,
  row.names = FALSE,
  col.names = c(
    "Model",
    "GARCH order (q)",
    "ARCH order (p)",
    "Asymmetry order (o)"
  ),
  caption = paste(
    "Orders selected automatically by BIC for the",
    "donor and target variance models."
  )
)

## ----inspect-iyg-donor-adjustments--------------------------------------------
iyg_donor_summary <- data.frame(
  Donor = as.character(iyg_restricted$donor_events),
  Weight = as.numeric(
    iyg_fit$linear_combinations
  ),
  Estimated_adjustment = as.numeric(
    iyg_fit$meta$omega_vec
  )
)

knitr::kable(
  iyg_donor_summary,
  row.names = FALSE,
  digits = 6,
  col.names = c(
    "Donor episode",
    "DBW weight",
    "Estimated variance adjustment"
  ),
  caption = paste(
    "Donor weights and estimated variance adjustments",
    "for the restricted IYG specification."
  )
)

## ----inspect-iyg-combined-adjustment------------------------------------------
iyg_adjustment_table <- data.frame(
  Unadjusted_forecast = as.numeric(
    iyg_fit$predictions$unadjusted[1]
  ),
  Combined_adjustment = as.numeric(
    iyg_fit$meta$combined_omega
  ),
  Adjusted_forecast = as.numeric(
    iyg_fit$predictions$adjusted[1]
  )
)

knitr::kable(
  iyg_adjustment_table,
  row.names = FALSE,
  digits = 6,
  col.names = c(
    "Unadjusted variance forecast",
    "Combined donor adjustment",
    "DBW-adjusted variance forecast"
  ),
  caption = "Construction of the DBW-adjusted IYG variance forecast."
)

## ----compare-iyg-variance-forecasts-------------------------------------------
iyg_proxy_index <-
  iyg_restricted$shock_time_vec[1] + 1L

iyg_observed_proxy <- as.numeric(
  iyg_restricted$Y_series_list[[1]][iyg_proxy_index]
)^2

iyg_forecast_values <- c(
  as.numeric(
    iyg_fit$predictions$unadjusted[1]
  ),
  as.numeric(
    iyg_fit$predictions$arithmetic_mean[1]
  ),
  as.numeric(
    iyg_fit$predictions$adjusted[1]
  )
)

iyg_forecast_table <- data.frame(
  Method = c(
    "Observed squared-return proxy",
    "Unadjusted",
    "Arithmetic mean",
    "DBW adjusted"
  ),
  Variance = c(
    iyg_observed_proxy,
    iyg_forecast_values
  ),
  Absolute_error = c(
    NA_real_,
    abs(
      iyg_forecast_values -
        iyg_observed_proxy
    )
  )
)

knitr::kable(
  iyg_forecast_table,
  row.names = FALSE,
  digits = 6,
  col.names = c(
    "Method",
    "Variance",
    "Absolute error"
  ),
  caption = paste(
    "One-step variance forecast comparison for IYG",
    "under the restricted donor-pool specification."
  )
)

