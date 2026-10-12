# postshock 0.2.1

## Fixes

* Forward explicit PCA settings from `SynthPrediction()` to `dbw()`,
  including the MultiPool dispatcher. Validate dimensions against the
  available SVD directions. Explicit PCA failures, unsuccessful or unknown
  DBW status, and invalid weights now produce errors instead of silently
  substituting equal weights; the default NULL and disabled-DBW paths remain.
* Forward all shared model, matching, penalty, seasonal and control settings
  through the MultiPool dispatcher. Reject conflicting settings while
  preserving separate baseline/pool DBW switches and unnormalized pool weights.
* Require complete, finite, aligned inputs when additional donor ARIMAX
  covariates are selected. Report invalid rows and available dates rather
  than deleting covariate rows and misaligning responses. Other input-policy
  branches retain their existing behavior.
* Restrict automatic GARCH-X selection to fits with numerical convergence
  code zero, finite nonempty coefficients and finite BIC. Preserve the first
  eligible fit on exact BIC ties; report an error if none qualifies.
  Backcast refits must meet the same conditions and strictly improve BIC.
  Non-finite standard errors alone do not disqualify a point fit.
* Pass the specified or detected per-window period to automatic ARIMA's
  response frequency, with strict validation of explicit periods. Retain
  the non-seasonal ARIMA(1,1,1) donor failure fallback.
* Correct and align documentation for GARCH-X order `c(q, p, o)`:
  `q` is GARCH, `p` is ARCH, and `o` is the asymmetry order.

## Reproducibility

* Add `aapl_s1_input`, derived from the unchanged frozen Apple objects,
  with evaluation observations separate from fitting and matching inputs.
* Install the Apple input builder, shared analysis and version/function-body
  identity checks under `replication/`. The Technical Report uses these
  installed resources and executable results without an adjacent manuscript.
* Hash replication function bodies as UTF-8 bytes with fixed LF line endings,
  avoiding false implementation-mismatch errors on Windows. Retain the
  version, implementation-identity and numerical checks.

## Compatibility

Public function signatures and existing frozen data are unchanged. Previously
ignored explicit settings may now change forecasts; invalid donor ARIMAX inputs,
failed explicit PCA requests and ineligible automatic GARCH-X fits can now
produce errors. These corrections preserve the article's existing main
application specifications.
