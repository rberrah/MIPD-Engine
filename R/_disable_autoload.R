# Shiny >= 1.5 sources all R/ scripts alphabetically before executing app.R.
# MIPD's modules have explicit dependency ordering in app.R (model catalogue,
# TDM drug profiles, then runtime presets). Autoloading this directory would
# run drug_profile_runtime_presets.R before tdm_drug_profiles.R and crash
# production before the first session, producing an HTTP 500.
#
# Shiny recognizes this exact filename and skips implicit R/ autoloading.
# app.R explicitly sources every required module, including the appropriate
# enabled/disabled Pmetrics backend. Do not remove this marker without
# implementing and testing a dependency-safe autoload architecture.
