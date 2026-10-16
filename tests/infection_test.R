args <- commandArgs(FALSE)
path <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(path), ".."), winslash = "/")
setwd(APP_ROOT)
source("app.R", local = TRUE)
fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
cfg <- infection_config()
flat <- data.frame(time = c(0, 12, 24), concentration = 8)
stopifnot(infection_index(flat, 1, cfg, c(0, 12)) == 100)
cfg$fu <- .5
stopifnot(infection_index(flat, 4, cfg, c(0, 12)) == 0)
cfg$metric <- "auc"
stopifnot(infection_index(flat, 2, cfg, c(0, 24)) == 48, is.na(infection_index(flat, 2, cfg, c(0, 12))))
cfg$basis <- "total"
stopifnot(infection_index(flat, 2, cfg, c(0, 24)) == 96)
cfg$metric <- "peak"
stopifnot(infection_index(flat, 2, cfg, c(0, 12)) == 4)
cfg$metric <- "time"
sloped <- data.frame(time = c(0, 4), concentration = c(0, 8))
stopifnot(infection_index(sloped, 4, cfg, c(0, 4)) == 50)
bolus <- data.frame(time = c(0, 0, 4, 4, 8), concentration = c(0, 8, 0, 8, 0))
stopifnot(infection_index(bolus, 4, cfg, c(0, 8)) == 50)
fails(infection_index(flat, 0, cfg, c(0, 12)))
fails(infection_config(list(fu = 1.01)))
fails(infection_config(list(mic = ">8")))
fails(infection_config(list(replicates = 50.5)))
fails(infection_config(list(infusion = 25)))
fails(infection_config(list(interval = 168, delta = .01, replicates = 1000)))

html <- '<html><select id="search_antibiotic"><option value="-1">Select</option><option value="2">Synthetic</option></select><table id="search-results-table"><thead><tr><th></th><th>0.5</th><th>1</th><th>2</th><th>Distributions</th><th>Observations</th><th>(T)ECOFF</th><th>Confidence interval</th></tr></thead><tbody><tr><td><a href="/search/show-registration/123?back=x">Synthetic species</a></td><td>2</td><td>5</td><td>3</td><td>3</td><td>10</td><td>(1)</td><td>0.5 - 2</td></tr></tbody></table></html>'
document <- xml2::read_html(html)
stopifnot(unname(eucast_choices(document)) == "2")
parsed <- eucast_parse(document)$rows[[1]]
stopifnot(parsed$observations == 10, parsed$ecoff == "(1)", sum(parsed$data$count) == 10)
fails(eucast_parse(xml2::read_html(sub('<td>10</td>', '<td>11</td>', html, fixed = TRUE))))
fails(eucast_parse(xml2::read_html('<html>Maintenance</html>')))
fails(eucast_get("https://example.org"))
stopifnot(length(eucast_load("2", "Synthetic", get = function(...) document)$rows) == 1)

workdir <- tempfile("infection-test-"); dir.create(workdir)
code <- '$PARAM TVV=10, TVCL=1, ETA1=0, ETA2=0\n$CMT CENT\n$OMEGA 0.09 0.09\n$SIGMA 0.1\n$MAIN\ndouble V=TVV*exp(ETA1+ETA(1));\ndouble CL=TVCL*exp(ETA2+ETA(2));\n$ODE\ndxdt_CENT=-CL/V*CENT;\n$TABLE\ndouble CP=CENT/V;\ndouble DV=CP*(1+EPS(1));\n$CAPTURE CP DV'
model <- mrgsolve::mcode("infection_fixture", code, soloc = workdir, quiet = TRUE)
fit <- list(id = "synthetic", label = "Synthetic", model = model, estimate = NULL, contract = list(adm_cmt = 1L), current_covariates = list())
cfg <- infection_config(list(metric = "auc", target = 200, dose = 100, dose2 = 200, infusion = 0, infusion2 = 0,
  interval = 12, interval2 = 12, replicates = 50, delta = .1, basis = "total"))
a <- infection_simulate(list(synthetic = fit), c(synthetic = 1), cfg, "IV")
b <- infection_simulate(list(synthetic = fit), c(synthetic = 1), cfg, "IV")
stopifnot(identical(a$curves, b$curves), !a$posterior,
  all(diff(a$curves$pta[a$curves$regimen == "1"]) <= 0),
  all(a$curves$pta[a$curves$regimen == "2"] >= a$curves$pta[a$curves$regimen == "1"]))
cfg$dose2 <- cfg$dose
same <- infection_simulate(list(synthetic = fit), c(synthetic = 1), cfg, "IV")
stopifnot(identical(same$curves$pta[same$curves$regimen == "1"], same$curves$pta[same$curves$regimen == "2"]))
dayfit <- fit; dayfit$model <- mrgsolve::param(model, TVCL = 24); dayfit$hours_per_model_time <- 24
day <- infection_simulate(list(synthetic = dayfit), c(synthetic = 1), cfg, "IV")
stopifnot(max(abs(day$exposure$median - same$exposure$median)) < 1e-4)
fixed <- fit; fixed$model <- mrgsolve::zero_re(model)
fails(infection_simulate(list(synthetic = fixed), c(synthetic = 1), cfg, "IV"))
fails(infection_simulate(list(synthetic = fit), c(synthetic = NA), cfg, "IV"))
bad <- cfg; bad$infusion <- 1
fails(infection_simulate(list(synthetic = fit), c(synthetic = 1), bad, "Oral"))
continuous <- fit; continuous$administration_modes <- "IV_CONTINUOUS"
fails(infection_simulate(list(synthetic = continuous), c(synthetic = 1), cfg, "IV"))
continuous_cfg <- cfg; continuous_cfg$infusion <- continuous_cfg$interval; continuous_cfg$infusion2 <- continuous_cfg$interval2
stopifnot(nrow(infection_simulate(list(synthetic = continuous), c(synthetic = 1), continuous_cfg, "IV")$curves) > 0)
posterior <- fit
posterior$estimate <- list(final_eta = list(c(ETA1 = 0, ETA2 = 0)), covariance = list(diag(.01, 2)))
post <- infection_simulate(list(synthetic = posterior), c(synthetic = 1), cfg, "IV")
stopifnot(post$posterior)
badpost <- posterior; badpost$estimate$covariance <- list(NULL)
fails(infection_simulate(list(synthetic = badpost), c(synthetic = 1), cfg, "IV"))
mix <- infection_simulate(list(first = fit, second = fit), c(first = .25, second = .75), cfg, "IV")
stopifnot(all(mix$curves$n == 50), identical(unname(mix$weights), c(.25, .75)))
payload <- list(type = "pk-workbench", id = "infection-test", version = 1L, view = "infection", config = cfg,
  models = list(list(source = "library", id = "vanco_pkjust", route = "IV")))
stopifnot(workbench_payload(payload)$view == "infection")
grid_cfg <- cfg
grid_cfg$grid <- list(min = 100, max = 200, step = 100, intervals = c(8,12), infusion = 0)
gr <- infection_simulate(list(synthetic = fit), c(synthetic = 1), grid_cfg, "IV")
ranked <- infection_rank(gr)
stopifnot(nrow(ranked) == 5, ranked$regimen[1] == "current", all(is.finite(ranked$pta)),
  all(ranked$target_met == (ranked$pta >= grid_cfg$pta_target)))
same_row <- gr$regimens$regimen[gr$regimens$dose == grid_cfg$dose & gr$regimens$interval == grid_cfg$interval]
stopifnot(length(same_row) == 2, identical(gr$curves$pta[gr$curves$regimen == same_row[1]], gr$curves$pta[gr$curves$regimen == same_row[2]]))
badgrid <- grid_cfg; badgrid$grid$step <- .001
fails(infection_config(badgrid))
badgrid$grid <- grid_cfg$grid; badgrid$grid$intervals <- numeric()
fails(infection_config(badgrid))
badgrid$grid <- grid_cfg$grid; badgrid$grid$infusion <- 24
fails(infection_config(badgrid))
badgrid <- grid_cfg; badgrid$replicates <- 250; badgrid$delta <- .01
badgrid$grid <- list(min = 100, max = 100, step = 1, intervals = 168, infusion = 0)
fails(infection_config(badgrid))
grid_cfg$grid$continuous <- TRUE
stopifnot(all(infection_regimens(grid_cfg)$infusion[-1] == infection_regimens(grid_cfg)$interval[-1]))
bad <- payload; bad$models <- list()
fails(workbench_payload(bad))

store <- reactiveVal(list(fits = list(synthetic = posterior), weights = c(synthetic = 1), route = "IV"))
shiny::testServer(infection_server, args = list(analysis_store = store, report_plot_uri = function(...) "", soloc = workdir, cache = new.env()), {
  session$setInputs(source = "tdm", metric = "auc", target = 200, basis = "total", tdm_scale = 1, dose = 100, dose2 = 100,
    interval = 12, interval2 = 12, infusion = 0, infusion2 = 0, replicates = 50, delta = .1, accept = TRUE)
  session$setInputs(simulate = 1)
  stopifnot(!is.null(simulation()), simulation()$posterior)
  first <- simulation()
  session$setInputs(pk_concentration = "DV")
  stopifnot(identical(simulation(), first))
  session$setInputs(mic = 2)
  stopifnot(is.null(simulation()))
})
if ("--eucast-live" %in% commandArgs(TRUE)) {
  live <- eucast_load("2", "Amoxicillin")
  stopifnot("Escherichia coli" %in% names(live$rows), sum(live$rows[["Escherichia coli"]]$data$count) == live$rows[["Escherichia coli"]]$observations)
  cat("Live EUCAST lookup verified; no distribution saved.\n")
}
for (dll in list.files(workdir, pattern = "\\.dll$", recursive = TRUE, full.names = TRUE)) try(dyn.unload(dll), silent = TRUE)
unlink(workdir, recursive = TRUE)
cat("Infectiology: indices, units, seeded population/posterior PTA, averaging, invalidation, bridge and EUCAST parser passed.\n")
