args <- commandArgs(FALSE)
path <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(path), ".."), winslash = "/")
setwd(APP_ROOT)
source("app.R", local = TRUE)
workdir <- tempfile("pd-pk-test-"); dir.create(workdir)
code <- "$PARAM KA=1, CL=2, V1=10, Q=3, V2=20\n$CMT DEPOT CENT PERI\n$ODE\ndxdt_DEPOT=-KA*DEPOT;\ndxdt_CENT=KA*DEPOT-(CL+Q)*CENT/V1+Q*PERI/V2;\ndxdt_PERI=Q*CENT/V1-Q*PERI/V2;\n$TABLE double CP=CENT/V1;\n$CAPTURE CP"
model <- mrgsolve::mcode("oral_two_cmt", code, soloc = workdir, quiet = TRUE)
context <- list(model = model, code = code, route = "Oral", adm_cmt = 1, concentration = "CP", concentration_scale = 1, time_unit = "h", label = "Synthetic oral two compartment")
doses <- data.frame(time = c(0, 24), amount = 100, infusion = 0)
profile <- pd_pk_profile(context, doses, 48)
stopifnot(min(profile$time) == 0, max(profile$time) == 48, profile$concentration[1] == 0, max(profile$concentration) > 3)
day_context <- context
day_context$model <- mrgsolve::param(model, KA = 24, CL = 48, Q = 72)
day_context$time_unit <- "day"
daily <- pd_pk_profile(day_context, doses, 48)
stopifnot(max(abs(profile$concentration - approx(daily$time, daily$concentration, profile$time)$y)) < 1e-4)
scaled_context <- context; scaled_context$concentration_scale <- 0.001
scaled <- pd_pk_profile(scaled_context, doses, 48)
stopifnot(max(abs(profile$concentration / 1000 - scaled$concentration)) < 1e-10)
bad <- doses; bad$infusion <- 1
stopifnot(inherits(try(pd_pk_profile(context, bad, 48), silent = TRUE), "try-error"))
config <- onco_config(list(free_pk = TRUE, decision = 2, horizon = 6, interval = 2, infusion = 0, history = data.frame(time = 0, amount = 100, infusion = 0)))
result <- onco_compare(config, pk_context = context)
stopifnot(all(c("untreated", "maintain", "change") %in% names(result$curves)), all(result$curves$untreated$concentration == 0))
stopifnot(tail(result$curves$maintain$TUMOR, 1) < tail(result$curves$untreated$TUMOR, 1))
past <- result$curves$maintain$time < config$decision
stopifnot(max(abs(result$curves$maintain$TUMOR[past] - approx(result$curves$change$time, result$curves$change$TUMOR, result$curves$maintain$time[past])$y)) < 1e-4)
empty <- pd_pk_profile(context, doses[FALSE, ], 48)
stopifnot(all(empty$concentration == 0))
json <- list(type = "pk-workbench", version = 1, id = "test-file", view = "onco", config = config,
  models = list(list(source = "code", code = code, route = "Oral", time_unit = "h", concentration_scale = 1)))
stopifnot(workbench_payload(json)$config$free_pk)
stopifnot(grepl("PD ONLY", onco_model_code(config)), grepl("mrgsolve::mcode", onco_export_script(config, pk_context = context)))
imported <- shiny::reactiveVal(NULL)
original_compile <- compile_model
compile_model <- function(...) model
shiny::testServer(pd_pk_server, args = list(soloc = workdir, cache = new.env(), imported = imported), {
  session$setInputs(source = "code", code = code, custom_route = "Oral", time_unit = "h", scale = 1, load = 1)
  session$setInputs(param_1_CL = 9)
  stopifnot(as.list(mrgsolve::param(context()$model))$CL == 9)
  imported(list(source = "code", code = code, route = "Oral", concentration_scale = .001)); session$flushReact()
  stopifnot(generation() == 2L, as.list(mrgsolve::param(context()$model))$CL == 2)
  stopifnot(grepl('value="0.001" selected', output$units$html, fixed = TRUE))
  session$setInputs(scale = .001)
  session$setInputs(load = 2)
  stopifnot(generation() == 3L, as.list(mrgsolve::param(context()$model))$CL == 2, context()$concentration_scale == .001, is.null(pending()))
})
compile_model <- original_compile
for (dll in list.files(workdir, pattern = "\\.dll$", recursive = TRUE, full.names = TRUE)) try(dyn.unload(dll), silent = TRUE)
unlink(workdir, recursive = TRUE)
cat("Free PK: oral two compartments, time/concentration conversion, oncology coupling, untreated control and JSON validation passed.\n")
