# MAP-BE position and model-aware diagnostics: deterministic, no patient files.
suppressPackageStartupMessages({
  library(mrgsolve)
  library(mapbayr)
  library(ggplot2)
  library(jsonlite)
})
`%||%` <- function(x,y) if (is.null(x) || !length(x)) y else x
APP_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = FALSE)
source(file.path(APP_ROOT, "R", "engine.R"), local = FALSE)
source(file.path(APP_ROOT, "R", "map_be_diagnostics.R"), local = FALSE)

omega <- matrix(c(.25, .06, .06, .09), 2L, 2L, dimnames=list(c("CL", "V"), c("CL", "V")))
rows <- map_be_eta_position(c(.5, -.3), omega, labels=rownames(omega),
                            selected=c(1L, 2L), posterior_covariance=diag(c(.01, .04)))
stopifnot(nrow(rows)==2L, identical(rows$parameter,c("CL","V")),
          isTRUE(all.equal(rows$prior_sd,c(.5,.3))),
          isTRUE(all.equal(rows$z,c(1,-1))),
          isTRUE(all.equal(rows$percentile,100*stats::pnorm(c(1,-1)))),
          isTRUE(all.equal(rows$posterior_sd,c(.1,.2))),
          all(rows$posterior_95>rows$eta),
          all(rows$posterior_5<rows$eta))
# Non-diagonal Omega: percentile is explicitly a *marginal* reference.
stopifnot(isTRUE(all.equal(rows$percentile[[1]],100*stats::pnorm(1))))
# Supplied selected indices are not equivalent to all model parameters.
partial <- map_be_eta_position(c(.5), omega, c("CL", "V"), selected=1L)
stopifnot(nrow(partial)==1L, identical(partial$parameter,"CL"))
fixed <- omega; fixed[2,]<-0; fixed[,2]<-0
stopifnot(nrow(map_be_eta_position(c(.5,-.3),fixed,c("CL","V")))==1L)
stopifnot(is.null(map_be_eta_position(c(.1),omega,c("CL","V"),selected=c(1L,2L))))
invalid <- omega; invalid[1,2] <- .2
stopifnot(is.null(map_be_eta_position(c(.1,.2),invalid)))
bad_cov <- matrix(c(1, 2, 2, 1),2,2)
no_ci <- map_be_eta_position(c(.5,-.3),omega,posterior_covariance=bad_cov)
stopifnot(all(is.na(no_ci$posterior_sd)), all(is.na(no_ci$posterior_5)))

prior <- data.frame(ke=c(.1,.2,.3),v=c(30,40,50),prob=c(.2,.5,.3))
posterior <- prior; posterior$prob <- c(.05,.9,.05)
support <- map_be_pmetrics_position(prior,posterior,c("ke","v"))
stopifnot(nrow(support$table)==2L,nrow(support$points)==12L,
          isTRUE(all.equal(support$table$patient[[1]],.2)),
          identical(unique(support$points$distribution), c("Prior", "Posterior")),
          isTRUE(all.equal(map_be_weighted_quantile(prior$ke,prior$prob,.8),.3)))
stopifnot(is.null(map_be_pmetrics_position(prior,
  data.frame(prob=c(.4,.6)),c("ke"))))
stopifnot(inherits(map_be_diagnostic_plot(list(type="mapbayr",table=rows)),"ggplot"))
stopifnot(inherits(map_be_diagnostic_plot(c(list(type="pmetrics"),support)),"ggplot"))

# Integrated mapbayr fit: map estimated ETA to *real* model OMEGA.
record <- model_record("vanco_roberts")
spec <- list(id="vanco_roberts", label="Roberts 2011", code=NULL,
  route="IV", mode=model_administration_mode(record,"IV","IV_CONTINUOUS"),
  adm_cmt_name=model_administration_cmt(record,"IV"))
doses <- data.frame(time=0,amount=1000,interval=12,count=4L,
                    infusion=12,ss=0L,time_uncertainty=0)
obs <- data.frame(time=47.5,concentration=20,time_uncertainty=0)
fit <- fit_one_model(spec,doses,obs,
  covariates=list(CRCL=90,WT=74),allow_custom=FALSE)
stopifnot(!is.null(fit$estimate))
data <- map_be_diagnostic(fit)
stopifnot(identical(data$type,"mapbayr"),
          nrow(data$table)>0L,
          nrow(data$table)<=nrow(as.matrix(mrgsolve::omat(fit$model))),
          all(data$table$percentile>0 & data$table$percentile<100),
          all(is.finite(data$table$eta)),
          all(is.finite(data$table$prior_sd)))
without_observations <- fit
without_observations$estimate <- NULL
stopifnot(identical(map_be_diagnostic(without_observations)$reason,"population"))
# The view must use genuine MAP ETAs and must never use experimental ML overrides.
fit$ml_eta_override <- stats::setNames(rep(100,nrow(data$table)), paste0("ETA",seq_len(nrow(data$table))))
stopifnot(isTRUE(all.equal(map_be_diagnostic(fit)$table$eta,data$table$eta)))

app <- paste(readLines(file.path(APP_ROOT,"app.R"),encoding="UTF-8",warn=FALSE),collapse="\n")
pos <- regexpr('              "Posologies",',app,fixed=TRUE)[[1]]
map <- regexpr('              "MAP-BE",',app,fixed=TRUE)[[1]]
model <- regexpr('              "Modèle",\n              value = "model",',app,fixed=TRUE)[[1]]
stopifnot(pos>0, map>pos, model>map,
          grepl('output$map_be_parameter_plot <- renderPlot',app,fixed=TRUE),
          grepl('output$map_be_parameter_table <- DT::renderDT',app,fixed=TRUE),
          grepl('output$download_map_be_parameters <- downloadHandler',app,fixed=TRUE),
          grepl('source(file.path(APP_ROOT, "R", "map_be_diagnostics.R")',app,fixed=TRUE))
cat("MAP-BE model-specific marginal population positioning, uncertainty guards and Pmetrics support: PASS\n")
