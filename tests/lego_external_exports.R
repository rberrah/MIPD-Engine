suppressPackageStartupMessages({library(jsonlite); library(mrgsolve)})
cases <- fromJSON("test-results/lego-external-exports.json", simplifyVector=FALSE)
work <- tempfile("lego-external-")
dir.create(work)
run_checks <- function() {
  on.exit(unlink(work, recursive=TRUE), add=TRUE)
  for (id in unique(vapply(cases, `[[`, character(1), "id"))) {
    variants <- Filter(function(x) x$id == id, cases)
    baseline <- NULL
    for (case in variants) {
      mod <- mcode(paste0(id, "_", case$format), case$code, soloc=work, quiet=TRUE)
      spec <- case$spec
      dosed <- Filter(function(n) n$dose > 0, spec$nodes)
      cmt <- if ("LEGO_INPUT" %in% names(mod@init)) match("LEGO_INPUT", names(mod@init)) else match(dosed[[1]]$name, names(mod@init))
      horizon <- if (id %in% c("samtani", "tjollyn")) 5000 else 48
      times <- seq(0, horizon, length.out=501)
      events <- ev(amt=100, cmt=cmt)
      out <- as.data.frame(mrgsim(zero_re(mod), events=events, tgrid=times, obsonly=TRUE, recsort=3))
      stopifnot(all(is.finite(out$DV)), max(out$DV)>0)
      if (is.null(baseline)) baseline <- out$DV
      stopifnot(max(abs(out$DV-baseline)) < 1e-7 * max(1, abs(baseline)))
      if (id == "samtani") {
        t <- out$time
        ka <- 0.000488; cl <- 4.95; volume <- 391; duration <- 319; fraction <- 0.168
        k <- cl/volume
        fast <- 100*fraction/(duration*cl)*(1-exp(-k*pmin(t,duration)))*exp(-k*pmax(t-duration,0))
        delayed <- pmax(t-duration,0)
        slow <- 100*(1-fraction)*ka/(ka-k)*(exp(-k*delayed)-exp(-ka*delayed))/volume
        stopifnot(max(abs(out$DV-fast-slow)) < 1e-7)
      }
      if (!is.null(case$native) && case$format == "mrgsolve") {
        native <- mcode(paste0(id,"_native"), case$native, soloc=work, quiet=TRUE)
        for (comparison in c(FALSE,TRUE)) {
          covs <- setNames(lapply(spec$covariates,function(c) if(comparison) c$comparison else c$reference), vapply(spec$covariates, `[[`, character(1), "name"))
          covs <- covs[!duplicated(names(covs))]
          a <- param(mod,covs); b <- param(native,covs)
          x <- as.data.frame(mrgsim(zero_re(a),events=events,tgrid=times,obsonly=TRUE,recsort=3))
          y <- as.data.frame(mrgsim(zero_re(b),events=ev(amt=100,cmt=1),tgrid=times,obsonly=TRUE,recsort=3))
          stopifnot(max(abs(x$DV-y$DV)) < 1e-6 * max(1,abs(y$DV)))
        }
      }
      message("PASS external native simulation ",id," ",case$format)
    }
  }
}
run_checks()
