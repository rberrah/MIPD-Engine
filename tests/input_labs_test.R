suppressPackageStartupMessages({ library(mrgsolve); library(jsonlite) })
APP_ROOT <- normalizePath("tdm-engine", winslash = "/", mustWork = TRUE)
source("tdm-engine/R/model_library.R")
run <- function() {
  directory <- tempfile("pk-input-lab-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fixture <- file.path(directory, "synthetic.json")
  stopifnot(system2("node", c("scripts/test_input_labs.mjs", shQuote(fixture))) == 0)
  fixtures <- jsonlite::fromJSON(fixture, simplifyVector = FALSE)
  reference <- mrgsolve::mcode("input_reference", '
$PARAM CL=5, V=25, KA=0.8, BIO=0.7
$CMT CENT GUT LOSS ELIM
$ODE
dxdt_GUT = -KA*GUT;
dxdt_CENT = BIO*KA*GUT - CL*CENT/V;
dxdt_LOSS = (1-BIO)*KA*GUT;
dxdt_ELIM = CL*CENT/V;
$TABLE
capture CP = CENT/V;
', soloc = directory, quiet = TRUE)
  for (i in seq_along(fixtures)) {
    f <- fixtures[[i]]; p <- f$parameters
    expected <- do.call(rbind, lapply(f$points, as.data.frame))
    oral <- f$lab == "absorption"
    mod <- mrgsolve::param(reference, CL = p$cl, V = p$vc, KA = if (oral) p$ka else .8, BIO = if (oral) p$f else 1)
    event <- mrgsolve::ev(amt = p$dose, cmt = if (oral) 2 else 1, rate = if (oral) 0 else p$dose/p$duration)
    result <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(mod,event), tgrid=expected$t, obsonly=TRUE, recsort=3, atol=1e-11, rtol=1e-11))
    stopifnot(max(abs(result$CP-expected$c)) < 1e-6,
              max(abs(result$ELIM-expected$eliminated)) < 1e-6,
              max(abs(result$LOSS-expected$lost)) < 1e-6)
    # Verify the actual Lego export as well, including its routed infusion event.
    code <- safe_lego_model_code(specification = f$spec)
    exported <- mrgsolve::zero_re(mrgsolve::mcode(paste0("input_lego",i),code,soloc=directory,quiet=TRUE))
    result <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(exported,mrgsolve::ev(amt=p$dose,cmt=if(oral) 2 else 1)),tgrid=expected$t,obsonly=TRUE,recsort=3,atol=1e-11,rtol=1e-11))
    stopifnot(max(abs(result$DV-expected$c)) < 1e-6)
  }
  message(length(fixtures), " input scenarios: independent ODE and Lego-generated mrgsolve comparisons passed.")
}
run()
