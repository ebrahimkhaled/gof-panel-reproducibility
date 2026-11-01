# run_extra.R -- the extra simulations asked for by the referee pass of the panel paper.
#
#   C  second corruption mechanisms (n = 500; n = 1000 dropped: about 70 s per data set): the correct model of the archive, with k = 1, 5, 10 recorded
#      records corrupted in one of two new ways: the response is flipped (misclassified y), or x2 is multiplied by 4.
#      The clean data sets are the same as in run_corrupt.R (same seeds), so k = 0 is corrupt/corrupt_k0_n*.csv.
#   A  two more null designs (n = 200, 500): "rare" (same covariates, intercept -2.6, about 13% events) and "six"
#      (six covariates, one of them binary). For each: 2,000 reference data sets (to calibrate) and 2,000 test data
#      sets (to measure the false-alarm rate).
#   B  calibration as an analyst would do it: for each data set, a parametric bootstrap from the FITTED model
#      (B replicates, refit, all 14 tests) gives the bootstrap p-value of four combination rules on two panels.
#      Size: null model, n = 200 and 500, 500 data sets, B = 199. Power: three alternatives at n = 500,
#      200 data sets each, B = 99 (the observed p-values are stored too, so the fixed reference calibration of the
#      paper can be applied to the same data sets).
#
# Every replicate is appended to disk as it finishes and the run resumes where it stopped. The 14 tests are computed
# by the same eval_cheap() and the same frozen package source (version 2.4.1, commit 27d5ef2) as the archive.
# Usage: Rscript run_extra.R [smoke]
args <- commandArgs(trailingOnly = TRUE)
SMOKE <- length(args) && args[1] == "smoke"
HERE <- "."
source("simulation/config.R")
JASA$pkg_tree <- "package_2.4.1_27d5ef2"
source(file.path(JASA$study, "R", "harness_core.R"))
OUT <- file.path(HERE, "pvalues", if (SMOKE) "extra_smoke" else "extra")
dir.create(OUT, showWarnings = FALSE)

TESTS <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel",
           "Stukel", "Osius-Rojek", "le-Cessie", "A.poly", "A.spline")
CLASSIC <- c("HL", "Osius-Rojek", "le-Cessie", "Stukel", "Tsiatis", "Xie")
seed_of <- function(tag) strtoi(substr(digest::digest(tag), 1, 7), 16L)

# ---- data generators -------------------------------------------------------------------------------------------------
gen_null_design <- function(design, n) {
  if (design == "rare") {
    x1 <- runif(n, -3, 3); x2 <- rnorm(n)
    d <- data.frame(x1 = x1, x2 = x2); eta <- -2.6 + 0.8 * x1 - 0.5 * x2
  } else {                                   # "six": six covariates, x3 binary
    d <- data.frame(x1 = runif(n, -3, 3), x2 = rnorm(n), x3 = rbinom(n, 1, 0.4), x4 = rnorm(n),
                    x5 = runif(n, -2, 2), x6 = rnorm(n))
    eta <- with(d, 0.2 + 0.6 * x1 - 0.4 * x2 + 0.7 * x3 + 0.3 * x4 - 0.2 * x5 + 0.25 * x6)
  }
  d$y <- rbinom(n, 1, plogis(eta))
  d
}
fit_all <- function(d) suppressWarnings(glm(y ~ ., data = d, family = binomial()))

# ---- combination rules (larger = more misfit); used inside the bootstrap ---------------------------------------------
rule_stats <- function(p) {
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  c(cauchy = mean(tan((0.5 - p) * pi)), hmp = mean(1 / p), fisher = -2 * sum(log(p)), minp = -min(p))
}
panel_stats <- function(pv) c(setNames(rule_stats(pv[TESTS]), paste0("full_", names(rule_stats(pv[TESTS])))),
                              setNames(rule_stats(pv[CLASSIC]), paste0("classic_", names(rule_stats(pv[CLASSIC])))))

# ---- one replicate of each study ---------------------------------------------------------------------------------------
rep_C <- function(i, n, k, mech) {
  s <- seed_of(paste("corrupt", n, sep = "|")) + i            # the clean data of run_corrupt.R
  set.seed(s)
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
  set.seed(s + 104729L)
  idx <- sample.int(n, k)
  if (mech == "flipy") y[idx] <- 1L - y[idx] else x2[idx] <- 4 * x2[idx]
  fit <- suppressWarnings(glm(y ~ x1 + x2, data = data.frame(y = y, x1 = x1, x2 = x2), family = binomial()))
  c(rep = i, seed = s, k = k, eval_cheap(fit))
}
rep_A <- function(i, design, role, n) {
  s <- seed_of(paste("extraA", design, role, n, sep = "|")) + i
  set.seed(s)
  d <- gen_null_design(design, n)
  c(rep = i, seed = s, events = mean(d$y), eval_cheap(fit_all(d)))
}
rep_B <- function(i, scn, n, B) {
  s <- seed_of(paste("extraB", scn, n, sep = "|")) + i
  set.seed(s)
  g <- gen_scenario(scn, n)
  pv <- eval_cheap(g$fit)
  t0 <- panel_stats(pv)
  ph <- as.numeric(fitted(g$fit)); dat <- g$data
  tb <- matrix(NA_real_, B, length(t0), dimnames = list(NULL, names(t0)))
  for (b in seq_len(B)) {
    set.seed(seed_of(paste("extraB", scn, n, i, b, sep = "|")))   # reseed every replicate: some tests touch the RNG
    db <- dat; db$y <- rbinom(length(ph), 1L, ph)
    fb <- suppressWarnings(glm(g$formula, data = db, family = binomial()))
    tb[b, ] <- panel_stats(eval_cheap(fb))
  }
  ok <- stats::complete.cases(tb)
  pboot <- sapply(names(t0), function(r) (1 + sum(tb[ok, r] >= t0[[r]])) / (sum(ok) + 1))
  c(rep = i, seed = s, B_used = sum(ok), pv[TESTS], setNames(pboot, paste0("pboot_", names(t0))))
}

# ---- runner -------------------------------------------------------------------------------------------------------------
cl <- parallel::makeCluster(if (SMOKE) 2L else max(1L, parallel::detectCores() - 1L))
parallel::clusterExport(cl, c("JASA", "TESTS", "CLASSIC", "seed_of", "gen_null_design", "fit_all", "rule_stats",
                              "panel_stats", "rep_C", "rep_A", "rep_B"))
invisible(parallel::clusterEvalQ(cl, { source(file.path(JASA$study, "R", "harness_core.R")); jasa_init(); NULL }))
jasa_init()

run_cell <- function(file, reps, fun, chunk, ...) {
  f <- file.path(OUT, file)
  have <- if (file.exists(f)) utils::read.csv(f)$rep else integer(0)
  todo <- setdiff(seq_len(reps), have)
  message(sprintf("%s: %d done, %d to go", file, length(have), length(todo)))
  for (ch in split(todo, ceiling(seq_along(todo) / chunk))) {
    t0 <- Sys.time()
    rows <- parallel::parLapplyLB(cl, ch, fun, ...)
    df <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
    utils::write.table(df, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
    message(sprintf("  %s: +%d (up to rep %d) %.0f s", file, nrow(df), max(ch), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

R_C <- if (SMOKE) 3L else 1000L; R_A <- if (SMOKE) 3L else 2000L
R_BS <- if (SMOKE) 2L else 500L; R_BP <- if (SMOKE) 2L else 200L
B_S <- if (SMOKE) 5L else 199L; B_P <- if (SMOKE) 5L else 99L
CH <- if (SMOKE) 3L else 92L; CHB <- if (SMOKE) 2L else 23L

for (n in c(500L)) for (mech in c("flipy", "x2")) for (k in c(1L, 5L, 10L))
  run_cell(sprintf("corrupt_%s_k%d_n%d.csv", mech, k, n), R_C, rep_C, CH, n = n, k = k, mech = mech)
for (design in c("rare", "six")) for (n in c(200L, 500L)) for (role in c("reference", "test"))
  run_cell(sprintf("design_%s_%s_n%d.csv", design, role, n), R_A, rep_A, CH, design = design, role = role, n = n)
for (scn in c("quad", "cloglog", "interaction"))
  run_cell(sprintf("boot_%s_n500.csv", scn), R_BP, rep_B, CHB, scn = scn, n = 500L, B = B_P)
for (n in c(500L, 200L))
  run_cell(sprintf("boot_null_n%d.csv", n), R_BS, rep_B, CHB, scn = "null", n = n, B = B_S)
parallel::stopCluster(cl)
message("done")
