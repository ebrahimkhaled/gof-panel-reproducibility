# run_extra2.R -- second round of extra simulations, asked for by the referee pass of 2026-11-01.
#
#   D  bootstrap calibration in the two other null designs of run_extra.R ("rare" and "six"), n = 200:
#      500 data sets each, B = 199. In the rare-event design Stukel's test is left out of both panels, as in the
#      analysis of the fixed-reference results (it fails to compute in 17% of data sets at n = 200).
#   E  bootstrap power at n = 200 on three alternatives of moderate power (quad, interaction, quad2):
#      300 data sets each, B = 199. The observed p-values are stored, so the fixed reference calibration of the
#      paper can be applied to the same data sets.
#
# Same tests, same frozen package source and same seeding scheme as run_extra.R, whose functions are reused.
# Every replicate is appended to disk as it finishes and the run resumes where it stopped.
# Usage: Rscript run_extra2.R [smoke]
args <- commandArgs(trailingOnly = TRUE)
SMOKE <- length(args) && args[1] == "smoke"
HERE <- "."
source("simulation/config.R")
JASA$pkg_tree <- "package_2.4.1_27d5ef2"
source(file.path(JASA$study, "R", "harness_core.R"))
OUT <- file.path(HERE, "pvalues", if (SMOKE) "extra_smoke" else "extra")
dir.create(OUT, showWarnings = FALSE)

# reuse the definitions of run_extra.R (tests, generators, combination rules, rep_B) without running it
defs <- parse(file.path("simulation", "run_extra.R"))
keep <- vapply(defs, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                 as.character(e[[2]]) %in% c("TESTS", "CLASSIC", "seed_of", "gen_null_design", "fit_all",
                                             "rule_stats", "panel_stats", "rep_B"), logical(1))
for (e in defs[keep]) eval(e, globalenv())

# panel statistics on a given pair of test sets (the rare-event design drops Stukel)
panel_stats_on <- function(pv, full, classic) {
  a <- rule_stats(pv[full]); b <- rule_stats(pv[classic])
  c(setNames(a, paste0("full_", names(a))), setNames(b, paste0("classic_", names(b))))
}

rep_D <- function(i, design, n, B) {
  full <- if (design == "rare") setdiff(TESTS, "Stukel") else TESTS
  classic <- if (design == "rare") setdiff(CLASSIC, "Stukel") else CLASSIC
  s <- seed_of(paste("extraD", design, n, sep = "|")) + i
  set.seed(s)
  d <- gen_null_design(design, n)
  fit <- fit_all(d)
  pv <- eval_cheap(fit)
  t0 <- panel_stats_on(pv, full, classic)
  ph <- as.numeric(fitted(fit))
  tb <- matrix(NA_real_, B, length(t0), dimnames = list(NULL, names(t0)))
  for (b in seq_len(B)) {
    set.seed(seed_of(paste("extraD", design, n, i, b, sep = "|")))   # reseed every replicate: some tests touch the RNG
    db <- d; db$y <- rbinom(length(ph), 1L, ph)
    tb[b, ] <- panel_stats_on(eval_cheap(fit_all(db)), full, classic)
  }
  ok <- stats::complete.cases(tb)
  pboot <- if (any(is.na(t0))) rep(NA_real_, length(t0)) else
    sapply(names(t0), function(r) (1 + sum(tb[ok, r] >= t0[[r]])) / (sum(ok) + 1))
  c(rep = i, seed = s, events = mean(d$y), B_used = sum(ok), pv[TESTS], setNames(pboot, paste0("pboot_", names(t0))))
}

cl <- parallel::makeCluster(if (SMOKE) 2L else max(1L, parallel::detectCores() - 1L))
parallel::clusterExport(cl, c("JASA", "TESTS", "CLASSIC", "seed_of", "gen_null_design", "fit_all", "rule_stats",
                              "panel_stats", "panel_stats_on", "rep_B", "rep_D"))
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

R_D <- if (SMOKE) 2L else 500L; R_E <- if (SMOKE) 2L else 300L
B <- if (SMOKE) 5L else 199L; CHB <- if (SMOKE) 2L else 23L
for (design in c("six", "rare"))
  run_cell(sprintf("bootD_%s_n200.csv", design), R_D, rep_D, CHB, design = design, n = 200L, B = B)
for (scn in c("quad", "interaction", "quad2"))
  run_cell(sprintf("bootE_%s_n200.csv", scn), R_E, rep_B, CHB, scn = scn, n = 200L, B = B)
parallel::stopCluster(cl)
message("done")
