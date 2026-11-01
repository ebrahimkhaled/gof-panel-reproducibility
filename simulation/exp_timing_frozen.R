## exp_timing_frozen.R -- measure the frozen rule's ACTUAL wall clock (referee build item).
## Two honest routes, timed separately on fresh n=500 datasets:
##   (1) scoring cost given the battery p-values (the "milliseconds" claim);
##   (2) battery evaluation for one dataset (what any user must pay);
##   (3) the full real-data pipeline: battery + B=199 parametric-bootstrap calibration
##       (what the applications actually ran). This is the number to print against
##       BAGofT's 427.7 s -- NOT milliseconds.
## Output: results/jasa_stage2/timing_frozen_n500.csv + printed medians.

source("simulation/config.R")
source(file.path(JASA$study, "R", "methods_registry.R"))
source(file.path(JASA$study, "R", "scenarios_registry.R"))
source(file.path(JASA$study, "R", "harness_core.R"))
jasa_init()

REPS <- 10L; B <- 199L; N <- 500L
W <- c(0.0098,0.0258,0.0099,0.1776,0.1994,0.0104,0.0409,0.0401,0.0346,0.2255,0.2260)
tan_of <- function(p) tan((0.5 - pmin(pmax(p, 1e-12), 1-1e-12)) * pi)

rows <- lapply(seq_len(REPS), function(i) {
  set.seed(77000 + i)
  g <- gen_scenario("null", N)
  t1 <- system.time(pv <- eval_cheap(g$fit))["elapsed"]              # battery, one dataset
  t0 <- system.time(for (k in 1:1000) s <- sum(W * tan_of(pv[1:11])))["elapsed"] / 1000  # scoring
  ph <- as.numeric(fitted(g$fit)); dat <- g$data; frm <- g$formula
  t2 <- system.time({                                                 # full bootstrap route
    for (b in seq_len(B)) {
      db <- dat; db$y <- rbinom(length(ph), 1L, ph)
      fb <- suppressWarnings(glm(frm, data = db, family = binomial()))
      invisible(eval_cheap(fb))
    }
  })["elapsed"]
  c(rep = i, scoring_ms = 1000 * t0, battery_s = t1, boot_pipeline_s = t1 + t2)
})
df <- as.data.frame(do.call(rbind, rows))
utils::write.csv(df, file.path(JASA$results, "jasa_stage2", "timing_frozen_n500.csv"), row.names = FALSE)
cat(sprintf("frozen-rule timing over %d fresh n=%d datasets (serial, B=%d):\n", REPS, N, B))
cat(sprintf("  scoring given battery p-values : median %.3f ms\n", median(df$scoring_ms)))
cat(sprintf("  battery evaluation, one dataset: median %.2f s\n",  median(df$battery_s)))
cat(sprintf("  FULL bootstrap pipeline        : median %.1f s   (vs BAGofT measured 427.7 s)\n",
    median(df$boot_pipeline_s)))
