# run_corrupt.R -- Result 4 of the panel paper: false alarms of single tests and of combined decisions when a few recorded
# covariate values are corrupted. The data come from the correct model of the LEGofT study (x1 ~ U(-3,3), x2 ~ N(0,1),
# eta0 = 0.3 + 0.8 x1 - 0.5 x2); in the RECORDED data, x1 is multiplied by 4 in k random records (the corruption model of
# the EDGE paper), and the model is fitted to the recorded data. k = 0, 1, 5, 10; the same clean data sets are reused for
# every k (paired). The 14 tests are those of the LEGofT archive, computed by the same eval_cheap() and the same frozen
# package commit (27d5ef2, checked out separately), so the results are comparable with that archive.
# Every replicate is appended to disk as it finishes; the run resumes where it stopped.
# Usage: Rscript run_corrupt.R [smoke]
args <- commandArgs(trailingOnly = TRUE)
SMOKE <- length(args) && args[1] == "smoke"
HERE <- "./legoft_panel"
source("./jasa_study/config.R")
JASA$pkg_tree <- "package_2.4.1_27d5ef2"
source(file.path(JASA$study, "R", "harness_core.R"))

NS <- c(500L, 1000L)
KS <- c(0L, 1L, 5L, 10L)
REPS <- if (SMOKE) c(`0` = 4L, `1` = 3L, `5` = 3L, `10` = 3L) else c(`0` = 2000L, `1` = 1000L, `5` = 1000L, `10` = 1000L)
OUT <- file.path(HERE, if (SMOKE) "corrupt_smoke" else "corrupt")
dir.create(OUT, showWarnings = FALSE)
CHUNK <- if (SMOKE) 3L else 70L

one_rep <- function(i, n, k) {
  s <- strtoi(substr(digest::digest(paste("corrupt", n, sep = "|")), 1, 7), 16L) + i   # same clean data for every k
  set.seed(s)
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
  x1r <- x1
  if (k > 0) {
    set.seed(s + 7919L)
    idx <- sample.int(n, k)
    x1r[idx] <- 4 * x1[idx]
  }
  dat <- data.frame(y = y, x1 = x1r, x2 = x2)
  fit <- suppressWarnings(glm(y ~ x1 + x2, data = dat, family = binomial()))
  c(rep = i, seed = s, k = k, eval_cheap(fit))
}

cl <- parallel::makeCluster(if (SMOKE) 2L else max(1L, parallel::detectCores() - 1L))
parallel::clusterExport(cl, c("JASA", "one_rep"))
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(JASA$study, "R", "harness_core.R")); jasa_init(); NULL
}))
jasa_init()
for (n in NS) for (k in KS) {
  f <- file.path(OUT, sprintf("corrupt_k%d_n%d.csv", k, n))
  done <- if (file.exists(f)) nrow(utils::read.csv(f)) else 0L
  todo <- setdiff(seq_len(REPS[[as.character(k)]]), if (done) utils::read.csv(f)$rep else integer(0))
  message(sprintf("n=%d k=%d: %d done, %d to go", n, k, done, length(todo)))
  for (ch in split(todo, ceiling(seq_along(todo) / CHUNK))) {
    rows <- parallel::parLapplyLB(cl, ch, one_rep, n = n, k = k)
    df <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
    utils::write.table(df, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
    message(sprintf("  n=%d k=%d: +%d (rep %d)", n, k, nrow(df), max(ch)))
  }
}
parallel::stopCluster(cl)
message("done")
