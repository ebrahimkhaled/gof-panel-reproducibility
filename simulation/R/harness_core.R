# =============================================================================
# harness_core.R -- init, seeding, per-rep evaluation, cell runner, cluster.
#
# Seeding contract (worker-count invariant):
#   cell_seed(role, family, n) = strtoi(first 7 hex of digest("role|family|n"))
#   per-rep seed               = cell_seed + rep
# Every rep calls set.seed(per-rep seed), generates its dataset, then evaluates
# (the cheap battery is deterministic given the data). Results are therefore
# IDENTICAL for any worker count and any scheduling order. The seeds manifest
# (results/jasa_seeds_manifest.csv) records every cell seed + the pinned
# package commit + sessionInfo.
# =============================================================================

# ---- init: source the pinned package tree + registries ----------------------
jasa_init <- function() {
  if (!exists("JASA", mode = "environment"))
    stop("source config.R before harness_core.R")
  for (f in JASA$pkg_files) source(file.path(JASA$pkg_tree, "R", f))
  source(file.path(JASA$base, "R", "scenarios.R"))           # original 12 DGPs
  source(file.path(JASA$study, "R", "scenarios_registry.R"))
  source(file.path(JASA$study, "R", "methods_registry.R"))
  invisible(TRUE)
}

check_pinned_commit <- function() {
  h <- tryCatch(system2("git", c("-C", shQuote(JASA$pkg_tree), "rev-parse", "HEAD"),
                        stdout = TRUE, stderr = TRUE)[1], error = function(e) NA_character_)
  if (is.na(h) || !nzchar(h)) { warning("could not read the package git HEAD"); return(NA_character_) }
  if (!identical(h, JASA$pkg_commit))
    warning(sprintf("package tree HEAD %s != pinned %s -- results are NOT comparable; checkout the pin or update config.R",
                    substr(h, 1, 9), substr(JASA$pkg_commit, 1, 9)))
  h
}

# ---- seeding ----------------------------------------------------------------
cell_seed <- function(role, family, n)
  strtoi(substr(digest::digest(paste(role, family, n, sep = "|")), 1, 7), 16L)

# ---- per-rep evaluation (cheap battery) -------------------------------------
rep_eval_cheap <- function(i, role, family, n) {
  s <- cell_seed(role, family, n) + i
  set.seed(s)
  g  <- gen_scenario(family, n)
  pv <- eval_cheap(g$fit)
  c(rep = i, seed = s, pv)
}

# ---- cell runner (one CSV per role x family x n; resumable) -----------------
run_cell_cheap <- function(role, family, n, reps, cl = NULL, force = JASA$FORCE) {
  out_csv <- file.path(JASA$stage1, sprintf("%s_%s_n%d_pvalues.csv", role, family, n))
  if (file.exists(out_csv) && !force) {
    nr <- tryCatch(nrow(utils::read.csv(out_csv, check.names = FALSE)), error = function(e) 0L)
    if (nr >= reps) { message("  skip (done): ", basename(out_csv)); return(invisible(out_csv)) }
  }
  t0 <- Sys.time()
  rows <- if (is.null(cl)) lapply(seq_len(reps), rep_eval_cheap, role = role, family = family, n = n)
          else parallel::parLapplyLB(cl, seq_len(reps), rep_eval_cheap,
                                     role = role, family = family, n = n)
  df <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
  utils::write.csv(df, out_csv, row.names = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  %-42s %4d reps  %7.1f s  (%.2f s/rep wall)",
                  basename(out_csv), reps, el, el / reps))
  invisible(out_csv)
}

read_cell <- function(role, family, n) {
  f <- file.path(JASA$stage1, sprintf("%s_%s_n%d_pvalues.csv", role, family, n))
  if (!file.exists(f)) stop("missing stage-1 cell: ", f)
  utils::read.csv(f, check.names = FALSE)
}

# ---- cluster helpers --------------------------------------------------------
make_workers <- function(nw = JASA$WORKERS) {
  if (nw < 2L) return(NULL)
  cl <- parallel::makePSOCKcluster(nw)
  parallel::clusterExport(cl, "JASA", envir = globalenv())
  parallel::clusterCall(cl, function() {
    source(file.path(JASA$study, "R", "harness_core.R")); jasa_init(); NULL
  })
  cl
}

# ---- seeds manifest ---------------------------------------------------------
write_seeds_manifest <- function(cells, commit_hash) {
  # cells: data.frame(role, family, n, reps)
  cells$cell_seed <- mapply(cell_seed, cells$role, cells$family, cells$n)
  cells$pkg_commit <- commit_hash
  cells$r_version  <- as.character(getRversion())
  cells$date       <- format(Sys.Date())
  f <- file.path(JASA$results, "jasa_seeds_manifest.csv")
  utils::write.csv(cells, f, row.names = FALSE)
  si <- file.path(JASA$results, "jasa_sessionInfo.txt")
  writeLines(c(capture.output(sessionInfo()),
               paste("pinned package commit:", commit_hash)), si)
  invisible(f)
}
