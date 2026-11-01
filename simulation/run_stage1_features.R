# =============================================================================
# run_stage1_features.R -- the cheap-method compute campaign (studies S2 + S4
# data layer). For every (role, family, n) cell it writes one per-rep p-value
# CSV covering ALL cheap methods (battery + Stukel + Osius-Rojek + corrected
# le-Cessie + CCT + A.poly + A.spline). Stage 2 consumes these to build the
# meta-learners; the classical summary tables come straight from these files.
#
# Roles: test      -> per-family evaluation reps (size / power)
#        train     -> meta-training features (M per family, pool P + null)
#        conformal -> out-of-sample null reference (Mc, null only)
#
# Resumable: finished cells are skipped (JASA_FORCE=1 overwrites).
# Run:   Rscript run_stage1_features.R          (full grid)
#        JASA_SMOKE=1 Rscript run_stage1_features.R   (tiny)
# =============================================================================

source("simulation/config.R")
source(file.path(JASA$study, "R", "harness_core.R"))
jasa_init()
commit <- check_pinned_commit()

# ---- build the cell list (Evidence Plan sections 2-3) -----------------------
test_reps <- function(family) {
  if (family == "null") JASA$R_SIZE
  else if (family %in% CONFIRM_ALT) JASA$R_POWER
  else JASA$R_RICHER                       # new5 + deployment D + t8link
}
cells <- do.call(rbind, c(
  # test cells: every family at every n (D + t8link only at n >= 500 in full runs)
  lapply(JASA$NS, function(nn) {
    fams <- if (JASA$SMOKE) ALL_FAMILIES else {
      if (nn < 500) c("null", POOL_P) else ALL_FAMILIES }
    data.frame(role = "test", family = fams, n = nn,
               reps = vapply(fams, test_reps, numeric(1)))
  }),
  # train cells: null + pool P at every n
  lapply(JASA$NS, function(nn)
    data.frame(role = "train", family = c("null", POOL_P), n = nn, reps = JASA$M_TRAIN)),
  # conformal null reference at every n
  lapply(JASA$NS, function(nn)
    data.frame(role = "conformal", family = "null", n = nn, reps = JASA$MC))
))
rownames(cells) <- NULL
cat(sprintf("Stage 1: %d cells, %d total reps, %d workers, smoke=%s\n",
            nrow(cells), sum(cells$reps), JASA$WORKERS, JASA$SMOKE))

write_seeds_manifest(cells, commit)

cl <- make_workers()
t0 <- Sys.time()
for (k in seq_len(nrow(cells)))
  run_cell_cheap(cells$role[k], cells$family[k], cells$n[k], cells$reps[k], cl = cl)
if (!is.null(cl)) parallel::stopCluster(cl)
cat(sprintf("Stage 1 DONE in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
