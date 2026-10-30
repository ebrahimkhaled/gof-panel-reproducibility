# =============================================================================
# config.R -- JASA study harness configuration (single source of truth).
# Every runner sources this first. All knobs are env-var overridable so the
# same scripts serve smoke tests (JASA_SMOKE=1) and the full campaign.
# =============================================================================

JASA <- new.env()

JASA$base    <- "."
JASA$study   <- file.path(JASA$base, "jasa_study")
JASA$results <- file.path(JASA$base, "results")
JASA$stage1  <- file.path(JASA$results, "jasa_stage1")   # per-rep cheap-method p-values + seeds
JASA$stage2  <- file.path(JASA$results, "jasa_stage2")   # per-rep meta p-values + summaries
JASA$bagoft  <- file.path(JASA$results, "jasa_bagoft")   # per-rep BAGofT p-values (append-mode)
for (d in c(JASA$stage1, JASA$stage2, JASA$bagoft)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- pinned package (v2.4.1 main; le-Cessie transpose fix). Decision D2. ------
# The harness SOURCES the git working tree (NOT the installed CRAN 2.4.0, whose
# le-Cessie is the buggy pre-2.4.1 version). check_pinned_commit() warns if the
# tree has moved off the pinned commit.
JASA$pkg_tree   <- "package_2.4.1_27d5ef2"
JASA$pkg_commit <- "27d5ef20613dc7d5cf0dade807988d8b77b0fed9"
JASA$pkg_files  <- c("ebrahim_farrington_test.R", "def_gof.R", "def_ensemble_gof.R",
                     "install_suggests.R", "run_all_gof.R", "cdef_gof.R", "edge_gof.R")

# --- reps / grids (env-overridable; JASA_SMOKE=1 flips all to tiny) -----------
.smoke <- identical(Sys.getenv("JASA_SMOKE", "0"), "1")
.genv  <- function(nm, def) as.integer(Sys.getenv(nm, as.character(def)))
JASA$SMOKE    <- .smoke
JASA$NS       <- as.integer(strsplit(Sys.getenv("JASA_NS", if (.smoke) "300" else "200,500,1000"), ",")[[1]])
JASA$R_SIZE   <- .genv("JASA_R_SIZE",   if (.smoke) 6L else 2000L)  # null (size) reps per n
JASA$R_POWER  <- .genv("JASA_R_POWER",  if (.smoke) 6L else 1000L)  # confirmatory-alt reps
JASA$R_RICHER <- .genv("JASA_R_RICHER", if (.smoke) 6L else 500L)   # new5 / deployment / t8link reps
JASA$M_TRAIN  <- .genv("JASA_M_TRAIN",  if (.smoke) 8L  else 200L)  # training features per family
JASA$MC       <- .genv("JASA_MC",       if (.smoke) 40L else 2000L) # conformal null reference (D3)
JASA$WORKERS  <- .genv("JASA_WORKERS",  max(1L, parallel::detectCores() - 1L))
JASA$FORCE    <- identical(Sys.getenv("JASA_FORCE", "0"), "1")      # 1 = overwrite finished cells

# --- BAGofT block (cost-managed; own flag -- Evidence Plan section 4) ---------
JASA$RUN_BAGOFT     <- identical(Sys.getenv("JASA_RUN_BAGOFT", "0"), "1")
JASA$BAG_NSPLITS    <- .genv("JASA_BAG_NSPLITS", 100L)   # published default
JASA$BAG_NSIM       <- .genv("JASA_BAG_NSIM",    100L)   # published default
JASA$BAG_NTREE      <- .genv("JASA_BAG_NTREE",    60L)   # parRF() published default
JASA$BAG_REPS_V500  <- .genv("JASA_BAG_REPS_V500",  400L)
JASA$BAG_REPS_V1000 <- .genv("JASA_BAG_REPS_V1000", 300L)
JASA$BAG_REPS_H500  <- .genv("JASA_BAG_REPS_H500",  300L)
JASA$BAG_REPS_H1000 <- .genv("JASA_BAG_REPS_H1000", 150L)
JASA$BAG_REPS_K     <- .genv("JASA_BAG_REPS_K",      40L)

# --- seeding ------------------------------------------------------------------
# Per-cell seed = first 7 hex digits of digest(role|family|n); per-rep seed =
# cell seed + rep. Every rep is fully self-seeded, so results are invariant to
# the worker count and to load-balancing order. The BAGofT block regenerates
# the IDENTICAL test datasets by reusing the "test" role's cell seed; BAGofT's
# own internal RNG is then re-seeded with a fixed offset.
JASA$METHOD_RNG_OFFSET <- 500000000L
