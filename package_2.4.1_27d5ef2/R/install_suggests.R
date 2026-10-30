# Optional ("Suggests") packages used by the run.all.gof() battery, and a
# user-consented installer for them. CRAN forbids a package from installing
# anything silently, so nothing here ever installs without either an interactive
# yes/no confirmation (the default) or the caller explicitly setting ask = FALSE,
# and run.all.gof() only ever reaches the installer in an interactive session.

# Which optional package each test needs to run (only tests that are skipped
# entirely when the package is absent; Stukel still runs without statmod).
.GOF_TEST_PKGS <- list(
  "GiViTI"          = c("givitiR", "callr"),
  "GiViTI-external" = c("givitiR", "callr"),
  "HL-GAM"          = "mgcv",
  "PR-GAM"          = "mgcv",
  "Xie-GAM"         = "mgcv",
  "BAGofT"          = "BAGofT",
  "Lai-Liu-HL"      = "ResourceSelection"
)

# The full optional set, including the ones that only improve (not gate) a test.
.GOF_ALL_SUGGESTS <- c("givitiR", "callr", "mgcv", "BAGofT",
                       "ResourceSelection", "statmod", "CompQuadForm")

#' Install the optional packages used by \code{run.all.gof()}
#'
#' The slow tests in \code{\link{run.all.gof}} rely on optional packages that
#' live in \code{Suggests} (\pkg{givitiR} and \pkg{callr} for the GiViTI
#' calibration test, \pkg{mgcv} for the GAM tests, \pkg{BAGofT} for the adaptive
#' test, and \pkg{ResourceSelection} for the Lai-Liu test). Per CRAN policy the
#' package never installs them on its own; this helper installs the missing ones
#' for you, asking first.
#'
#' @param pkgs Optional character vector of package names to consider. Defaults
#'   to the full optional set used by the battery.
#' @param ask Logical; when \code{TRUE} (the default in interactive sessions)
#'   you are shown the list of packages to be installed/updated and asked to
#'   confirm first. Set \code{ask = FALSE} to proceed without a prompt (e.g. in a
#'   setup script you control).
#' @param update Logical; when \code{FALSE} (the default) only the \emph{missing}
#'   packages are installed and anything already present is left untouched. When
#'   \code{TRUE} the function also checks (via \code{\link[utils]{old.packages}})
#'   which of the present packages are out of date and offers to update those
#'   too. The update check contacts your CRAN mirror, so it is a little slower.
#' @return Invisibly, the character vector of packages that were installed or
#'   updated (empty if nothing was needed).
#' @examples
#' \dontrun{
#' # install whatever optional packages are missing, after confirming:
#' gof_install_suggests()
#'
#' # also update any that are out of date:
#' gof_install_suggests(update = TRUE)
#'
#' # just the GiViTI dependencies, no prompt:
#' gof_install_suggests(c("givitiR", "callr"), ask = FALSE)
#' }
#' @seealso \code{\link{run.all.gof}}
#' @export
gof_install_suggests <- function(pkgs = NULL, ask = interactive(),
                                 update = FALSE) {
  want    <- if (is.null(pkgs)) .GOF_ALL_SUGGESTS else pkgs
  present <- vapply(want, requireNamespace, logical(1), quietly = TRUE)
  missing <- want[!present]

  # installed-but-out-of-date ones, only when asked (old.packages() hits CRAN)
  outdated <- character(0)
  if (isTRUE(update)) {
    old <- tryCatch(rownames(utils::old.packages()), error = function(e) NULL)
    outdated <- intersect(want[present], old)
  }

  todo <- unique(c(missing, outdated))
  if (!length(todo)) {
    message("All optional packages for run.all.gof() are already installed",
            if (isTRUE(update)) " and up to date." else ".")
    return(invisible(character(0)))
  }

  if (isTRUE(ask)) {
    if (!interactive()) {
      message("Optional package(s) to install/update: ",
              paste(todo, collapse = ", "),
              ". Run gof_install_suggests() interactively, or with ask = FALSE.")
      return(invisible(todo))
    }
    parts <- character(0)
    if (length(missing))
      parts <- c(parts, sprintf("install %s", paste(missing, collapse = ", ")))
    if (length(outdated))
      parts <- c(parts, sprintf("update %s", paste(outdated, collapse = ", ")))
    ans <- readline(sprintf("For run.all.gof(), %s? [y/N] ",
                            paste(parts, collapse = "; ")))
    if (!tolower(trimws(ans)) %in% c("y", "yes")) {
      message("Skipped. You can install them later with gof_install_suggests().")
      return(invisible(todo))
    }
  }

  utils::install.packages(todo)
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still))
    warning("Could not install: ", paste(still, collapse = ", "))
  invisible(todo)
}

# Called by run.all.gof(): when interactive and tests in the run need optional
# packages that aren't installed, offer to install them (consent-gated). Never
# does anything in a non-interactive session, so R CMD check / scripts / CRAN are
# untouched and the library is never modified without the user's say-so.
.gof_maybe_install <- function(sel, include_slow, install = "ask") {
  install <- match.arg(install, c("ask", "no", "yes"))
  if (identical(install, "no") || !interactive()) return(invisible())
  eff <- sel
  if (!isTRUE(include_slow))
    eff <- eff[!vapply(eff, function(nm) isTRUE(.GOF_REGISTRY[[nm]]$slow), logical(1))]
  needed  <- unique(unlist(.GOF_TEST_PKGS[intersect(eff, names(.GOF_TEST_PKGS))],
                           use.names = FALSE))
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible())
  gof_install_suggests(missing, ask = identical(install, "ask"))
  invisible()
}
