# check_reproduce_row.R -- the frozen commit (27d5ef2, checked out separately) must reproduce stored rows of the LEGofT
# archive exactly: regenerate rows from their stored seeds and compare all 15 p-values.
source("./jasa_study/config.R")
JASA$pkg_tree <- "package_2.4.1_27d5ef2"
source(file.path(JASA$study, "R", "harness_core.R")); jasa_init()
st <- utils::read.csv(file.path(JASA$stage1, "test_null_n500_pvalues.csv"), check.names = FALSE)
for (i in 1:3) {
  set.seed(st$seed[i]); g <- gen_scenario("null", 500); pv <- eval_cheap(g$fit)
  d <- max(abs(pv - unlist(st[i, names(pv)])), na.rm = TRUE)
  cat(sprintf("row %d seed %d: max |recomputed - stored| = %.3g\n", i, st$seed[i], d))
}
