# check_reproduce_extra.R -- run from the root of this folder: rebuild stored rows of the extra simulations (pvalues/extra/)
# from their seeds with the functions of run_extra.R, and compare all p-values with the stored ones.
source("simulation/config.R")
JASA$pkg_tree <- "package_2.4.1_27d5ef2"
source(file.path(JASA$study, "R", "harness_core.R")); jasa_init()
ex <- parse("simulation/run_extra.R")
keep <- vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
    as.character(e[[2]]) %in% c("TESTS", "CLASSIC", "seed_of", "gen_null_design", "fit_all", "rule_stats", "panel_stats",
                                "rep_C", "rep_A", "rep_B")
}, logical(1))
for (e in ex[keep]) eval(e, globalenv())
cmp <- function(file, row) {
  st <- utils::read.csv(file.path("pvalues", "extra", file), check.names = FALSE)
  st <- st[st$rep == row[["rep"]], names(row), drop = FALSE]
  cat(sprintf("%-32s rep %d: max |recomputed - stored| = %.3g\n", file, row[["rep"]],
              max(abs(unlist(st) - row), na.rm = TRUE)))
}
cmp("corrupt_flipy_k5_n500.csv", rep_C(7, 500L, 5L, "flipy"))
cmp("corrupt_x2_k10_n500.csv", rep_C(3, 500L, 10L, "x2"))
cmp("design_rare_test_n200.csv", rep_A(11, "rare", "test", 200L))
cmp("design_six_reference_n500.csv", rep_A(2, "six", "reference", 500L))
cmp("boot_null_n200.csv", rep_B(5, "null", 200L, 199L))
# second round (run_extra2.R): bootstrap in the other designs and at n = 200
ex2 <- parse("simulation/run_extra2.R")
keep2 <- vapply(ex2, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                  as.character(e[[2]]) %in% c("panel_stats_on", "rep_D"), logical(1))
for (e in ex2[keep2]) eval(e, globalenv())
cmp("bootD_rare_n200.csv", rep_D(4, "rare", 200L, 199L))
cmp("bootE_quad2_n200.csv", rep_B(9, "quad2", 200L, 199L))
