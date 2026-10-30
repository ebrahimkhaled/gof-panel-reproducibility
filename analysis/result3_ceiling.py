"""result3_ceiling.py -- Result 3 of the panel paper: how much can ANY weighting of a fixed battery gain over equal weights?
From the LEGofT archive (no simulation). Each test's p-value becomes a normal score z = qnorm(1 - p).

(a) Theoretical ceiling (local, normal approximation). With null covariance S of the scores and mean shift d under an
    alternative, the best linear combination has noncentrality d' S^-1 d; equal weights give (1'd)^2 / (1'S1). Their ratio
    is the efficiency of equal weights; power at level 0.05 follows from Phi(sqrt(NC) - 1.645).
(b) Empirical ceiling with an honest split. The ORACLE knows the alternative: its weights w = S^-1 d are estimated on that
    family's TRAINING replicates (200) and its power measured on the separate TEST replicates (1,000), with the critical
    value from the independent conformal null set. Compared with equal-weight linear (Stouffer) and Cauchy-MC.
(c) For reference, the stored held-out results of the learned LEGofT rule and equal-weight Cauchy over the same 11 inputs
    (results/jasa_stage2/FINAL_heldout_n*.csv).
Outputs: result3_ceiling.csv."""
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import norm
FAMILY = lambda fname: re.sub(r"^(test|train|conformal)_|_n[0-9]+_pvalues[.]csv$", "", os.path.basename(fname))   # family names can contain underscores

HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(HERE, "..", "pvalues", "archive")
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel",
         "Osius-Rojek", "le-Cessie", "A.poly", "A.spline"]
A, ZA = 0.05, norm.ppf(0.95)


def Z(kind, fam, n):
    f = os.path.join(S1, f"{kind}_{fam}_n{n}_pvalues.csv")
    if not os.path.exists(f):
        return None
    d = pd.read_csv(f)[TESTS].dropna()
    return norm.ppf(1 - d.clip(1e-12, 1 - 1e-12).to_numpy())


def cauchy_stat(Zm):
    P = 1 - norm.cdf(Zm)
    return np.tan((0.5 - np.clip(P, 1e-12, 1 - 1e-12)) * np.pi).mean(axis=1)


rows = []
for n in [500, 1000]:
    ref, tnull, trnull = Z("conformal", "null", n), Z("test", "null", n), Z("train", "null", n)
    mu0, S = ref.mean(axis=0), np.cov(ref, rowvar=False)
    Sinv = np.linalg.pinv(S)
    one = np.ones(len(TESTS))
    crit_eq = np.quantile(ref @ one, 1 - A)
    crit_cct = np.quantile(cauchy_stat(ref), 1 - A)
    fams = sorted({FAMILY(f) for f in os.listdir(S1) if f.startswith("train_") and f.endswith(f"_n{n}_pvalues.csv")} - {"null"})
    for fam in fams:
        tr, te = Z("train", fam, n), Z("test", fam, n)
        if tr is None or te is None:
            continue
        # (a) theory, from the test replicates (population-level shift)
        d = te.mean(axis=0) - mu0
        nc_opt = float(d @ Sinv @ d)
        nc_eq = float((one @ d) ** 2 / (one @ S @ one))
        # (b) oracle weights from the TRAINING replicates of this family, evaluated on TEST
        w = Sinv @ (tr.mean(axis=0) - trnull.mean(axis=0))
        crit_w = np.quantile(ref @ w, 1 - A)
        rows.append({"n": n, "family": fam,
                     "eff_equal": nc_eq / nc_opt if nc_opt > 0 else np.nan,
                     "theory_power_opt": norm.cdf(np.sqrt(nc_opt) - ZA), "theory_power_eq": norm.cdf(np.sqrt(nc_eq) - ZA),
                     "oracle_power": float((te @ w >= crit_w).mean()),
                     "equal_linear_power": float((te @ one >= crit_eq).mean()),
                     "cauchy_mc_power": float((cauchy_stat(te) >= crit_cct).mean())})
out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, "result3_ceiling.csv"), index=False)
pd.set_option("display.width", 220)
print(out.round(3).to_string(index=False))
print()
print(out.groupby("n")[["eff_equal", "theory_power_opt", "theory_power_eq", "oracle_power", "equal_linear_power", "cauchy_mc_power"]].mean().round(3).to_string())
print()
for n in [500, 1000]:
    f = os.path.join(HERE, "..", "heldout_and_timing", f"FINAL_heldout_n{n}.csv")
    if os.path.exists(f):
        h = pd.read_csv(f)
        print(f"stored held-out results n={n}: columns {list(h.columns)[:12]}")
        print(h.head(14).round(3).to_string(index=False))
