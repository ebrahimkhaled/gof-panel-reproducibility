"""result6_extra.py -- analysis of the extra simulations of run_extra.R (stored in extra/).

C  other corruption mechanisms (n = 500): flipped responses, or x2 multiplied by 4, in k = 1, 5, 10 records. Single tests
   at their nominal 5%; combined decisions = Cauchy combination with the critical value of the archive's null
   reference set (as in Result 4). k = 0 is replicates 1-1000 of corrupt/corrupt_k0_n500.csv (the same clean data).
A  two more null designs (rare events; six covariates with one binary). False-alarm rate of the four decisions, with
   the Cauchy critical value taken (i) from the design's own reference set and (ii) from the reference set of the
   original design (what happens if one calibrates on the wrong design). Effective number of tests.
B  bootstrap calibration from the fitted model: rejection rate of the bootstrap p-value at 5% for four rules on two
   panels; under the null (size) and under three alternatives (power), where the same data sets are also judged with
   the fixed reference calibration of the paper (paired comparison).
Outputs: result6_corrupt.csv, result6_designs.csv, result6_boot_size.csv, result6_boot_power.csv."""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
X = os.path.join(HERE, "..", "pvalues", "extra")
S1 = os.path.join(HERE, "..", "pvalues", "archive")
A = 0.05
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel",
         "Osius-Rojek", "le-Cessie", "A.poly", "A.spline"]
PANELS = {"FULL": TESTS, "CLASSIC": ["HL", "Osius-Rojek", "le-Cessie", "Stukel", "Tsiatis", "Xie"],
          "ONE-PER-GROUP": ["HL", "Xie", "DEF.poly3", "Osius-Rojek", "le-Cessie"],
          "PROTECTED": ["HL", "EF", "Pigeon-Heyse", "le-Cessie"]}
rd = lambda f: pd.read_csv(os.path.join(X, f)) if os.path.exists(os.path.join(X, f)) else None
cstat = lambda P: np.tan((0.5 - np.clip(P, 1e-15, 1 - 1e-15)) * np.pi).mean(axis=1)
cp = lambda P: 0.5 - np.arctan(cstat(P)) / np.pi


def meff(R):
    lam = np.abs(np.linalg.eigvalsh(R))
    return float(np.sum((lam >= 1) + (lam - np.floor(lam))))


def crit_from(ref, cols):
    return np.quantile(cstat(ref[cols].dropna().to_numpy()), 1 - A)


# ---- C: other corruption mechanisms ----------------------------------------------------------------------------------
ref500 = pd.read_csv(os.path.join(S1, "conformal_null_n500_pvalues.csv"))
k0 = pd.read_csv(os.path.join(HERE, "..", "pvalues", "corruption", "corrupt_k0_n500.csv"))
k0 = k0[k0.rep <= 1000]
rows = []
for mech in ["flipy", "x2"]:
    for k in [0, 1, 5, 10]:
        d = k0 if k == 0 else rd(f"corrupt_{mech}_k{k}_n500.csv")
        if d is None:
            continue
        r = {"mechanism": mech, "k": k, "reps": len(d)}
        for t in TESTS:
            r[t] = (d[t] <= A).mean()
        for p, cols in PANELS.items():
            P = d[cols].dropna().to_numpy()
            r["panel:" + p] = (cstat(P) > crit_from(ref500, cols)).mean()
        rows.append(r)
C = pd.DataFrame(rows)
C.to_csv(os.path.join(HERE, "result6_corrupt.csv"), index=False)

# ---- A: other null designs -----------------------------------------------------------------------------------------------
rows = []
for design in ["rare", "six"]:
    for n in [200, 500]:
        ref, tst = rd(f"design_{design}_reference_n{n}.csv"), rd(f"design_{design}_test_n{n}.csv")
        if ref is None or tst is None or len(tst) < 50:
            continue
        orig_ref = pd.read_csv(os.path.join(S1, f"conformal_null_n{n}_pvalues.csv"))
        na = {t: float(tst[t].isna().mean()) for t in TESTS}
        for p in ["CLASSIC", "FULL"]:
            cols = [t for t in PANELS[p] if na[t] <= 0.02]           # a test that fails in > 2% of data sets is left out
            P = tst[cols].dropna().to_numpy()
            Z = norm.isf(np.clip(P, 1e-15, 1 - 1e-15))
            Rm = pd.DataFrame(Z).corr(method="spearman").to_numpy()
            single = (P <= A).mean(axis=0)
            rows.append({"design": design, "n": n, "panel": p, "tests": len(cols), "left_out": ",".join(sorted(set(PANELS[p]) - set(cols))),
                         "reps": len(P), "events": tst.events.mean(), "ANY": (P.min(axis=1) <= A).mean(),
                         "Bonferroni": (P.min(axis=1) <= A / len(cols)).mean(), "Cauchy": (cp(P) <= A).mean(),
                         "Cauchy-MC": (cstat(P) > crit_from(ref, cols)).mean(),
                         "Cauchy-MC_wrong_design": (cstat(P) > crit_from(orig_ref, cols)).mean(),
                         "single_min": single.min(), "single_max": single.max(), "M_eff": meff(Rm),
                         "max_na": max(na[t] for t in PANELS[p])})
Dz = pd.DataFrame(rows)
Dz.to_csv(os.path.join(HERE, "result6_designs.csv"), index=False)

# ---- B: bootstrap calibration ---------------------------------------------------------------------------------------------
RULES = [c[6:] for c in ["pboot_full_cauchy", "pboot_full_hmp", "pboot_full_fisher", "pboot_full_minp",
                         "pboot_classic_cauchy", "pboot_classic_hmp", "pboot_classic_fisher", "pboot_classic_minp"]]
rows = []
for n in [200, 500]:
    d = rd(f"boot_null_n{n}.csv")
    if d is None:
        continue
    for r in RULES:
        v = d["pboot_" + r].dropna()
        rate = (v <= A).mean()
        rows.append({"n": n, "rule": r, "reps": len(v), "size": rate, "se": np.sqrt(rate * (1 - rate) / max(len(v), 1)),
                     "median_B": d.B_used.median()})
BS = pd.DataFrame(rows)
BS.to_csv(os.path.join(HERE, "result6_boot_size.csv"), index=False)
rows = []
for scn in ["quad", "cloglog", "interaction"]:
    d = rd(f"boot_{scn}_n500.csv")
    if d is None:
        continue
    for p, key in [("FULL", "full"), ("CLASSIC", "classic")]:
        cols = PANELS[p]
        P = d[cols].to_numpy(); ok = ~np.isnan(P).any(axis=1)
        fixed = cstat(P[ok]) > crit_from(ref500, cols)
        boot = (d["pboot_" + key + "_cauchy"].to_numpy()[ok] <= A)
        rows.append({"alternative": scn, "panel": p, "reps": int(ok.sum()), "power_bootstrap": boot.mean(),
                     "power_fixed_reference": fixed.mean(), "agree": (boot == fixed).mean()})
BP = pd.DataFrame(rows)
BP.to_csv(os.path.join(HERE, "result6_boot_power.csv"), index=False)

pd.set_option("display.width", 250)
for name, T in [("C corruption", C), ("A designs", Dz), ("B bootstrap size", BS), ("B bootstrap power", BP)]:
    print("\n" + name)
    print(T.round(3).to_string(index=False) if len(T) else "(no data yet)")
