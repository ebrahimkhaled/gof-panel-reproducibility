"""result7_extra2.py -- analysis of the second round of extra simulations (run_extra2.R, stored in extra/).

D  bootstrap calibration in the two other null designs at n = 200: false-alarm rate of the bootstrap p-value at 5% for
   four rules on two panels (Stukel left out in the rare-event design).
E  bootstrap against the fixed reference calibration on the same alternative data sets, for the Cauchy combination:
   n = 200 (bootE_*, B = 199) and, for completeness, the n = 500 runs of run_extra.R (boot_*, B = 99). Reported as the
   two powers, the counts of data sets rejected by one calibration only, and a 95% interval for the paired difference.
Outputs: result7_boot_designs.csv, result7_boot_paired.csv."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
X = os.path.join(HERE, "..", "pvalues", "extra")
S1 = os.path.join(HERE, "..", "pvalues", "archive")
A = 0.05
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel",
         "Osius-Rojek", "le-Cessie", "A.poly", "A.spline"]
PANELS = {"FULL": TESTS, "CLASSIC": ["HL", "Osius-Rojek", "le-Cessie", "Stukel", "Tsiatis", "Xie"]}
cstat = lambda P: np.tan((0.5 - np.clip(P, 1e-15, 1 - 1e-15)) * np.pi).mean(axis=1)

# ---- D: bootstrap size in the other designs -------------------------------------------------------------------------
rows = []
for design in ["six", "rare"]:
    f = os.path.join(X, f"bootD_{design}_n200.csv")
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f)
    for c in [c for c in d.columns if c.startswith("pboot_")]:
        v = d[c].dropna()
        rate = (v <= A).mean()
        rows.append({"design": design, "n": 200, "rule": c[6:], "reps": len(v), "size": rate,
                     "se": np.sqrt(rate * (1 - rate) / len(v)), "events": d.events.mean(), "median_B": d.B_used.median()})
BD = pd.DataFrame(rows)
BD.to_csv(os.path.join(HERE, "result7_boot_designs.csv"), index=False)

# ---- E: paired comparison of the two calibrations ---------------------------------------------------------------------
rows = []
for n, prefix, B, scns in [(200, "bootE", 199, ["quad", "interaction", "quad2"]),
                           (500, "boot", 99, ["quad", "cloglog", "interaction"])]:
    ref = pd.read_csv(os.path.join(S1, f"conformal_null_n{n}_pvalues.csv"))
    for scn in scns:
        f = os.path.join(X, f"{prefix}_{scn}_n{n}.csv")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        for p, key in [("FULL", "full"), ("CLASSIC", "classic")]:
            cols = PANELS[p]
            P = d[cols].to_numpy(); ok = ~np.isnan(P).any(axis=1)
            crit = np.quantile(cstat(ref[cols].dropna().to_numpy()), 1 - A)
            fixed = cstat(P[ok]) > crit
            boot = d["pboot_" + key + "_cauchy"].to_numpy()[ok] <= A
            diff = boot.astype(float) - fixed.astype(float)
            se = diff.std(ddof=1) / np.sqrt(len(diff))
            rows.append({"n": n, "B": B, "alternative": scn, "panel": p, "reps": int(ok.sum()),
                         "power_bootstrap": boot.mean(), "power_fixed": fixed.mean(),
                         "bootstrap_only": int((boot & ~fixed).sum()), "fixed_only": int((~boot & fixed).sum()),
                         "diff": diff.mean(), "lo": diff.mean() - 1.96 * se, "hi": diff.mean() + 1.96 * se})
BE = pd.DataFrame(rows)
BE.to_csv(os.path.join(HERE, "result7_boot_paired.csv"), index=False)

pd.set_option("display.width", 250)
for name, T in [("D bootstrap size, other designs", BD), ("E bootstrap vs fixed reference", BE)]:
    print("\n" + name)
    print(T.round(3).to_string(index=False) if len(T) else "(no data yet)")
