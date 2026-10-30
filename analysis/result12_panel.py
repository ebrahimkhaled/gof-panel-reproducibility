"""result12_panel.py -- Results 1 and 2 of the panel paper, from the LEGofT archive (learned-def/results/jasa_stage1/), where
one call evaluates every test on the same data set (alignment check: the stored Cauchy column equals the Cauchy combination
recomputed from the same row's three basis p-values, to 3e-15). No simulation.

Result 1 -- one decision from a panel. Panels: CLASSIC (HL, Osius-Rojek, le Cessie, Stukel, Tsiatis, Xie: tests an analyst
commonly runs side by side) and FULL (the 14 tests of the archive). Decisions at alpha = 0.05: ANY (reject if any member
rejects), Bonferroni, Cauchy (equal-weight Cauchy combination, analytic p-value), Cauchy-MC (the Cauchy statistic compared
with the 95th percentile of its values on the independent 'conformal' null set: calibration by simulation). Size is
measured on the separate 'test' null set (2,000 replicates per n); power on the 23 'test' alternative families.

Result 2 -- redundancy. Normal scores z = qnorm(1 - p); Spearman correlation across replicates; effective number of
independent tests M_eff (Li and Ji 2005); clusters by average linkage on 1 - |r| (cut at |r| = 0.5).
Outputs: result1_size.csv, result1_power.csv, result2_meff.csv, result2_corr_null_n1000.csv, result2_clusters.txt."""
import glob, os
import re
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import norm
FAMILY = lambda fname: re.sub(r"^(test|train|conformal)_|_n[0-9]+_pvalues[.]csv$", "", os.path.basename(fname))   # family names can contain underscores

HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(HERE, "..", "pvalues", "archive")
A = 0.05
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel",
         "Osius-Rojek", "le-Cessie", "A.poly", "A.spline"]
PANELS = {"CLASSIC": ["HL", "Osius-Rojek", "le-Cessie", "Stukel", "Tsiatis", "Xie"], "FULL": TESTS,
          "ONE-PER-GROUP": ["HL", "Xie", "DEF.poly3", "Osius-Rojek", "le-Cessie"]}   # one test from each cluster of result 2
NS = [200, 500, 1000]   # n = 300 holds only 6 null replicates (smoke runs): excluded


def load(kind, fam, n):
    f = os.path.join(S1, f"{kind}_{fam}_n{n}_pvalues.csv")
    return pd.read_csv(f) if os.path.exists(f) else None


def cauchy_stat(P):
    return np.tan((0.5 - np.clip(P, 1e-15, 1 - 1e-15)) * np.pi).mean(axis=1)


def cauchy_p(P):
    return 0.5 - np.arctan(cauchy_stat(P)) / np.pi


def rates(P, k, crit_mc):
    return {"ANY": (P.min(axis=1) <= A).mean(), "Bonferroni": (P.min(axis=1) <= A / k).mean(),
            "Cauchy": (cauchy_p(P) <= A).mean(), "Cauchy-MC": (cauchy_stat(P) >= crit_mc).mean()}


fams = sorted({FAMILY(f) for f in glob.glob(os.path.join(S1, "test_*_pvalues.csv"))} - {"null"})
size_rows, pow_rows = [], []
for n in NS:
    ref, tnull = load("conformal", "null", n), load("test", "null", n)
    if ref is None or tnull is None:
        continue
    for pname, cols in PANELS.items():
        k = len(cols)
        P_ref = ref[cols].dropna().to_numpy()
        crit_mc = np.quantile(cauchy_stat(P_ref), 1 - A)
        P0 = tnull[cols].dropna().to_numpy()
        r0 = rates(P0, k, crit_mc)
        single = {c: (P0[:, j] <= A).mean() for j, c in enumerate(cols)}
        size_rows.append({"n": n, "panel": pname, "reps": len(P0), **r0,
                          "single_min": min(single.values()), "single_max": max(single.values())})
        # size-adjusted thresholds for the single tests (from the reference null set) to find the best single test fairly
        crit_single = {c: np.quantile(P_ref[:, j], A) for j, c in enumerate(cols)}
        for fam in fams:
            d = load("test", fam, n)
            if d is None:
                continue
            P = d[cols].dropna().to_numpy()
            r = rates(P, k, crit_mc)
            singles = {c: (P[:, j] <= crit_single[c]).mean() for j, c in enumerate(cols)}
            best = max(singles, key=singles.get)
            pow_rows.append({"n": n, "panel": pname, "family": fam, "reps": len(P), **r,
                             "best_single": best, "best_single_power": singles[best], **{"single:" + c: v for c, v in singles.items()}})
size = pd.DataFrame(size_rows)
power = pd.DataFrame(pow_rows)
size.to_csv(os.path.join(HERE, "result1_size.csv"), index=False)
power.to_csv(os.path.join(HERE, "result1_power.csv"), index=False)


# ---- result 2
def meff(R):
    lam = np.abs(np.linalg.eigvalsh(R))
    return float(np.sum((lam >= 1) + (lam - np.floor(lam))))


m_rows, clus = [], []
for n in [500, 1000]:
    for kind, fam in [("test", "null")] + [("test", f) for f in fams]:
        d = load(kind, fam, n)
        if d is None:
            continue
        Z = d[TESTS].dropna().apply(lambda c: norm.ppf(1 - c.clip(1e-15, 1 - 1e-15)))
        R = Z.corr(method="spearman")
        if R.isna().any().any():
            continue
        Rc = R.loc[PANELS["CLASSIC"], PANELS["CLASSIC"]]
        m_rows.append({"n": n, "family": fam, "reps": len(Z), "M_eff_FULL14": meff(R.to_numpy()),
                       "M_eff_CLASSIC6": meff(Rc.to_numpy()),
                       "median_abs_r": float(np.median(np.abs(R.to_numpy()[np.triu_indices(14, 1)])))})
        if fam == "null":
            R.round(3).to_csv(os.path.join(HERE, f"result2_corr_null_n{n}.csv"))
            lab = fcluster(linkage(squareform(1 - np.abs(R.to_numpy()), checks=False), "average"), t=0.5, criterion="distance")
            g = {}
            for name, c in zip(R.columns, lab):
                g.setdefault(c, []).append(name)
            clus.append(f"null n={n}: " + " | ".join(", ".join(v) for v in g.values()))
meffs = pd.DataFrame(m_rows)
meffs.to_csv(os.path.join(HERE, "result2_meff.csv"), index=False)
open(os.path.join(HERE, "result2_clusters.txt"), "w", encoding="utf-8").write("\n".join(clus) + "\n")

pd.set_option("display.width", 220)
print("RESULT 1 -- SIZE on the independent test null set (nominal 0.05)")
print(size.round(3).to_string(index=False))
print("\nRESULT 1 -- POWER, mean over the 23 families (raw for ANY/Bonferroni/Cauchy; Cauchy-MC calibrated)")
print(power.groupby(["panel", "n"])[["ANY", "Bonferroni", "Cauchy", "Cauchy-MC", "best_single_power"]].mean().round(3).to_string())
print("\nbest FIXED single test (same test for every family), size-adjusted, mean over families:")
for (pn, n), g in power.groupby(["panel", "n"]):
    sc = {c[7:]: g[c].mean() for c in g.columns if c.startswith("single:") and g[c].notna().all()}
    b = max(sc, key=sc.get)
    print(f"  {pn:14s} n={n:5d}: {b} {sc[b]:.3f}   Cauchy-MC {g['Cauchy-MC'].mean():.3f}   families where Cauchy-MC >= fixed best: {(g['Cauchy-MC'] >= g['single:' + b] - 1e-12).sum()}/{len(g)}")
print("\nRESULT 2 -- effective number of independent tests")
print(meffs.groupby("n")[["M_eff_FULL14", "M_eff_CLASSIC6", "median_abs_r"]].describe().round(2).T.to_string())
print("\n" + "\n".join(clus))
