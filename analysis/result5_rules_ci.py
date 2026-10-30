"""result5_rules_ci.py -- referee re-analysis of the stored p-values (no new simulation).

(a) More combination rules, every one calibrated the same way as Cauchy-MC: the statistic is compared with its 95th
    percentile on the independent 'conformal' null set. Rules: Cauchy (Liu and Xie 2020), minP (= Bonferroni made
    exact by simulation), Fisher (-2 sum log p), Stouffer (sum of normal scores), harmonic mean p-value (Wilson 2019)
    and the quadratic form z' S^-1 z (S = correlation of the normal scores on the reference set). Raw Bonferroni is
    kept for comparison.
(b) Uncertainty. 95% bootstrap intervals (B = 1000) for the size of each rule, for average power, and for paired
    power differences against Cauchy-MC and against the best fixed single test (size-adjusted). Each bootstrap
    replicate resamples the reference set (so the critical values move) and the rows of every test set.
(c) Effective number of tests and the five clusters: bootstrap over the null data sets (B = 500); interval for M_eff
    and the share of replicates that give the same five groups.
Outputs: result5_size.csv, result5_power.csv, result5_diff.csv, result5_meff_ci.csv."""
import glob, os, re
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import norm, spearmanr

FAMILY = lambda fname: re.sub(r"^(test|train|conformal)_|_n[0-9]+_pvalues[.]csv$", "", os.path.basename(fname))
HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(HERE, "..", "pvalues", "archive")
A, B, B_MEFF = 0.05, 1000, 500
rng = np.random.default_rng(20261030)
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel",
         "Osius-Rojek", "le-Cessie", "A.poly", "A.spline"]
PANELS = {"CLASSIC": ["HL", "Osius-Rojek", "le-Cessie", "Stukel", "Tsiatis", "Xie"], "FULL": TESTS,
          "ONE-PER-GROUP": ["HL", "Xie", "DEF.poly3", "Osius-Rojek", "le-Cessie"]}
NS = [200, 500, 1000]


def load(kind, fam, n):
    f = os.path.join(S1, f"{kind}_{fam}_n{n}_pvalues.csv")
    return pd.read_csv(f) if os.path.exists(f) else None


def clip(P):
    return np.clip(P, 1e-15, 1 - 1e-15)


def statistics(P, Sinv):
    """Every rule as 'larger = more misfit'."""
    P = clip(P)
    Z = norm.isf(P)
    return {"Cauchy": np.tan((0.5 - P) * np.pi).mean(axis=1),
            "minP": -P.min(axis=1),
            "Fisher": -2 * np.log(P).sum(axis=1),
            "Stouffer": Z.sum(axis=1),
            "HMP": (1 / P).mean(axis=1),
            "Quadratic": np.einsum("ij,jk,ik->i", Z, Sinv, Z)}


RULES = ["Cauchy", "minP", "Fisher", "Stouffer", "HMP", "Quadratic"]
fams = sorted({FAMILY(f) for f in glob.glob(os.path.join(S1, "test_*_pvalues.csv"))} - {"null"})
size_rows, pow_rows, diff_rows = [], [], []

for n in NS:
    ref, tnull = load("conformal", "null", n), load("test", "null", n)
    for pname, cols in PANELS.items():
        k = len(cols)
        Pr = ref[cols].dropna().to_numpy()
        Sinv = np.linalg.inv(np.corrcoef(norm.isf(clip(Pr)), rowvar=False))
        st_ref = statistics(Pr, Sinv)
        P0 = tnull[cols].dropna().to_numpy()
        st0 = statistics(P0, Sinv)
        # data for every family: statistics of each rule + raw Bonferroni + size-adjustable single-test p-values
        cells = []
        for fam in fams:
            d = load("test", fam, n)
            if d is None:
                continue
            P = d[cols].dropna().to_numpy()
            cells.append((fam, statistics(P, Sinv), P.min(axis=1) <= A / k, P))

        def evaluate(ir, i0, irows):
            """Size of each rule, average power of each rule and of every single test (size-adjusted)."""
            crit = {r: np.quantile(st_ref[r][ir], 1 - A) for r in RULES}
            size = {r: (st0[r][i0] > crit[r]).mean() for r in RULES}
            size["Bonferroni"] = (P0[i0].min(axis=1) <= A / k).mean()
            cs = np.quantile(Pr[ir], A, axis=0)
            pw = {r: [] for r in RULES + ["Bonferroni"]}
            single = []
            for (fam, st, bon, P), ix in zip(cells, irows):
                for r in RULES:
                    pw[r].append((st[r][ix] > crit[r]).mean())
                pw["Bonferroni"].append(bon[ix].mean())
                single.append((P[ix] <= cs).mean(axis=0))
            return size, {r: np.array(v) for r, v in pw.items()}, np.array(single)

        full = (np.arange(len(Pr)), np.arange(len(P0)), [np.arange(len(c[3])) for c in cells])
        size, pw, single = evaluate(*full)
        best = int(np.argmax(single.mean(axis=0)))       # best FIXED single test, chosen on the full data
        boots = []
        for b in range(B):
            ib = (rng.integers(0, len(Pr), len(Pr)), rng.integers(0, len(P0), len(P0)),
                  [rng.integers(0, len(c[3]), len(c[3])) for c in cells])
            s, p, sg = evaluate(*ib)
            boots.append((s, {r: v.mean() for r, v in p.items()}, sg.mean(axis=0)[best]))
        lo_hi = lambda x: (np.quantile(x, 0.025), np.quantile(x, 0.975))
        for r in RULES + ["Bonferroni"]:
            sl, sh = lo_hi([bt[0][r] for bt in boots])
            pl, ph = lo_hi([bt[1][r] for bt in boots])
            size_rows.append({"n": n, "panel": pname, "rule": r, "size": size[r], "size_lo": sl, "size_hi": sh})
            pow_rows.append({"n": n, "panel": pname, "rule": r, "families": len(cells), "power": pw[r].mean(),
                             "power_lo": pl, "power_hi": ph,
                             "wins_vs_cauchy": int((pw[r] > pw["Cauchy"] + 1e-12).sum()) if r != "Cauchy" else None})
            if r != "Cauchy":
                dd = [bt[1][r] - bt[1]["Cauchy"] for bt in boots]
                dl, dh = lo_hi(dd)
                diff_rows.append({"n": n, "panel": pname, "comparison": f"{r} - Cauchy",
                                  "diff": pw[r].mean() - pw["Cauchy"].mean(), "lo": dl, "hi": dh})
        dd = [bt[1]["Cauchy"] - bt[2] for bt in boots]
        dl, dh = lo_hi(dd)
        diff_rows.append({"n": n, "panel": pname, "comparison": f"Cauchy - best single ({cols[best]})",
                          "diff": pw["Cauchy"].mean() - single.mean(axis=0)[best], "lo": dl, "hi": dh})
        print(f"n={n} {pname} done", flush=True)

pd.DataFrame(size_rows).to_csv(os.path.join(HERE, "result5_size.csv"), index=False)
pd.DataFrame(pow_rows).to_csv(os.path.join(HERE, "result5_power.csv"), index=False)
pd.DataFrame(diff_rows).to_csv(os.path.join(HERE, "result5_diff.csv"), index=False)


# ---- (c) M_eff and cluster stability
def meff(R):
    lam = np.abs(np.linalg.eigvalsh(R))
    return float(np.sum((lam >= 1) + (lam - np.floor(lam))))


def groups(R):
    lab = fcluster(linkage(squareform(1 - np.abs(R), checks=False), "average"), t=0.5, criterion="distance")
    return frozenset(frozenset(np.flatnonzero(lab == c)) for c in np.unique(lab))


m_rows = []
ic = [TESTS.index(t) for t in PANELS["CLASSIC"]]
for n in [500, 1000]:
    Z = norm.isf(clip(load("test", "null", n)[TESTS].dropna().to_numpy()))
    R = spearmanr(Z).correlation
    g0 = groups(R)
    mf, mc, same = [], [], 0
    for b in range(B_MEFF):
        Rb = spearmanr(Z[rng.integers(0, len(Z), len(Z))]).correlation
        mf.append(meff(Rb)); mc.append(meff(Rb[np.ix_(ic, ic)]))
        same += groups(Rb) == g0
    m_rows.append({"n": n, "M_eff_full": meff(R), "full_lo": np.quantile(mf, .025), "full_hi": np.quantile(mf, .975),
                   "M_eff_classic": meff(R[np.ix_(ic, ic)]), "classic_lo": np.quantile(mc, .025),
                   "classic_hi": np.quantile(mc, .975), "groups": len(g0), "same_groups_share": same / B_MEFF})
pd.DataFrame(m_rows).to_csv(os.path.join(HERE, "result5_meff_ci.csv"), index=False)

pd.set_option("display.width", 220)
print(pd.DataFrame(size_rows).round(3).to_string(index=False))
print(pd.DataFrame(pow_rows).round(3).to_string(index=False))
print(pd.DataFrame(diff_rows).round(4).to_string(index=False))
print(pd.DataFrame(m_rows).round(3).to_string(index=False))
