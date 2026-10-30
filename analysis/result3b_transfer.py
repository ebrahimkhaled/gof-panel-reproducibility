"""result3b_transfer.py -- why the large reweighting ceiling of result 3 is not reached by learning across alternatives.
(1) Transfer matrix: oracle weights fitted on family A's training replicates, power measured on family B's test replicates.
(2) Pooled weights: weights fitted on the pooled training shift of all OTHER families (leave-one-family-out), as a learner
    that does not know the alternative must do; power on the left-out family's test replicates.
(3) Angles between the oracle weight vectors of different families (in the S-metric): if they point in different
    directions, the pooled direction is close to an average that no single family wants.
LEGofT archive only; no simulation. Outputs: result3b_transfer.csv, result3b_pooled.csv, result3b_angles.csv."""
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
A = 0.05


def Z(kind, fam, n):
    f = os.path.join(S1, f"{kind}_{fam}_n{n}_pvalues.csv")
    if not os.path.exists(f):
        return None
    d = pd.read_csv(f)[TESTS].dropna()
    return norm.ppf(1 - d.clip(1e-12, 1 - 1e-12).to_numpy())


def cauchy_stat(Zm):
    P = 1 - norm.cdf(Zm)
    return np.tan((0.5 - np.clip(P, 1e-12, 1 - 1e-12)) * np.pi).mean(axis=1)


T_rows, P_rows, G_rows = [], [], []
for n in [500, 1000]:
    ref, trnull = Z("conformal", "null", n), Z("train", "null", n)
    S = np.cov(ref, rowvar=False); Sinv = np.linalg.pinv(S)
    fams = sorted({FAMILY(f) for f in os.listdir(S1) if f.startswith("train_") and f.endswith(f"_n{n}_pvalues.csv")} - {"null"})
    tr = {f: Z("train", f, n) for f in fams}; te = {f: Z("test", f, n) for f in fams}
    fams = [f for f in fams if tr[f] is not None and te[f] is not None]
    shift = {f: tr[f].mean(axis=0) - trnull.mean(axis=0) for f in fams}
    W = {f: Sinv @ shift[f] for f in fams}
    power = lambda w, Zt: float((Zt @ w >= np.quantile(ref @ w, 1 - A)).mean())
    crit_c = np.quantile(cauchy_stat(ref), 1 - A)
    for a in fams:
        for b in fams:
            T_rows.append({"n": n, "weights_from": a, "tested_on": b, "power": power(W[a], te[b])})
    for b in fams:
        pooled = Sinv @ np.mean([shift[f] for f in fams if f != b], axis=0)
        P_rows.append({"n": n, "family": b, "oracle": power(W[b], te[b]), "pooled_LOFO": power(pooled, te[b]),
                       "equal_linear": power(np.ones(len(TESTS)), te[b]),
                       "cauchy_mc": float((cauchy_stat(te[b]) >= crit_c).mean())})
    for i, a in enumerate(fams):
        for b in fams[i + 1:]:
            cos = W[a] @ S @ W[b] / np.sqrt((W[a] @ S @ W[a]) * (W[b] @ S @ W[b]))
            G_rows.append({"n": n, "a": a, "b": b, "cos": float(cos)})
T, Pp, G = pd.DataFrame(T_rows), pd.DataFrame(P_rows), pd.DataFrame(G_rows)
T.to_csv(os.path.join(HERE, "result3b_transfer.csv"), index=False)
Pp.to_csv(os.path.join(HERE, "result3b_pooled.csv"), index=False)
G.to_csv(os.path.join(HERE, "result3b_angles.csv"), index=False)
pd.set_option("display.width", 220)
print(Pp.groupby("n")[["oracle", "pooled_LOFO", "equal_linear", "cauchy_mc"]].mean().round(3).to_string())
for n in [500, 1000]:
    t = T[T.n == n]
    own = t[t.weights_from == t.tested_on].power.mean()
    other = t[t.weights_from != t.tested_on].power.mean()
    print(f"n={n}: oracle weights on their OWN family {own:.3f}; on OTHER families {other:.3f}")
    g = G[G.n == n].cos
    print(f"        cosine between oracle weight directions: median {g.median():.2f}, quartiles {g.quantile(.25):.2f} to {g.quantile(.75):.2f}, share below 0.5: {(g < 0.5).mean():.2f}")
