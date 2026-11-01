"""result8_legoft.py -- the shipped LEGofT rule (legoft() in ebrahim.gof: frozen weights on eleven tests, weighted Cauchy
statistic) against the equal-weight Cauchy combination of the same eleven tests, on the stored archive.

Both statistics are calibrated on the archive's null reference set (95th percentile), as everywhere in the paper.
Reported: false-alarm rate on the null test set, average power over the alternatives, and false-alarm rate with
k = 0, 1, 5, 10 corrupted covariate values (corruption study). Stored p-values only; no simulation.
Output: result8_legoft.csv."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(HERE, "..", "pvalues", "archive")
A = 0.05
# frozen weights shipped in ebrahim.gof (.LEGOFT_WEIGHTS)
W = {"HL": 0.0098, "HL-equalwidth": 0.0258, "Pigeon-Heyse": 0.0099, "EF": 0.0104, "DEF.poly2": 0.0409,
     "DEF.poly3": 0.0401, "DEF.stukel": 0.0346, "Tsiatis": 0.1776, "Xie": 0.1994, "A.poly": 0.2255, "A.spline": 0.2260}
COLS = list(W)
wv = np.array([W[c] for c in COLS])
cauchy = lambda P: np.tan((0.5 - np.clip(P, 1e-15, 1 - 1e-15)) * np.pi)
RULES = {"LEGofT (frozen weights)": lambda P: cauchy(P) @ wv, "Equal weights": lambda P: cauchy(P).mean(axis=1)}


def P_of(f):
    return pd.read_csv(f)[COLS].dropna().to_numpy()


rows = []
for n in [500, 1000]:
    ref = P_of(os.path.join(S1, f"conformal_null_n{n}_pvalues.csv"))
    crit = {r: np.quantile(fn(ref), 1 - A) for r, fn in RULES.items()}
    null = P_of(os.path.join(S1, f"test_null_n{n}_pvalues.csv"))
    alts = sorted(f for f in os.listdir(S1) if f.startswith("test_") and f.endswith(f"_n{n}_pvalues.csv")
                  and not f.startswith("test_null_"))
    for r, fn in RULES.items():
        row = {"n": n, "rule": r, "size": (fn(null) > crit[r]).mean(),
               "power": np.mean([(fn(P_of(os.path.join(S1, f))) > crit[r]).mean() for f in alts]), "alternatives": len(alts)}
        for k in [0, 1, 5, 10]:
            d = pd.read_csv(os.path.join(HERE, "..", "pvalues", "corruption", f"corrupt_k{k}_n{n}.csv"))
            if k == 0:
                d = d[d.rep <= 1000]
            row[f"corrupt_k{k}"] = (fn(d[COLS].dropna().to_numpy()) > crit[r]).mean()
        rows.append(row)
T = pd.DataFrame(rows)
T.to_csv(os.path.join(HERE, "result8_legoft.csv"), index=False)
print(T.round(3).to_string(index=False))
