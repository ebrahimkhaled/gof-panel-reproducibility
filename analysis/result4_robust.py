"""result4_robust.py -- Result 4 of the panel paper: is a combined goodness-of-fit decision only as robust as its most
fragile input? Data: corrupt/corrupt_k{k}_n{n}.csv from run_corrupt.R (correct model; x1 multiplied by 4 in k recorded
values; same clean data sets for every k), plus the clean LEGofT archive for power.

Calibration of every combined decision: its critical value comes from the archive's independent conformal null set; false
alarms are measured on replicates 1-1000 (k = 0: clean size; k > 0: the same data sets with k corrupted values), and
clean size also on all 2,000 clean replicates (column 0all).
Panels: FULL (14 tests), ONE-PER-GROUP (HL, Xie, DEF.poly3, Osius-Rojek, le Cessie), PROTECTED (fixed in advance from the
EF and EDGE papers: HL, EF, Pigeon-Heyse, le Cessie), and LEARNED-ROBUST (data-driven: drop every test whose false-alarm
rate at k = 5 on replicates 1-500 exceeds 0.10; the rule is then judged on replicates 501-1000 only).
Power: each panel's calibrated Cauchy combination on the clean archive alternatives (critical value from the archive's
conformal null set). Outputs: result4_single.csv, result4_panels.csv, result4_power.csv."""
import os
import re
import numpy as np
import pandas as pd
FAMILY = lambda fname: re.sub(r"^(test|train|conformal)_|_n[0-9]+_pvalues[.]csv$", "", os.path.basename(fname))   # family names can contain underscores

HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(HERE, "..", "pvalues", "archive")
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel",
         "Osius-Rojek", "le-Cessie", "A.poly", "A.spline"]
PANELS = {"FULL": TESTS, "ONE-PER-GROUP": ["HL", "Xie", "DEF.poly3", "Osius-Rojek", "le-Cessie"],
          "PROTECTED": ["HL", "EF", "Pigeon-Heyse", "le-Cessie"]}
A = 0.05


def cstat(P):
    return np.tan((0.5 - np.clip(P, 1e-12, 1 - 1e-12)) * np.pi).mean(axis=1)


def load(k, n):
    return pd.read_csv(os.path.join(HERE, "..", "pvalues", "corruption", f"corrupt_k{k}_n{n}.csv"))


single, panels, power = [], [], []
for n in [m for m in (500, 1000) if all(os.path.exists(os.path.join(HERE, "..", "pvalues", "corruption", f"corrupt_k{k}_n{m}.csv")) and len(load(k, m)) >= (2000 if k == 0 else 1000) for k in (0, 1, 5, 10))]:
    cells = {k: load(k, n).query("rep <= 1000").sort_values("rep") for k in [0, 1, 5, 10]}
    cells["0all"] = load(0, n)   # clean size on all 2,000 clean replicates
    for k, d in cells.items():
        for t in TESTS:
            single.append({"n": n, "k": k, "test": t, "false_alarm": float((d[t].dropna() <= A).mean())})
    # learned-robust panel: chosen on reps 1-500 at k = 5
    tr = cells[5].query("rep <= 500")
    learned = [t for t in TESTS if (tr[t].dropna() <= A).mean() <= 0.10]
    allp = dict(PANELS, **{"LEARNED-ROBUST": learned})
    for pname, cols in allp.items():
        # critical value from the archive's independent conformal null set (2,000 data sets of the same model)
        cref = pd.read_csv(os.path.join(S1, f"conformal_null_n{n}_pvalues.csv"))[cols].dropna().to_numpy()
        crit_a = crit = np.quantile(cstat(cref), 1 - A)
        for k, d in cells.items():
            dd = d.query("rep > 500") if pname == "LEARNED-ROBUST" else d
            P = dd[cols].dropna().to_numpy()
            panels.append({"n": n, "k": k, "panel": pname, "members": len(cols), "reps": len(P),
                           "ANY": float((P.min(axis=1) <= A).mean()), "Cauchy-MC": float((cstat(P) >= crit).mean())})
        # power on the clean archive
        for f in sorted(os.listdir(S1)):
            if f.startswith("test_") and f.endswith(f"_n{n}_pvalues.csv") and "_null_" not in f:
                P = pd.read_csv(os.path.join(S1, f))[cols].dropna().to_numpy()
                if len(P):
                    power.append({"n": n, "panel": pname, "family": FAMILY(f), "power": float((cstat(P) >= crit_a).mean())})
    print(f"n={n}: learned-robust panel keeps {learned}")
S, Pn, Pw = pd.DataFrame(single), pd.DataFrame(panels), pd.DataFrame(power)
S.to_csv(os.path.join(HERE, "result4_single.csv"), index=False)
Pn.to_csv(os.path.join(HERE, "result4_panels.csv"), index=False)
Pw.to_csv(os.path.join(HERE, "result4_power.csv"), index=False)
pd.set_option("display.width", 220)
print("\nSINGLE TESTS: false-alarm rate by number of corrupted records")
print(S.pivot_table(index=["test"], columns=["n", "k"], values="false_alarm").round(3).to_string())
print("\nPANELS")
print(Pn.pivot_table(index=["panel"], columns=["n", "k"], values="Cauchy-MC").round(3).to_string())
print(Pn.pivot_table(index=["panel"], columns=["n", "k"], values="ANY").round(3).to_string())
print("\nPOWER on clean data (calibrated Cauchy), mean over families")
print(Pw.groupby(["panel", "n"]).power.mean().unstack().round(3).to_string())
