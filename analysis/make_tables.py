"""make_tables.py -- writes every table of the panel paper as a LaTeX fragment in paper/tables/, straight from the result
files. No number in the manuscript tables is typed by hand; re-run this after any change to the results."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)
R = lambda f: pd.read_csv(os.path.join(HERE, f))
f3 = lambda x: f"{x:.3f}"
PANEL_NAME = {"CLASSIC": "Classical (6)", "FULL": "Full (14)", "ONE-PER-GROUP": "One per group (5)",
              "PROTECTED": "Protected (4)", "LEARNED-ROBUST": "Screened (7 or 8)"}
TEST_NAME = {"HL": "Hosmer--Lemeshow", "HL-equalwidth": "HL, equal-width groups", "Pigeon-Heyse": "Pigeon--Heyse",
             "Tsiatis": "Tsiatis", "Xie": "Xie et al.", "EF": "EF", "DEF.poly2": "EDGE, quadratic",
             "DEF.poly3": "EDGE, cubic", "DEF.stukel": "EDGE, Stukel shape", "Stukel": "Stukel",
             "Osius-Rojek": "Osius--Rojek", "le-Cessie": "le Cessie--van Houwelingen",
             "A.poly": "Covariate, polynomial", "A.spline": "Covariate, spline"}


def write(name, body):
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
        fh.write(body)


# ---- Table: size of the four decisions -----------------------------------------------------------------------------
s = R("result1_size.csv")
rows = []
for p in ["CLASSIC", "ONE-PER-GROUP", "FULL"]:
    for _, r in s[s.panel == p].sort_values("n").iterrows():
        rows.append(f"{PANEL_NAME[p]} & {int(r.n)} & {f3(r.ANY)} & {f3(r.Bonferroni)} & {f3(r.Cauchy)} & {f3(r['Cauchy-MC'])} & "
                    f"{f3(r.single_min)}--{f3(r.single_max)} \\\\")
    rows.append("\\addlinespace")
write("tab_size.tex", "\n".join([
    "\\begin{tabular}{lrccccc}", "\\toprule",
    "Panel & $n$ & Any test & Bonferroni & Cauchy, analytic & Cauchy, simulated & Single tests \\\\", "\\midrule",
    *rows[:-1], "\\bottomrule", "\\end{tabular}"]))

# ---- Table: power of the calibrated combination against the best single test -------------------------------------------
pw = R("result1_power.csv")
rows, bests = [], set()
for p in ["CLASSIC", "ONE-PER-GROUP", "FULL"]:
    for n in sorted(pw.n.unique()):
        g = pw[(pw.panel == p) & (pw.n == n)]
        sc = {c[7:]: g[c].mean() for c in g.columns if c.startswith("single:") and g[c].notna().all()}
        b = max(sc, key=sc.get); bests.add(b)
        wins = int((g["Cauchy-MC"] >= g["single:" + b] - 1e-12).sum())
        rows.append(f"{PANEL_NAME[p]} & {n} & {f3(g['Cauchy-MC'].mean())} & {f3(g.Bonferroni.mean())} & "
                    f"{f3(sc[b])} & {wins}/{len(g)} \\\\")
    rows.append("\\addlinespace")
# the best fixed single test is le Cessie--van Houwelingen in every row, so it is named in the caption, not in a column
assert bests == {"le-Cessie"}, bests
write("tab_power.tex", "\n".join([
    "\\begin{tabular}{lrcccc}", "\\toprule",
    " & & \\multicolumn{2}{c}{Combination} & Best single & Cauchy at least \\\\",
    "\\cmidrule(lr){3-4}",
    "Panel & $n$ & Cauchy, simulated & Bonferroni & test & as good \\\\",
    "\\midrule", *rows[:-1], "\\bottomrule", "\\end{tabular}"]))

# ---- Table: redundancy ---------------------------------------------------------------------------------------------
m = R("result2_meff.csv"); m["family"] = m.family.fillna("null")   # pandas reads the word null as missing
rows = []
for n in sorted(m.n.unique()):
    g = m[m.n == n]; nul = g[g.family == "null"].iloc[0]; alt = g[g.family != "null"]
    rows.append(f"{n} & {nul.M_eff_FULL14:.0f} & {alt.M_eff_FULL14.min():.0f}--{alt.M_eff_FULL14.max():.0f} & "
                f"{nul.M_eff_CLASSIC6:.0f} & {alt.M_eff_CLASSIC6.min():.0f}--{alt.M_eff_CLASSIC6.max():.0f} & {nul.median_abs_r:.2f} \\\\")
write("tab_meff.tex", "\n".join([
    "\\begin{tabular}{rccccc}", "\\toprule",
    " & \\multicolumn{2}{c}{Full panel (14 tests)} & \\multicolumn{2}{c}{Classical panel (6 tests)} & \\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
    "$n$ & Null & Alternatives & Null & Alternatives & Median $|r|$, null \\\\", "\\midrule",
    *rows, "\\bottomrule", "\\end{tabular}"]))

# ---- Table: the ceiling of reweighting ---------------------------------------------------------------------------------
c = R("result3_ceiling.csv"); pl = R("result3b_pooled.csv"); an = R("result3b_angles.csv"); tr = R("result3b_transfer.csv")
rows = []
for n in sorted(c.n.unique()):
    g, h = c[c.n == n], pl[pl.n == n]
    t = tr[(tr.n == n) & (tr.weights_from != tr.tested_on)]
    a = an[an.n == n].cos
    rows.append(f"{n} & {len(g)} & {f3(g.eff_equal.mean())} & {f3(g.oracle_power.mean())} & {f3(t.power.mean())} & "
                f"{f3(h.pooled_LOFO.mean())} & {f3(g.equal_linear_power.mean())} & {f3(g.cauchy_mc_power.mean())} & {a.median():.2f} \\\\")
write("tab_ceiling.tex", "\n".join([
    "\\begin{tabular}{rcccccccc}", "\\toprule",
    " & & & \\multicolumn{3}{c}{Weights chosen from alternatives} & \\multicolumn{2}{c}{Equal weights} & \\\\",
    "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}",
    "$n$ & Families & Efficiency & Same & Another & All others & Linear & Cauchy & Median cosine \\\\", "\\midrule",
    *rows, "\\bottomrule", "\\end{tabular}"]))

# ---- Tables: robustness to corrupted records -----------------------------------------------------------------------------
if os.path.exists(os.path.join(HERE, "result4_single.csv")):
    S, P, W = R("result4_single.csv"), R("result4_panels.csv"), R("result4_power.csv")
    ns = sorted(S.n.unique()); ks = [0, 1, 5, 10]
    head = " & ".join(f"\\multicolumn{{{len(ks)}}}{{c}}{{$n={n}$}}" for n in ns)
    cmid = "".join(f"\\cmidrule(lr){{{2 + i * len(ks)}-{1 + (i + 1) * len(ks)}}}" for i in range(len(ns)))
    kk = " & ".join(f"$k={k}$" for _ in ns for k in ks)
    order = S[(S.n == ns[0]) & (S.k.astype(str) == "10")].sort_values("false_alarm").test.tolist()
    rows = [TEST_NAME[t] + " & " + " & ".join(f3(S[(S.n == n) & (S.k.astype(str) == str(k)) & (S.test == t)].false_alarm.iloc[0])
                                             for n in ns for k in ks) + " \\\\" for t in order]
    write("tab_corrupt_single.tex", "\n".join([
        "\\begin{tabular}{l" + "c" * len(ns) * len(ks) + "}", "\\toprule", " & " + head + " \\\\", cmid,
        "Test & " + kk + " \\\\", "\\midrule", *rows, "\\bottomrule", "\\end{tabular}"]))
    rows = []
    for p in ["FULL", "ONE-PER-GROUP", "PROTECTED", "LEARNED-ROBUST"]:
        cells = []
        for n in ns:
            g = P[(P.n == n) & (P.panel == p)]
            cells += [f3(g[g.k.astype(str) == str(k)]["Cauchy-MC"].iloc[0]) for k in ks]
            cells.append(f3(W[(W.n == n) & (W.panel == p)].power.mean()))
        rows.append(PANEL_NAME[p] + " & " + " & ".join(cells) + " \\\\")
    k5 = len(ks) + 1
    head = " & ".join(f"\\multicolumn{{{k5}}}{{c}}{{$n={n}$}}" for n in ns)
    cmid = "".join(f"\\cmidrule(lr){{{2 + i * k5}-{1 + (i + 1) * k5}}}" for i in range(len(ns)))
    kk = " & ".join(" & ".join(f"$k={k}$" for k in ks) + " & Power" for _ in ns)
    write("tab_corrupt_panels.tex", "\n".join([
        "\\begin{tabular}{l" + "c" * len(ns) * k5 + "}", "\\toprule", " & " + head + " \\\\", cmid,
        "Panel & " + kk + " \\\\", "\\midrule", *rows, "\\bottomrule", "\\end{tabular}"]))
# ---- Table: other combination rules, with bootstrap intervals (result5_rules_ci.py) ------------------------------------
if os.path.exists(os.path.join(HERE, "result5_power.csv")):
    S5, P5, D5 = R("result5_size.csv"), R("result5_power.csv"), R("result5_diff.csv")
    RULE_NAME = {"Cauchy": "Cauchy", "HMP": "Harmonic mean $p$", "minP": "Minimum $p$", "Fisher": "Fisher",
                 "Stouffer": "Stouffer", "Quadratic": "Quadratic $z^\\top S^{-1}z$", "Bonferroni": "Bonferroni (not calibrated)"}
    pts = lambda x: f"{100 * x:+.1f}".replace("-", "$-$")
    rows = []
    for p in ["CLASSIC", "FULL"]:
        rows.append(f"\\multicolumn{{7}}{{l}}{{\\emph{{{PANEL_NAME[p]}}}}} \\\\")
        for rule in RULE_NAME:
            cells = []
            for n in [500, 1000]:
                sz = S5[(S5.n == n) & (S5.panel == p) & (S5.rule == rule)].iloc[0]
                pw = P5[(P5.n == n) & (P5.panel == p) & (P5.rule == rule)].iloc[0]
                if rule == "Cauchy":
                    df = "--"
                else:
                    d = D5[(D5.n == n) & (D5.panel == p) & (D5.comparison == f"{rule} - Cauchy")].iloc[0]
                    df = f"{pts(d['diff'])} ({pts(d.lo)}, {pts(d.hi)})"
                cells += [f3(sz["size"]), f3(pw.power), df]
            rows.append(RULE_NAME[rule] + " & " + " & ".join(cells) + " \\\\")
        rows.append("\\addlinespace")
    write("tab_rules.tex", "\n".join([
        "\\begin{tabular}{lcclccl}", "\\toprule",
        " & \\multicolumn{3}{c}{$n=500$} & \\multicolumn{3}{c}{$n=1000$} \\\\",
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
        "Rule & Size & Power & Difference (95\\% CI) & Size & Power & Difference (95\\% CI) \\\\", "\\midrule",
        *rows[:-1], "\\bottomrule", "\\end{tabular}"]))
print("tables written to", OUT, ":", sorted(os.listdir(OUT)))
