"""make_tables.py -- writes every table of the panel paper as a LaTeX fragment in paper/tables/, straight from the result
files. No number in the manuscript tables is typed by hand; re-run this after any change to the results."""
import os
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "tables")
os.makedirs(OUT, exist_ok=True)
R = lambda f: pd.read_csv(os.path.join(HERE, f))
# one rounding rule everywhere: half up on the decimal value, so 0.0455 prints as 0.046 (f"{x:.3f}" rounds the binary value)
f3 = lambda x: str(Decimal(repr(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))
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
    pts = lambda x: ("+" if x >= 0 else "$-$") + str(Decimal(repr(abs(100 * float(x)))).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
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
# ---- Table: power of each panel alternative by alternative (cost of the protected panel) --------------------------------
if os.path.exists(os.path.join(HERE, "result4_power.csv")):
    W4 = R("result4_power.csv")
    cols = ["FULL", "ONE-PER-GROUP", "PROTECTED", "LEARNED-ROBUST"]
    piv = {n: W4[W4.n == n].pivot(index="family", columns="panel", values="power") for n in [500, 1000]}
    fams = piv[500].assign(loss=piv[500].FULL - piv[500].PROTECTED).sort_values("loss", ascending=False).index
    rows = []
    for f in fams:
        cells = [f3(piv[n].loc[f, c]) for n in [500, 1000] for c in cols]
        rows.append(f"\\texttt{{{f.replace('_', chr(92) + '_')}}} & " + " & ".join(cells) + " \\\\")
    rows.append("\\midrule")
    rows.append("Average & " + " & ".join(f3(piv[n][c].mean()) for n in [500, 1000] for c in cols) + " \\\\")
    short = {"FULL": "Full", "ONE-PER-GROUP": "One/group", "PROTECTED": "Protected", "LEARNED-ROBUST": "Screened"}
    write("tab_panel_power.tex", "\n".join([
        "\\begin{tabular}{lcccccccc}", "\\toprule",
        " & \\multicolumn{4}{c}{$n=500$} & \\multicolumn{4}{c}{$n=1000$} \\\\", "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
        "Alternative & " + " & ".join(short[c] for _ in range(2) for c in cols) + " \\\\", "\\midrule",
        *rows, "\\bottomrule", "\\end{tabular}"]))

# ---- Tables: the extra simulations of run_extra.R (result6_extra.py) ----------------------------------------------------
if os.path.exists(os.path.join(HERE, "result6_designs.csv")):
    C6, D6 = R("result6_corrupt.csv"), R("result6_designs.csv")
    BS6, BP6 = R("result6_boot_size.csv"), R("result6_boot_power.csv")
    # other corruption mechanisms: single tests and panels, k = 0, 1, 5, 10 for each mechanism
    ks = [0, 1, 5, 10]
    order = C6[(C6.mechanism == "x2") & (C6.k == 10)].iloc[0][list(TEST_NAME)].sort_values().index.tolist()
    cell = lambda m, k, col: f3(C6[(C6.mechanism == m) & (C6.k == k)][col].iloc[0])
    rows = [TEST_NAME[t] + " & " + " & ".join(cell(m, k, t) for m in ["flipy", "x2"] for k in ks) + " \\\\" for t in order]
    rows.append("\\midrule")
    rows += [f"Panel: {PANEL_NAME[p]} & " + " & ".join(cell(m, k, "panel:" + p) for m in ["flipy", "x2"] for k in ks) + " \\\\"
             for p in ["FULL", "CLASSIC", "ONE-PER-GROUP", "PROTECTED"]]
    kk = " & ".join(f"$k={k}$" for _ in range(2) for k in ks)
    write("tab_extra_corrupt.tex", "\n".join([
        "\\begin{tabular}{lcccccccc}", "\\toprule",
        " & \\multicolumn{4}{c}{Flipped responses} & \\multicolumn{4}{c}{$x_2$ multiplied by 4} \\\\",
        "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}", "Test & " + kk + " \\\\", "\\midrule",
        *rows, "\\bottomrule", "\\end{tabular}"]))
    # other null designs
    DESIGN_NAME = {"rare": "Rare events", "six": "Six covariates"}
    rows = []
    for _, r in D6.iterrows():
        rows.append(f"{DESIGN_NAME[r.design]} & {int(r.n)} & {PANEL_NAME[r.panel].split(' (')[0]} ({int(r.tests)}) & "
                    f"{f3(r.ANY)} & {f3(r.Bonferroni)} & {f3(r.Cauchy)} & {f3(r['Cauchy-MC'])} & "
                    f"{f3(r['Cauchy-MC_wrong_design'])} & {f3(r.single_min)}--{f3(r.single_max)} & {r.M_eff:.0f} \\\\")
    write("tab_designs.tex", "\n".join([
        "\\begin{tabular}{lrlccccccc}", "\\toprule",
        " & & & & & \\multicolumn{3}{c}{Cauchy} & & \\\\", "\\cmidrule(lr){6-8}",
        "Design & $n$ & Panel & Any test & Bonferroni & analytic & own design & original design & Single tests & $M_{\\mathrm{eff}}$ \\\\",
        "\\midrule", *rows, "\\bottomrule", "\\end{tabular}"]))
    # bootstrap calibration: size of four rules on two panels, and power against the fixed reference calibration
    RN = {"cauchy": "Cauchy", "hmp": "Harmonic mean $p$", "fisher": "Fisher", "minp": "Minimum $p$"}
    # size in four settings: the original design at n = 200 and 500 (run_extra.R) and the two other designs at n = 200
    # (run_extra2.R, result7_extra2.py)
    BD7 = R("result7_boot_designs.csv") if os.path.exists(os.path.join(HERE, "result7_boot_designs.csv")) else None
    def boot_size(setting, rule):
        if setting.startswith("orig"):
            return BS6[(BS6.n == int(setting[4:])) & (BS6.rule == rule)]["size"].iloc[0]
        return BD7[(BD7.design == setting) & (BD7.rule == rule)]["size"].iloc[0]
    settings = ["orig200", "orig500"] + (["rare", "six"] if BD7 is not None else [])
    heads = {"orig200": "$n=200$", "orig500": "$n=500$", "rare": "Rare, 200", "six": "Six cov., 200"}
    rows = []
    for p in ["classic", "full"]:
        rows.append(f"\\multicolumn{{{len(settings) + 1}}}{{l}}{{\\emph{{{PANEL_NAME[p.upper()]}}}}} \\\\")
        for rule in RN:
            rows.append(RN[rule] + " & " + " & ".join(f3(boot_size(st, f"{p}_{rule}")) for st in settings) + " \\\\")
    write("tab_boot_size.tex", "\n".join([
        "\\begin{tabular}{l" + "c" * len(settings) + "}", "\\toprule",
        " & \\multicolumn{2}{c}{Original design}" + (" & \\multicolumn{2}{c}{Other designs}" if BD7 is not None else "") + " \\\\",
        "\\cmidrule(lr){2-3}" + ("\\cmidrule(lr){4-5}" if BD7 is not None else ""),
        "Rule & " + " & ".join(heads[st] for st in settings) + " \\\\", "\\midrule", *rows, "\\bottomrule", "\\end{tabular}"]))
    # power: the two calibrations on the same alternative data sets, as a paired comparison
    ALT_NAME = {"quad": "Missing $x_1^2$ term", "quad2": "Missing $x_2^2$ term", "cloglog": "Wrong link (cloglog)",
                "interaction": "Missing interaction"}
    BE7 = R("result7_boot_paired.csv")
    rows = []
    for (n, B), g in BE7.groupby(["n", "B"], sort=True):
        rows.append(f"\\multicolumn{{7}}{{l}}{{\\emph{{$n={n}$, $B={B}$, {g.reps.iloc[0]} data sets per alternative}}}} \\\\")
        for _, r in g.iterrows():
            rows.append(f"{ALT_NAME[r.alternative]} & {PANEL_NAME[r.panel]} & {f3(r.power_bootstrap)} & {f3(r.power_fixed)} & "
                        f"{int(r.bootstrap_only)} & {int(r.fixed_only)} & {pts(r['diff'])} ({pts(r.lo)}, {pts(r.hi)}) \\\\")
    write("tab_boot_power.tex", "\n".join([
        "\\begin{tabular}{llccccl}", "\\toprule",
        " & & \\multicolumn{2}{c}{Power} & \\multicolumn{2}{c}{Rejected only by} & \\\\", "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
        "Alternative & Panel & Bootstrap & Fixed & bootstrap & fixed & Difference (95\\% CI) \\\\", "\\midrule",
        *rows, "\\bottomrule", "\\end{tabular}"]))
print("tables written to", OUT, ":", sorted(os.listdir(OUT)))
