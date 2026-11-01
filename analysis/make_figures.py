"""make_figures.py -- every figure of the panel paper, drawn from the result files (no number typed by hand).
Vector PDFs in analysis/figures/, drawn at the journal text width (34 pc = 5.67 in) so that they are placed unscaled and
the lettering stays at 7-9 pt. Colours: a validated categorical palette (colour-blind safe on adjacent pairs); every
series also has its own marker and a direct label, so nothing depends on colour alone and the figures survive greyscale.
Re-run after any change to the results."""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
R = lambda f: pd.read_csv(os.path.join(HERE, f))

W = 5.67                                   # text width in inches
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e3df"
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.axisbelow": True, "lines.linewidth": 1.6, "lines.markersize": 4.5,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "figure.dpi": 150})
PANEL = {"CLASSIC": "Classical panel (6 tests)", "ONE-PER-GROUP": "One test per group (5)", "FULL": "Full panel (14 tests)"}
NAME = {"HL": "Hosmer–Lemeshow", "HL-equalwidth": "HL, equal-width groups", "Pigeon-Heyse": "Pigeon–Heyse",
        "Tsiatis": "Tsiatis", "Xie": "Xie et al.", "EF": "EF", "DEF.poly2": "EDGE, quadratic", "DEF.poly3": "EDGE, cubic",
        "DEF.stukel": "EDGE, Stukel shape", "Stukel": "Stukel", "Osius-Rojek": "Osius–Rojek",
        "le-Cessie": "le Cessie–van Houwelingen", "A.poly": "Covariate, polynomial", "A.spline": "Covariate, spline"}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)


def nominal(ax, se_n=2000):
    """Nominal level and the band of +-2 Monte Carlo standard errors around it."""
    se = np.sqrt(0.05 * 0.95 / se_n)
    ax.axhspan(0.05 - 2 * se, 0.05 + 2 * se, color=GRID, lw=0, zorder=0)
    ax.axhline(0.05, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)


# ---- Figure 1: false-alarm rate of the four decisions -------------------------------------------------------------------
s = R("result1_size.csv")
fig, axes = plt.subplots(1, 3, figsize=(W, 2.3), sharey=True)
DEC = [("ANY", "Any test rejects", ORANGE, "o"), ("Cauchy", "Cauchy, analytic", VIOLET, "s"),
       ("Bonferroni", "Bonferroni", AQUA, "^"), ("Cauchy-MC", "Cauchy, simulated", BLUE, "D")]
for ax, p in zip(axes, ["CLASSIC", "ONE-PER-GROUP", "FULL"]):
    g = s[s.panel == p].sort_values("n")
    nominal(ax)
    for col, lab, c, m in DEC:
        ax.plot(range(3), g[col], color=c, marker=m, label=lab, zorder=3)
    ax.set_xticks(range(3), ["200", "500", "1000"]); ax.set_xlabel("Sample size $n$")
    ax.set_title(PANEL[p], loc="left"); ax.set_ylim(0, 0.32); ax.set_xlim(-0.25, 2.25)
axes[0].set_ylabel("False-alarm rate")
h, l = axes[0].get_legend_handles_labels()
h.append(mpl.lines.Line2D([], [], color=MUTED, lw=0.8, ls=(0, (3, 2)))); l.append("nominal 0.05 ± 2 s.e.")
fig.legend(h, l, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.08), columnspacing=1.0, handlelength=1.8)
save(fig, "fig1_size.pdf")

# ---- Figure 2: power of the full combination against the best single test, alternative by alternative ------------------
pw = R("result1_power.csv")
fig, axes = plt.subplots(1, 3, figsize=(W, 2.15), sharex=True, sharey=True)
for ax, n in zip(axes, [200, 500, 1000]):
    g = pw[(pw.panel == "FULL") & (pw.n == n)]
    x, y = g["single:le-Cessie"].to_numpy(), g["Cauchy-MC"].to_numpy()
    ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    above = y >= x - 1e-12
    ax.scatter(x[above], y[above], s=18, color=BLUE, edgecolor="white", lw=0.6, zorder=3, label="combination at least as good")
    ax.scatter(x[~above], y[~above], s=18, color=ORANGE, marker="s", edgecolor="white", lw=0.6, zorder=3,
               label="single test better")
    ax.set_title(f"$n={n}$", loc="left")
    ax.text(0.97, 0.05, f"{above.sum()} of {len(g)} on or\nabove the line", ha="right", va="bottom", transform=ax.transAxes,
            fontsize=7, color=INK)
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xlabel("le Cessie–van Houwelingen")
axes[0].set_ylabel("Full panel, Cauchy simulated")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.08))
save(fig, "fig2_power_scatter.pdf")

# ---- Figure 3: null correlation of the fourteen tests, ordered by the five groups --------------------------------------
Rm = R("result2_corr_null_n1000.csv").set_index("Unnamed: 0")
GROUPS = [["HL", "HL-equalwidth", "Pigeon-Heyse", "EF"], ["Tsiatis", "Xie"], ["DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel"],
          ["Osius-Rojek"], ["le-Cessie", "A.poly", "A.spline"]]
order = [t for gr in GROUPS for t in gr]
M = Rm.loc[order, order].abs().to_numpy()
fig, ax = plt.subplots(figsize=(W * 0.78, W * 0.62))
im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1)
ax.grid(False)
for i in range(len(order)):
    for j in range(len(order)):
        ax.text(j, i, f"{M[i, j]:.2f}".lstrip("0") if round(M[i, j], 2) < 1 else "1", ha="center", va="center", fontsize=5.6,
                color="white" if M[i, j] > 0.6 else INK)
pos = 0
for gr in GROUPS:
    ax.add_patch(Rectangle((pos - 0.5, pos - 0.5), len(gr), len(gr), fill=False, ec=ORANGE, lw=1.4))
    pos += len(gr)
ax.set_xticks(range(len(order)), [NAME[t] for t in order], rotation=50, ha="right", fontsize=7)
ax.set_yticks(range(len(order)), [NAME[t] for t in order], fontsize=7)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("Absolute Spearman correlation under the null", fontsize=7); cb.ax.tick_params(labelsize=6.5); cb.outline.set_visible(False)
save(fig, "fig4_correlation.pdf")

# ---- Figure 4: other combination rules against the Cauchy combination, with 95% bootstrap intervals -------------------
D = R("result5_diff.csv")
RULES = [("HMP", "Harmonic mean $p$"), ("Fisher", "Fisher"), ("minP", "Minimum $p$"), ("Stouffer", "Stouffer"),
         ("Bonferroni", "Bonferroni, not calibrated"), ("Quadratic", "Quadratic $z^{\\top}S^{-1}z$")]
LIM = (-9, 6)
fig, axes = plt.subplots(1, 2, figsize=(W, 2.4), sharey=True)
for ax, p in zip(axes, ["CLASSIC", "FULL"]):
    ax.axvline(0, color=MUTED, lw=0.8)
    for off, n, c, m in [(-0.15, 500, AQUA, "o"), (0.15, 1000, BLUE, "D")]:
        for i, (r, _) in enumerate(RULES):
            d = D[(D.n == n) & (D.panel == p) & (D.comparison == f"{r} - Cauchy")].iloc[0]
            v, lo, hi = 100 * d["diff"], 100 * d.lo, 100 * d.hi
            yv = i + off
            if v < LIM[0]:
                ax.annotate(f"{v:.0f}", xy=(LIM[0] + 0.1, yv), xytext=(LIM[0] + 1.6, yv), va="center", fontsize=6.5,
                            color=c, arrowprops=dict(arrowstyle="->", color=c, lw=0.9))
                continue
            ax.plot([lo, hi], [yv, yv], color=c, lw=1.4, solid_capstyle="round")
            ax.plot(v, yv, marker=m, color=c, mec="white", mew=0.6, ms=5,
                    label=f"$n={n}$" if i == 0 else None, ls="none")
    ax.set_yticks(range(len(RULES)), [lab for _, lab in RULES])
    ax.set_xlim(*LIM); ax.set_xlabel("Power minus Cauchy combination (points)")
    ax.set_title(PANEL[p], loc="left"); ax.grid(axis="y", visible=False)
axes[0].invert_yaxis()
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.55, 1.06))
save(fig, "fig3_rules.pdf")

# ---- Figure 5: the limit of weighting ---------------------------------------------------------------------------------------
c3, pool, tr = R("result3_ceiling.csv"), R("result3b_pooled.csv"), R("result3b_transfer.csv")
rows = []
for n in [500, 1000]:
    off = tr[(tr.n == n) & (tr.weights_from != tr.tested_on)].power.mean()
    q = pool[pool.n == n]
    rows.append({"n": n, "Oracle weights (alternative known)": c3[c3.n == n].oracle_power.mean(),
                 "Equal weights, Cauchy": q.cauchy_mc.mean(), "Pooled from other alternatives": q.pooled_LOFO.mean(),
                 "Equal weights, linear sum": q.equal_linear.mean(), "Weights from one other alternative": off})
T = pd.DataFrame(rows).set_index("n")
labels = list(T.columns)
fig, ax = plt.subplots(figsize=(W, 1.95))
for i, lab in enumerate(labels):
    a, b = T.loc[500, lab], T.loc[1000, lab]
    ax.plot([a, b], [i, i], color=GRID, lw=3, solid_capstyle="round", zorder=1)
    ax.plot(a, i, "o", color=AQUA, mec="white", mew=0.6, ms=6, zorder=3, label="$n=500$" if i == 0 else None)
    ax.plot(b, i, "D", color=BLUE, mec="white", mew=0.6, ms=5.5, zorder=3, label="$n=1000$" if i == 0 else None)
    ax.text(a - 0.012, i, f"{a:.3f}", ha="right", va="center", fontsize=6.8, color=MUTED)
    ax.text(b + 0.012, i, f"{b:.3f}", ha="left", va="center", fontsize=6.8, color=MUTED)
ax.set_yticks(range(len(labels)), labels); ax.invert_yaxis(); ax.grid(axis="y", visible=False)
ax.set_xlim(0.40, 0.97); ax.set_xlabel("Power averaged over 17 alternatives")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, ncol=2)
save(fig, "fig5_weighting.pdf")

# ---- Figure 6: false alarms of combined decisions under corrupted covariate values ------------------------------------
P4 = R("result4_panels.csv")
PN = [("FULL", "Full", ORANGE, "o"), ("ONE-PER-GROUP", "One per group", VIOLET, "s"),
      ("LEARNED-ROBUST", "Screened", AQUA, "^"), ("PROTECTED", "Protected", BLUE, "D")]
fig, axes = plt.subplots(1, 2, figsize=(W, 2.3), sharey=True)
for ax, n in zip(axes, [500, 1000]):
    nominal(ax, 1000)
    ends = []
    for p, lab, c, m in PN:
        g = P4[(P4.n == n) & (P4.panel == p) & (P4.k.astype(str) != "0all")].copy()
        g["k"] = g.k.astype(int); g = g.sort_values("k")
        ax.plot(g.k, g["Cauchy-MC"], color=c, marker=m, zorder=3)
        ends.append([g["Cauchy-MC"].iloc[-1], lab])
    ends.sort(key=lambda e: e[0])
    for a, b in zip(ends, ends[1:]):          # keep at least 0.045 between neighbouring labels
        b[0] = max(b[0], a[0] + 0.045)
    for y, lab in ends:
        ax.text(10.9, y, lab, va="center", fontsize=7, color=INK, clip_on=False)
    ax.set_xticks([0, 1, 5, 10]); ax.set_xlim(-0.4, 10.3); ax.set_ylim(0, 0.7)
    ax.set_xlabel("Corrupted covariate values $k$"); ax.set_title(f"$n={n}$", loc="left")
axes[0].set_ylabel("False-alarm rate, Cauchy simulated")
for t in axes[1].texts:
    t.set_clip_on(False)
fig.subplots_adjust(right=0.84, wspace=0.62)
save(fig, "fig7_corrupt_panels.pdf")

# ---- Figure 7: false alarms of single tests under corruption (heat map) -------------------------------------------------
S4 = R("result4_single.csv")
S4 = S4[S4.k.astype(str) != "0all"].copy(); S4["k"] = S4.k.astype(int)
order = S4[(S4.n == 500) & (S4.k == 10)].sort_values("false_alarm").test.tolist()
cols = [(n, k) for n in [500, 1000] for k in [0, 1, 5, 10]]
M = np.array([[S4[(S4.n == n) & (S4.k == k) & (S4.test == t)].false_alarm.iloc[0] for n, k in cols] for t in order])
fig, ax = plt.subplots(figsize=(W * 0.80, 3.1))
im = ax.imshow(M, cmap="Oranges", vmin=0, vmax=0.65, aspect="auto")
ax.grid(False)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6, color="white" if M[i, j] > 0.4 else INK)
ax.axvline(3.5, color="white", lw=3)
ax.set_xticks(range(len(cols)), [f"$k={k}$" for _, k in cols], fontsize=7)
ax.set_yticks(range(len(order)), [NAME[t] for t in order], fontsize=7)
ax.tick_params(length=0)
ax.text(1.5, -1.2, "$n=500$", ha="center", fontsize=7.5); ax.text(5.5, -1.2, "$n=1000$", ha="center", fontsize=7.5)
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
cb.set_label("False-alarm rate (nominal 0.05)", fontsize=7); cb.ax.tick_params(labelsize=6.5); cb.outline.set_visible(False)
save(fig, "fig6_corrupt_single.pdf")

# ---- Figure 4: false alarms in two other null designs, full panel ------------------------------------------------------
D6 = R("result6_designs.csv")
fig, axes = plt.subplots(1, 2, figsize=(W, 2.5), sharey=True)
DEC6 = DEC[1:] + [("Cauchy-MC_wrong_design", "Cauchy, simulated on the original design", BLUE, "D")]
TITLE6 = {"rare": "Rare events (13% events, 13 tests)", "six": "Six covariates (59% events, 14 tests)"}
for ax, d in zip(axes, ["rare", "six"]):
    g = D6[(D6.design == d) & (D6.panel == "FULL")].sort_values("n")
    nominal(ax)
    for col, lab, c, m in DEC6:
        wrong = col.endswith("wrong_design")
        ax.plot(range(2), g[col], color=c, marker=m, label=lab, zorder=3, ls=(0, (2, 1.5)) if wrong else "-",
                mfc="white" if wrong else c)
    # the any-test rule is far off this scale: print its rates at the top, as Figure 3 does for off-scale values
    for x, v in enumerate(g["ANY"]):
        ax.annotate(f"{v:.2f}", xy=(x, 0.108), xytext=(x, 0.094), ha="center", fontsize=6.5, color=ORANGE,
                    arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9))
    ax.set_xticks(range(2), ["200", "500"]); ax.set_xlabel("Sample size $n$")
    ax.set_title(TITLE6[d], loc="left"); ax.set_ylim(0, 0.11); ax.set_xlim(-0.3, 1.3)
axes[0].set_ylabel("False-alarm rate")
h, l = axes[0].get_legend_handles_labels()
h.insert(0, mpl.lines.Line2D([], [], color=ORANGE, marker=r"$\uparrow$", ls="none", ms=6)); l.insert(0, "Any test rejects (off scale)")
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.12), columnspacing=1.2, handlelength=2.2)
save(fig, "fig8_designs.pdf")
print("figures written to", OUT, ":", sorted(os.listdir(OUT)))
