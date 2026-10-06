#!/usr/bin/env python3
"""Generate CBC submission figures from DEPOSITED data only (no re-simulation).
Fig1: ligand-protein minimal heavy-atom distance vs time (stay-in-pocket) for 7 systems.
Fig2: MM-PBSA dG (GB & PB) bar chart with population SD.
All inputs are files already in the v1.1.0 deposit (lig_dist.xvg, T_md_final_7_2026-10-05.csv).
"""
import os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(REPO, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# CID -> label (role noted in legend)
SYS = [
    ("25022668", "GSK650394 (pos. ctrl)"),
    ("5280445", "Luteolin"),
    ("445154", "Resveratrol"),
    ("5280863", "Kaempferol"),
    ("5281708", "Daidzein"),
    ("439246", "Naringenin"),
    ("164676", "Tanshinone IIA"),
]
# colors
COL = ["#111111", "#1f77b4", "#2ca02c", "#9467bd", "#8c564b", "#ff7f0e", "#d62728"]

# ---------------- Fig 1: lig_dist trajectory ----------------
plt.figure(figsize=(7.2, 4.6), dpi=300)
for (cid, lab), c in zip(SYS, COL):
    f = os.path.join(REPO, "analysis", "trajectories", cid, "lig_dist.xvg")
    t, y = [], []
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            p = line.split()
            t.append(float(p[0]) / 1000.0)  # ps -> ns
            y.append(float(p[1]))
    plt.plot(t, y, lw=1.1, color=c, label=lab)
plt.axhline(0.40, ls="--", lw=1.0, color="#888888", label="0.40 nm threshold")
plt.xlabel("Time (ns)")
plt.ylabel("Ligand–protein min heavy-atom distance (nm)")
plt.title("Stay-in-pocket: ligand–protein minimal distance over 100 ns MD")
plt.xlim(0, 100)
plt.ylim(0, 0.45)
plt.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
plt.tight_layout()
f1 = os.path.join(FIGDIR, "Fig1_lig_dist_trajectory.png")
plt.savefig(f1, dpi=300)
plt.close()
print("wrote", f1, os.path.getsize(f1), "bytes")

# ---------------- Fig 2: MM-PBSA dG bars ----------------
csvp = os.path.join(REPO, "results", "T_md_final_7_2026-10-05.csv")
rows = []
with open(csvp) as fh:
    for r in csv.DictReader(fh):
        rows.append(r)
# order by SYS
order = [c for c, _ in SYS]
rows.sort(key=lambda r: order.index(r["cid"]))
labels = [dict(SYS)[r["cid"]] for r in rows]
gb = [float(r["dg_gb_kcal_mol"]) for r in rows]
gb_sd = [float(r["dg_gb_sd"]) for r in rows]
pb = [float(r["dg_pb_kcal_mol"]) for r in rows]
pb_sd = [float(r["dg_pb_sd"]) for r in rows]

x = np.arange(len(labels))
w = 0.38
plt.figure(figsize=(7.2, 4.6), dpi=300)
b1 = plt.bar(x - w/2, gb, w, yerr=gb_sd, capsize=3, color="#1f77b4", label="GB")
b2 = plt.bar(x + w/2, pb, w, yerr=pb_sd, capsize=3, color="#ff7f0e", label="PB")
plt.axhline(0, color="#444444", lw=0.8)
plt.xticks(x, labels, rotation=30, ha="right", fontsize=8)
plt.xlabel("Complex")
plt.ylabel(r"MM-PBSA $\Delta G$ (kcal/mol, mean $\pm$ SD)")
plt.title("Binding free energy (99 snapshots): GB vs PB")
plt.legend(fontsize=8)
plt.tight_layout()
f2 = os.path.join(FIGDIR, "Fig2_mmpbsa_dg.png")
plt.savefig(f2, dpi=300)
plt.close()
print("wrote", f2, os.path.getsize(f2), "bytes")
print("OK")
