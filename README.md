# pnd-sgk1-tcm-md

**SGK1 as a computational target for perioperative neurocognitive disorders: in silico TCM screening and 100 ns molecular dynamics validation — reproducibility package**

This repository is the version-controlled data-and-code deposit cited in the accompanying manuscript
(`docs/manuscript_humanized.md`). It contains the docking/MD *input files*, the *trajectory metadata*
(every reported number traces to a deposited text artifact), and the *analysis scripts* required to
reproduce the study. Per the author's project convention, all data are publicly available here; there
is no "available on request".

> **Study type.** Purely computational. No experimental assay, cell, or animal data were generated.
> SGK1 (PDB 7PUE) is a *priority hypothesis* nominated by an upstream cross-species pipeline, not a
> demonstrated causal driver of PND. See the manuscript Discussion/Limitations for the full honesty statement.

---

## 1. File map

```
pnd-sgk1-tcm-md/
├── README.md                       # this file
├── LICENSE                         # MIT
├── CITATION.cff                    # citation metadata
├── MANIFEST                        # sha256 checksums of every deposited file (relative path)
├── GITHUB_DEPOSIT_SOP.md           # Chinese step-by-step deposit manual
├── author_verification_statement.md
├── .github/workflows/release.yml   # tag-driven GitHub release stub
├── input/
│   ├── receptor_md.pdb             # rebuilt SGK1(7PUE chain A) receptor used in ALL simulations (4 chains, 271 residues)
│   └── lig/
│       ├── 25022668.mol2           # GSK650394 (positive control, IC50 ~62 nM)
│       ├── 5280445.mol2            # luteolin      (flavonoid)
│       ├── 445154.mol2             # resveratrol   (stilbene)
│       ├── 5280863.mol2            # kaempferol    (flavonoid)
│       ├── 5281708.mol2            # daidzein      (flavonoid)
│       ├── 439246.mol2             # naringenin    (flavonoid)
│       └── 164676.mol2             # tanshinone IIA (diterpene)
├── scripts/
│   ├── rebuild_receptor.py         # PDBFixer-based receptor rebuild (adds missing intra-residue atoms only)
│   ├── run_md.sh                   # full MD pipeline: pdb2gmx -> solvate -> ionize -> EM -> NVT -> NPT -> 100 ns production
│   ├── analyze.sh                  # MM-PBSA trigger + stay-in-pocket / stability summary (fixed 2026-10-02)
│   ├── compute_missing_traj.py     # mdtraj re-computation of rg / independent-fit ligand RMSD / ligand–protein min distance
│   └── compute_convergence.py      # mdtraj re-computation of Rg + RMSF convergence
├── analysis/
│   ├── T_md_final_7_2026-10-02.csv  # master per-system MD summary (Table 1 source)
│   ├── T_mmpbsa_summary.csv         # MM-PBSA ΔG GB/PB per system (Table 2 source)
│   └── trajectories/
│       └── <CID>/                   # 7 systems; each contains:
│           ├── lig_dist.xvg         # ligand–protein minimal heavy-atom distance (STAY-IN-POCKET evidence)
│           ├── lig_rmsd.xvg         # independent-fit ligand RMSD (conformational drift)
│           ├── rmsd.xvg             # full-complex backbone RMSD (ARTIFACT for rebuilt 4-chain receptor — not a stability metric)
│           ├── rg.xvg               # radius of gyration (ARTIFACT for rebuilt 4-chain receptor — deposited for transparency, NOT a stability metric)
│           ├── rmsf.xvg             # per-residue RMSF (elevated by rebuilt-receptor chain ends — not a stability metric alone)
│           ├── rgyrate.xvg          # (present for 5/7 systems) gmx gyrate alternative output
│           ├── hbond.xvg            # (present for 6/7 systems) protein H-bond count
│           ├── sasa.xvg             # (present for 6/7 systems) solvent-accessible surface area
│           ├── FINAL_RESULTS_MMPBSA.dat   # gmx_MMPBSA end-state ΔG (GB & PB)
│           ├── FINAL_DECOMP_MMPBSA.dat    # (present for 5/7 systems) residue-level decomposition
│           ├── bind_energy.csv      # per-snapshot binding energies
│           └── md.log               # proves "Finished mdrun" (100 ns completed)
└── docs/
    ├── manuscript_humanized.md      # submitted manuscript (Vancouver refs, data-availability statement points here)
    └── cover_letter.md
```

> **Why some `.xvg` are missing for 5280445 / 445154.** Those two systems were pulled back from the
> compute instance before `gmx hbond`/`gmx sasa`/`gmx gyrate` extras were written. Their
> `lig_dist` / `lig_rmsd` / `rmsd` / `rg` / `rmsf` were instead recomputed from the local `md.xtc`
> with `compute_missing_traj.py` + `compute_convergence.py` (mdtraj), so the *stay-in-pocket* and
> *convergence* evidence is present for **all 7** systems. The missing `hbond`/`sasa`/`rgyrate` are
> secondary diagnostics and do not affect any manuscript number.

---

## 2. Software versions

| Tool | Version | Role |
|---|---|---|
| GROMACS | 2026.3 (conda-forge build) | MD engine (`mdrun`) |
| acpype | (conda-forge) | ligand topology, GAFF2 force field |
| AutoDock Vina | 1.2.x | consensus docking (D4) |
| gmx_MMPBSA | 1.6.x | MM-PBSA ΔG (GB & PB) |
| MDTraj | 1.10.x | trajectory re-analysis (rg, ligand RMSD, min distance) |
| RDKit | 2026.03.6 | BOILED-egg / Lipinski / Veber / Egan ADME filters (D5) |
| PDBFixer | (conda-forge) | receptor rebuild (missing intra-residue atoms only) |
| Python | 3.13 | analysis scripts |

Force field: **amber99sb-ildn** (protein) + **GAFF2** (ligands); water model **TIP3P**;
electrostatics **PME**; constraints **LINCS** on h-bonds.

---

## 3. Receptor preparation (important honesty note)

The receptor deposited here (`input/receptor_md.pdb`) is a **rebuilt** structure derived from
PDB 7PUE chain A (residues 82–376):

- Built with `scripts/rebuild_receptor.py` (PDBFixer): missing **intra-residue** atoms were added,
  native chain breaks were **preserved** (C(i)–N(i+1) > 2.0 Å → TER + split into chains A/B/C/D),
  and original residue numbering was restored. Result: 4 chains / 271 residues / ~2200 heavy atoms.
- This is NOT the raw crystal. The raw 7PUE chain A was used as the starting point; the rebuilt
  version is what was actually simulated.

A prior attempt used a wrongly-rebuilt receptor (`addMissingAtoms()` had fabricated a missing loop
and displaced a segment by ~22.5 Å), which caused an NVT explosion. **That receptor and the 7 MD
runs built on it were discarded and re-run on the clean rebuilt structure.** The runs archived in
`analysis/trajectories/` are the corrected ones (each `md.log` reports "Finished mdrun").

---

## 4. MD protocol (from `scripts/run_md.sh`)

| Stage | Integrator | Steps | Notes |
|---|---|---|---|
| Energy minimization | steepest descent | 50,000 | — |
| NVT equilibration | md | 50,000 | V-rescale thermostat, 300 K, tc-grps=System |
| NPT equilibration | md | 500,000 | Berendsen barostat, 1.0 bar, isotropic |
| Production | md | 50,000,000 | 100 ns; Parrinello–Rahman barostat, 1.0 bar; dt = 0.002 ps (2 fs) |

Global: PME electrostatics; `rcoulomb = rvdw = 1.0 nm`; `nstlist = 10`; LINCS on h-bonds;
genion `-neutral -conc 0.15` (0.15 M NaCl). All seven systems completed the full 100 ns.

---

## 5. MM-PBSA

`gmx_MMPBSA` with both GB and PB solvent models, averaged over 99 snapshots from the production
trajectory (single-trajectory protocol). Entropy (normal-mode / quasi-harmonic) was **not** computed,
so the reported ΔG is a *qualitative* binding estimate and must **not** be converted to K_d or IC_50.

---

## 6. Reproduction

```bash
# 1. (optional) rebuild receptor from 7PUE
python scripts/rebuild_receptor.py 7PUE.pdb input/receptor_md.pdb

# 2. run MD for one system (ligand CID, e.g. 5280445)
bash scripts/run_md.sh 5280445 input/lig/5280445.mol2

# 3. MM-PBSA + stay-in-pocket / stability summary
bash scripts/analyze.sh 5280445

# 4. (if xvg from gmx were not written) recompute from local md.xtc with mdtraj
python scripts/compute_missing_traj.py --all
python scripts/compute_convergence.py --all
```

Raw 100-ns `.xtc` trajectories (~GB each) are **not** deposited on GitHub (size limit). They are
regenerable from the deposited inputs + `run_md.sh`, and every number in the manuscript is recoverable
from the deposited text artifacts (`.xvg`, `.dat`, `.csv`, `md.log`).

---

## 7. Limitations (carry into any reuse)

1. Purely computational; no experimental validation of SGK1 inhibition or binding.
2. Receptor is rebuilt from 7PUE, not the raw crystal.
3. Single-trajectory MM-PBSA without entropy correction → qualitative ΔG.
4. Full-complex backbone RMSD, per-residue RMSF, and radius of gyration (Rg) all contain a systematic
   artifact from the rebuilt four-chain receptor (the separated chains drift apart during the run,
   inflating Rg and backbone RMSD by tens of percent). Stability is judged from the ligand–protein
   minimal distance only; Rg is excluded as a stability metric.
5. Upstream bioinformatics that nominated SGK1 was limited by pseudoreplication and an absent human
   epigenomic signal; SGK1 remains a hypothesis to be tested experimentally.

---

## 8. MANIFEST

`MANIFEST` lists `sha256  <relative-path>` for every deposited file. Verify with:

```bash
sha256sum -c MANIFEST      # GNU/Linux
# or (Git Bash / Windows):
while read h p; do printf '%s  %s\n' "$h" "$p"; done < MANIFEST | sha256sum -c
```

Tag: **v1.0.0**.
