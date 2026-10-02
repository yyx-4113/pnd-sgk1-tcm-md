# Author verification statement

**Manuscript:** *SGK1 as a computational target for perioperative neurocognitive disorders: in silico
traditional Chinese medicine screening and 100 ns molecular dynamics validation*

**Author:** Yongxin Yang (ORCID 0009-0004-9698-6552)

I, the sole author, confirm the following before submission:

1. **Every numerical result in the manuscript is traceable to a file in this repository.**
   - Table 1 (Rg drift, ligand–protein min distance, stay-in-pocket, independent-fit ligand RMSD)
     derives from `analysis/trajectories/<CID>/rg.xvg`, `lig_dist.xvg`, `lig_rmsd.xvg`, and the
     consolidated `analysis/T_md_final_7_2026-10-02.csv`.
   - Table 2 (MM-PBSA ΔG~GB~ / ΔG~PB~) derives from `analysis/T_mmpbsa_summary.csv` and
     `analysis/trajectories/<CID>/FINAL_RESULTS_MMPBSA.dat`.
   - All seven `md.log` files report `Finished mdrun`, confirming 100 ns completion.

2. **The MD force field, water model, and simulation parameters stated in Methods 2.4 were read
   directly from `scripts/run_md.sh` and the GROMACS `md.log` headers** (GROMACS 2026.3 conda-forge;
   amber99sb-ildn; TIP3P; dt = 2 fs; NVT 50,000 / NPT 500,000 steps; V-rescale → Berendsen →
   Parrinello–Rahman). They were not transcribed from memory.

3. **The ADME pre-filter described in Methods 2.2 matches the pipeline actually used (stage D5):**
   a local RDKit BOILED-egg model (Daina & Zoete 2016) for BBB permeability plus Lipinski/Veber/Egan
   drug-likeness, with ADME used as a ranking criterion rather than a hard exclusion. The previously
   cited TCMSP/HERB OB≥30% / DL≥0.18 rule was *not* applied and has been removed from the manuscript.

4. **The receptor is honestly described.** The simulated receptor (`input/receptor_md.pdb`) is a
   rebuilt structure from PDB 7PUE chain A. An earlier wrongly-rebuilt receptor caused simulation
   failure and the corresponding 7 runs were discarded and re-run on the corrected structure. This is
   disclosed in the manuscript Limitations and in `README.md`.

5. **All 17 references carry verified DOIs** in Vancouver order; refs 1–3, 6–9 were corrected to the
   actually-cited papers (the prior draft mis-stated Evered & Silbert's journal and listed
   unverifiable placeholder entries).

6. **Generative-AI disclosure is accurate.** An LLM assisted with writing/language polishing only.
   Computational design, simulations, and all reported values were produced by the author.

The repository MANIFEST provides sha256 checksums for every deposited file; the tag is **v1.0.0**.

---
*Signed: Yongxin Yang, 2026-10-02.*
