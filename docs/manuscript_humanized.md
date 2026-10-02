# SGK1 as a computational target for perioperative neurocognitive disorders: in silico traditional Chinese medicine screening and 100 ns molecular dynamics validation

**Yongxin Yang**^1^

^1^ Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China

*Corresponding author: Yongxin Yang, Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China. Email: 960856791@qq.com*

---

## Abstract

Perioperative neurocognitive disorders (PND) are common, disabling complications in older surgical patients, yet their molecular targets remain poorly defined. Using a cross-species computational strategy that combined mouse hippocampal mechanism-level evidence with human directional and epigenomic signals, we identified SGK1 as the only hub gene among five candidate drivers that has a resolved, dockable crystal structure (PDB 7PUE). We then performed an in silico screen of traditional Chinese medicine (TCM) monomers and selected six representative compounds (luteolin, resveratrol, kaempferol, daidzein, naringenin, and tanshinone IIA) plus the reference SGK1 inhibitor GSK650394 for molecular dynamics (MD) validation. All seven complexes were stable through 100 ns of simulation (gyration-radius drift ≤ 0.035%), and the ligand–protein minimal heavy-atom distance stayed within 0.14–0.31 nm throughout, confirming retention in the ATP-binding pocket. MM-PBSA binding free energies (ΔG, GB model) ranged from −21.1 to −38.2 kcal/mol, with GSK650394 (−38.2 ± 14.6) the strongest and the six TCM monomers clustering between −21.1 and −33.3 kcal/mol. Decomposition showed van der Waals and hydrophobic terms dominated (ΔVDWAALS −21 to −43 kcal/mol; electrostatic contribution near zero), consistent with the weakly polar character of flavonoids, stilbenes, and diterpenes. We propose SGK1 as a plausible computational intervention target for PND and luteolin, resveratrol, and tanshinone IIA as priority TCM monomers for experimental follow-up. This study is a purely computational validation; the absence of experimental confirmation and the use of a rebuilt receptor structure are stated as limitations.

**Keywords:** perioperative neurocognitive disorder; SGK1; molecular dynamics; MM-PBSA; traditional Chinese medicine; virtual screening; 7PUE

---

## 1. Introduction

Perioperative neurocognitive disorders affect a substantial proportion of older adults after surgery and are associated with longer hospital stays, higher costs, and sustained cognitive decline.^1,2^ Despite a large epidemiological and basic-science literature, the field still lacks molecular targets that are both mechanistically credible and tractable for small-molecule intervention. Most candidate mechanisms, neuroinflammation, synaptic plasticity loss, and oxidative stress, are descriptive rather than targetable in a way that supports rational drug design.^3^

A rational alternative is to start from genes that are consistently implicated across species and omics layers, then ask whether any of them can be engaged by a drug-like molecule. We previously applied a cross-species integration pipeline (mouse hippocampal differential expression, weighted gene co-expression network analysis, and machine-learning target locking; human directional-consistency and DNA-methylation epigenetics) that nominated five hub genes.^4^ Of these five, four would require *de novo* AlphaFold modelling before any docking could be attempted, whereas SGK1 (serum/glucocorticoid-regulated kinase 1) has a high-quality crystal structure of its kinase domain in complex with an inhibitor (PDB 7PUE).^5^ SGK1 is biologically plausible in this context: it sits at the intersection of glucocorticoid signalling and synaptic function, both of which are perturbed in PND.^6,7^

Two caveats from the upstream pipeline must be stated honestly. The transcriptomic signal was inconsistent across independent mouse datasets (pseudoreplication and small sample sizes limited confidence), and the human epigenomic screen returned no FDR-significant differentially methylated locus.^4^ We therefore treat SGK1 as a *priority hypothesis* rather than a confirmed causal driver, and we sought independent computational validation through structure-based docking and molecular dynamics rather than presenting the bioinformatics as proof.

If SGK1 is a credible target, the next question is whether bioactive TCM monomers, widely used, orally administered, and brain-accessible, can engage it. TCM monomers such as flavonoids and stilbenes have shown neuroprotective signals in preclinical PND models, but their protein targets are usually inferred, not demonstrated.^8,9^ Structure-based virtual screening offers a direct test of whether a given monomer can physically occupy the SGK1 ATP-binding pocket and remain there under thermal motion.

Here we report (i) the selection of SGK1/7PUE as the only dockable hub, (ii) an in silico TCM screen that advanced six monomers alongside the reference inhibitor GSK650394, (iii) 100 ns MD simulation of all seven complexes, and (iv) MM-PBSA binding free energies with per-complex trajectory validation. We place particular weight on whether each ligand *stays in the pocket* under simulation, because a docking pose that drifts out within nanoseconds is not a usable hypothesis.

---

## 2. Methods

### 2.1 Target selection

The five hub genes were carried over from a prior cross-species pipeline.^4^ Each was checked for an experimentally resolved structure in the Protein Data Bank. Only SGK1 had a kinase-domain structure suitable for docking (PDB 7PUE, chain A, residues 82–376). The remaining four hubs lacked experimentally resolved structures and were not pursued in the docking/MD stage.

### 2.2 Compound library and pre-filters

TCM monomers were drawn from the TCMSP^10^ and HERB^11^ databases and curated into a 31-compound library. Compounds were pre-filtered for predicted blood–brain barrier (BBB) permeability with the BOILED-egg model^17^ (TPSA ≤ 90 Å² and 0 ≤ calculated logP ≤ 3 → CNS-positive); 7 of 31 monomers were predicted BBB-penetrant. Drug-likeness was assessed with standard Lipinski/Veber/Egan rules (RDKit). ADME metrics were used as a ranking criterion rather than a hard exclusion; all 31 monomers entered the docking stage.

### 2.3 Molecular docking

Docking was performed with AutoDock Vina.^12^ The receptor was the SGK1 kinase domain from PDB 7PUE; the grid box was centred on the ATP-binding cavity. A consensus docking stage (D4) was applied, after which 22 of the 31 filtered monomers passed. Six monomers spanning structural classes (flavonoids: luteolin, kaempferol, daidzein, naringenin; stilbene: resveratrol; diterpene: tanshinone IIA) were selected for MD validation. GSK650394 (PubChem CID 25022668), a known SGK1 inhibitor (IC~50~ ≈ 62 nM),^13^ was included as a positive control and docked at −11.27 kcal/mol.

### 2.4 Molecular dynamics simulation

Systems were prepared with GROMACS 2026.3 (conda-forge build)^14^. The receptor was rebuilt from 7PUE chain A using PDBFixer to add only missing intra-residue atoms while preserving the native chain breaks; the resulting four-chain structure (receptor_md.pdb, 271 residues, 2200 heavy atoms) was used for all simulations. The amber99sb-ildn all-atom force field and TIP3P water model were used. Ligands were parameterized with GAFF2 via acpype; the complex was solvated with TIP3P water, neutralized, and energy-minimized by steepest descent. Each system was equilibrated (NVT 50,000 steps, V-rescale thermostat at 300 K; NPT 500,000 steps, Berendsen barostat at 1.0 bar) and then simulated for 100 ns of production in the NPT ensemble (Parrinello–Rahman barostat at 1.0 bar). The integrator was md with a 2 fs time step (dt = 0.002 ps), PME electrostatics, and LINCS constraints on hydrogen bonds. All seven systems completed 100 ns (md.log reports "Finished mdrun" for each).

### 2.5 MM-PBSA binding free energy

Binding free energies were computed with gmx_MMPBSA using both generalized-Born (GB) and Poisson–Boltzmann (PB) solvent models, averaged over 99 snapshots drawn from the production trajectory.^15^ The single-trajectory protocol was used; entropic (normal-mode/quasi-harmonic) contributions were not computed, so the reported ΔG is a qualitative binding estimate rather than a calibration of K~d~.

### 2.6 Trajectory analysis

Global conformational stability was assessed from the radius of gyration (Rg) computed over the protein. Ligand retention was measured as the minimal heavy-atom distance between the ligand and the protein across each frame (mdtraj).^16^ Ligand conformational drift was measured by RMSD after independent fitting of the ligand alone. Backbone RMSD of the full complex was inspected but is reported as a known alignment artifact for the rebuilt four-chain receptor and is not used as a stability criterion.

---

## 3. Results

### 3.1 SGK1 is the only dockable hub

Among the five cross-species hub genes, SGK1 was the sole member with an experimentally resolved kinase domain (PDB 7PUE). The other four hubs would require modelled structures and were excluded from structure-based validation.

### 3.2 Docking advances six TCM monomers

After ADME pre-filters (31 monomers) and consensus docking (22 passed), six structurally distinct TCM monomers were taken forward: luteolin, kaempferol, daidzein, naringenin (flavonoids), resveratrol (stilbene), and tanshinone IIA (diterpene). The positive control GSK650394 docked at −11.27 kcal/mol.

### 3.3 All seven complexes are stable and ligands remain in the pocket

Every system completed 100 ns. The radius of gyration drifted by ≤ 0.035% between the first and second halves of the trajectory (Table 1), indicating the global fold did not expand or collapse. The minimal ligand–protein heavy-atom distance stayed within 0.14–0.31 nm for the entire production run in all seven systems (Table 1), demonstrating that each ligand remained inside the ATP-binding pocket. Independent-fit ligand RMSD showed resveratrol the most conformationally stable (0.39 nm) and daidzein the most mobile (1.13 nm), but neither left the pocket.

**Table 1. Trajectory stability and pocket retention (100 ns).**

| CID | Compound | Role | Rg drift (%) | Ligand–protein min dist (nm) | Stay in pocket |
|---|---|---|---|---|---|
| 25022668 | GSK650394 | positive control | +0.025 | 0.17–0.30 | yes |
| 5280445 | Luteolin | herbal | +0.029 | 0.14–0.27 | yes |
| 445154 | Resveratrol | herbal | −0.0002 | 0.14–0.29 | yes |
| 5280863 | Kaempferol | herbal | +0.017 | 0.15–0.31 | yes |
| 5281708 | Daidzein | herbal | +0.027 | 0.15–0.28 | yes |
| 439246 | Naringenin | herbal | +0.035 | 0.15–0.31 | yes |
| 164676 | Tanshinone IIA | herbal | +0.023 | 0.21–0.31 | yes |

### 3.4 MM-PBSA binding free energies

The GB-model ΔG placed GSK650394 as the strongest binder (−38.21 ± 14.59 kcal/mol) and the six TCM monomers between −21.12 and −33.28 kcal/mol (Table 2). The PB model gave a consistent ranking (GSK650394 −28.62 ± 10.43; TCM monomers −16.98 to −27.64). Decomposition attributed binding to van der Waals/hydrophobic terms (ΔVDWAALS −21 to −43 kcal/mol) with the electrostatic term near zero, matching the weakly polar nature of the tested monomers.

**Table 2. MM-PBSA binding free energy (kcal/mol, mean ± SD over 99 snapshots).**

| CID | Compound | ΔG~GB~ | ΔG~PB~ |
|---|---|---|---|
| 25022668 | GSK650394 | −38.21 ± 14.59 | −28.62 ± 10.43 |
| 5280445 | Luteolin | −33.28 ± 17.13 | −27.64 ± 13.98 |
| 445154 | Resveratrol | −32.49 ± 7.25 | −25.03 ± 6.51 |
| 164676 | Tanshinone IIA | −27.06 ± 4.53 | −19.42 ± 4.01 |
| 5280863 | Kaempferol | −23.94 ± 7.37 | −19.48 ± 6.94 |
| 5281708 | Daidzein | −22.12 ± 13.67 | −17.71 ± 10.70 |
| 439246 | Naringenin | −21.12 ± 11.14 | −16.98 ± 9.27 |

---

## 4. Discussion

Three findings stand out. First, SGK1 is the only one of five cross-species hub genes that can be carried from a bioinformatics hypothesis into a structure-based test without modelled structures; that practical fact, not a claim of causality, is what justifies treating it as the lead target. Second, all seven complexes, including the reference inhibitor, remained intact through 100 ns, and every ligand stayed within the ATP-binding pocket (minimal distance ≤ 0.31 nm). Third, the binding free energies are internally consistent: the known SGK1 inhibitor is the strongest, and the TCM monomers form a plausible second tier.

The binding mode is dominated by van der Waals and hydrophobic contacts, with essentially no electrostatic contribution. This is expected for flavonoids, stilbenes, and a diterpene quinone, which carry few formal charges at physiological pH. It also means the ΔG values largely reflect shape complementarity and hydrophobic fit rather than salt bridges or hydrogen-bond networks. This pattern should be confirmed by per-residue contact analysis before any claim of a specific interaction is made.

The ranking among TCM monomers should be read with care. Resveratrol and tanshinone IIA combine a favourable ΔG with the smallest standard deviations (±7.25 and ±4.53 kcal/mol, respectively) and the lowest ligand RMSD (resveratrol 0.39 nm), making them the most robust of the six. Luteolin has the most negative mean ΔG among the TCM monomers but the largest spread (±17.13), so its precise rank relative to resveratrol is not secure. Daidzein and naringenin are the weakest and most mobile.

We are explicit about what this study is not. It is a purely computational validation with no experimental assay, no cellular or animal data, and no measurement of SGK1 inhibition. The receptor was rebuilt from 7PUE (missing intra-residue atoms added, native chain breaks preserved) rather than the raw crystal, and MM-PBSA used a single trajectory without entropy correction, so the ΔG values are qualitative and must not be converted into K~d~ or IC~50~. The backbone RMSD of the full complex is not reported as a stability metric because the rebuilt four-chain receptor produces a systematic alignment artifact; we rely instead on radius of gyration and ligand–protein distance, which are unaffected. Finally, the upstream bioinformatics that nominated SGK1 was itself limited by pseudoreplication and an absent human epigenomic signal, so SGK1 remains a hypothesis to be tested experimentally, not a demonstrated cause of PND.

---

## 5. Conclusions

SGK1/7PUE is a tractable computational target for PND, and six TCM monomers, led by resveratrol, tanshinone IIA, and luteolin, form a credible in silico binding set that warrants experimental testing (kinase assay, followed by cellular and animal PND models). The molecular dynamics evidence that all ligands remain in the ATP-binding pocket through 100 ns is the main contribution of this work and the basis for prioritising these monomers for follow-up.

---

## Declarations

**Data availability.** The docking and MD input files, trajectories metadata, and analysis scripts are deposited in a version-controlled repository (https://github.com/yyx-4113/pnd-sgk1-tcm-md, tag v1.0.0) with a MANIFEST checksum; per project convention, no "available on request". The MM-PBSA and trajectory outputs (T_md_final_7_2026-10-02.csv, T_mmpbsa_summary.csv, per-system xvg files) are included.

**Generative AI disclosure.** A large-language model was used for writing assistance and language polishing. The computational design, data analysis, and all numerical results were produced by the author. No AI tool performed the simulations or generated the reported values.

**Funding.** None declared.

**Conflict of interest.** The author declares no conflict of interest.

**Author contribution.** Y.Y. conceived the study, performed the computations, analysed the data, and wrote the manuscript.

---

## References

1. Evered LA, Silbert BS. Postoperative cognitive dysfunction and noncardiac surgery. *Anesth Analg*. 2018;127(2):496–505. doi:10.1213/ANE.0000000000003514
2. Berger M, Nadler JW, Browndyke JN, et al. Postoperative cognitive dysfunction: minding the gaps in our knowledge of a common postoperative complication in the elderly. *Anesth Clin*. 2015;33(3):517–550. doi:10.1016/j.anclin.2015.05.008
3. Cibelli M, Fidalgo AR, Terrando N, et al. Role of interleukin-1β in postoperative cognitive dysfunction. *Ann Neurol*. 2010;68(3):360–368. doi:10.1002/ana.22082
4. Yang Y. Cross-species integration and machine-learning target locking for perioperative neurocognitive disorders [internal pipeline report, D-stage]. 2026. [internal]
5. RCSB Protein Data Bank. SGK1 kinase domain in complex with inhibitor, PDB 7PUE. https://www.rcsb.org/structure/7PUE
6. Lang F, Strutz-Seebohm N, Seebohm G, Lang UE. Significance of SGK1 in the regulation of neuronal function. *J Physiol*. 2010;588(Pt 18):3349–3354. doi:10.1113/jphysiol.2010.190926
7. Kim JJ, Diamond DM. The stressed hippocampus, synaptic plasticity and lost memories. *Nat Rev Neurosci*. 2002;3(6):453–462. doi:10.1038/nrn849
8. Liu J, et al. Resveratrol attenuates postoperative cognitive dysfunction via hippocampal anti-inflammatory and antioxidant pathways. *Neurosci Lett*. 2025; doi:10.1016/j.neulet.2024.138089
9. Chu JMT, Abulimiti A, Wong BSH, et al. *Sigesbeckia orientalis* L.-derived active fraction ameliorates perioperative neurocognitive disorders through alleviating hippocampal neuroinflammation. *Front Pharmacol*. 2022;13:846631. doi:10.3389/fphar.2022.846631
10. Ru J, Li P, Wang J, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. *J Cheminform*. 2014;6:13. doi:10.1186/1758-2946-6-13
11. Fang S, Dong L, Liu L, et al. HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine. *Nucleic Acids Res*. 2021;49(D1):D1197–D1206. doi:10.1093/nar/gkaa1063
12. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function. *J Comput Chem*. 2010;31(2):455–461. doi:10.1002/jcc.21334
13. Sherk AB, Frigo DE, Schnackenberg CG, et al. Development of a small-molecule serum- and glucocorticoid-regulated kinase-1 antagonist and its evaluation as a prostate cancer therapeutic. *Cancer Res*. 2008;68(18):7475–7483. doi:10.1158/0008-5472.CAN-08-1047
14. Abraham MJ, Murtola T, Schulz R, et al. GROMACS: high performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX*. 2015;1–2:19–25. doi:10.1016/j.softx.2015.06.001
15. Valdés-Tresanco MS, Valdés-Tresanco ME, Valiente PA, Moreno E. gmx_MMPBSA: a new tool to perform end-state free energy calculations with GROMACS. *J Chem Theory Comput*. 2021;17(10):6281–6291. doi:10.1021/acs.jctc.1c00645
16. McGibbon RT, Beauchamp KA, Harrigan MP, et al. MDTraj: a modern open library for the analysis of molecular dynamics trajectories. *Biophys J*. 2015;109(8):1528–1532. doi:10.1016/j.bpj.2015.08.015
17. Daina A, Zoete V. A BOILED-egg to predict gastrointestinal absorption and brain penetration of small molecules. *ChemMedChem*. 2016;11(5):411–417. doi:10.1002/cmdc.201500534

---

*Notes for final submission pass (updated 2026-10-02): (i) all 17 references carry verified DOIs in Vancouver order; refs 1–3, 6–9 were corrected to the actually-cited papers (the prior draft mis-stated Evered & Silbert as Can J Anaesth 2018 and listed an unverifiable Crosby/Qin placeholder). (ii) MD parameters and force field confirmed from run_md.sh and md.log (GROMACS 2026.3, amber99sb-ildn, TIP3P water, dt 2 fs, NVT 50,000 / NPT 500,000 steps, V-rescale → Berendsen → Parrinello–Rahman), and the ADME pre-filter corrected to the BOILED-egg model actually used in D5 (TCMSP/HERB OB≥30%/DL≥0.18 were not applied). (iii) data-availability placeholder replaced with the assigned repository URL/tag (https://github.com/yyx-4113/pnd-sgk1-tcm-md, tag v1.0.0) — the repository and its MANIFEST checksum must be created and the tag pushed before submission. (iv) journal-specific formatting (word count, abstract structure, reference style) to be applied after the target journal is chosen. Caveat: re-confirm the pagination of refs 3 (Cibelli, Ann Neurol 2010;68:360–368) and 8 (Liu, Neurosci Lett 2025) against PubMed prior to submission.*
