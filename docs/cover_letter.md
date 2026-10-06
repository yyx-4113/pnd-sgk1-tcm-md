# Cover Letter

**To the Editor-in-Chief, Computational Biology and Chemistry**

Dear Editor,

Please find enclosed our manuscript entitled *"SGK1 as a computational candidate target for perioperative neurocognitive disorders: in silico traditional Chinese medicine screening and 100 ns molecular dynamics pose-retention assessment"* for consideration as a **Full-length article**.

Perioperative neurocognitive disorders lack molecular targets that are both mechanistically credible and tractable for small-molecule intervention. In prior cross-species work we nominated five hub genes; among them, SGK1 is the only one with an experimentally resolved kinase domain (PDB 7PUE; Halland et al., 2022), which makes it the sole candidate that can be carried from a bioinformatics hypothesis into a structure-based test without modelled structures. We report an in silico screen of traditional Chinese medicine monomers that advanced six compounds (luteolin, resveratrol, kaempferol, daidzein, naringenin, tanshinone IIA) alongside the reference SGK1 inhibitor GSK650394, and we assess pose retention for all seven complexes by 100 ns molecular dynamics. Every complex retained the ligand in close receptor contact (minimal ligand–protein heavy-atom distance 0.14–0.31 nm throughout; per-frame maxima ≤0.314 nm). We disclose transparently that the rebuilt receptor preserves three native crystallographic chain breaks as separate chains and that, with no inter-chain restraints in production MD, those fragments can separate; the global radius-of-gyration and full-complex backbone-RMSD metrics are therefore invalid as stability criteria and are excluded, and the absolute MM-PBSA ΔG is reported as qualitative. A continuous-chain receptor rebuild (three missing loops modelled by MODELLER) and the re-simulation that would make the stability claim definitive are provided as a ready-to-run pipeline in the repository. MM-PBSA binding free energies placed GSK650394 first (−38.2 ± 14.6 kcal/mol) with the six TCM monomers clustering at −21.1 to −33.3 kcal/mol, driven by van der Waals and hydrophobic terms. The results are presented as two tables and two figures (ligand–protein minimal-distance trajectories and the MM-PBSA free-energy comparison).

This work fits the Aims & Scope of *Computational Biology and Chemistry*: the journal states that protein modelling and molecular docking studies should be thoroughly validated, and that in the absence of experimental results, molecular dynamics simulations with detailed free energy calculations should be used as complementary techniques to support the major conclusions. Our study is exactly such a complementary pose-retention design — a purely computational assessment in which 100 ns MD and MM-PBSA serve as the pose-retention layer for docking-derived hypotheses, with all limitations (rebuilt four-chain receptor and the required continuous-chain re-simulation, single-trajectory MM-PBSA without entropy correction, no experimental assay) stated openly.

We have provided the input files, analysis scripts, the rebuilt four-chain receptor, the continuous-chain receptor rebuild script, the re-dock/re-MD pipeline, and output tables in a version-controlled repository (https://github.com/yyx-4113/pnd-sgk1-tcm-md, tag v1.1.2) with a MANIFEST checksum; data are publicly available, not "on request". The upstream five-hub nomination is deposited as docs/hub_nomination_report.md.

A large-language model was used for writing assistance and language polishing; the computational design, simulations, and all reported numerical results were produced by the author. No AI tool performed the calculations.

We confirm that this manuscript is original, not under consideration elsewhere, and that all authors have approved the submission. There is no funding to declare and no conflict of interest.

Thank you for your consideration.

Sincerely,

**Yongxin Yang, B.M.** (ORCID: 0009-0004-9698-6552)
Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China
Email: 960856791@qq.com
