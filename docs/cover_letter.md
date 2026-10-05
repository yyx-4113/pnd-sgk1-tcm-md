# Cover Letter

**To the Editor-in-Chief, Computational Biology and Chemistry**

Dear Editor,

Please find enclosed our manuscript entitled *"SGK1 as a computational target for perioperative neurocognitive disorders: in silico traditional Chinese medicine screening and 100 ns molecular dynamics validation"* for consideration as a **Full-length article**.

Perioperative neurocognitive disorders lack molecular targets that are both mechanistically credible and tractable for small-molecule intervention. In prior cross-species work we nominated five hub genes; among them, SGK1 is the only one with an experimentally resolved kinase domain (PDB 7PUE; Halland et al., 2022), which makes it the sole candidate that can be carried from a bioinformatics hypothesis into a structure-based test without modelled structures. We report an in silico screen of traditional Chinese medicine monomers that advanced six compounds (luteolin, resveratrol, kaempferol, daidzein, naringenin, tanshinone IIA) alongside the reference SGK1 inhibitor GSK650394, and we validate all seven complexes by 100 ns molecular dynamics. Every complex retained its ligand inside the ATP-binding pocket (minimal ligand–protein heavy-atom distance 0.14–0.31 nm throughout); global radius-of-gyration and full-complex backbone-RMSD metrics are unreliable for the rebuilt four-chain receptor (same artifact class) and were excluded as stability criteria. MM-PBSA binding free energies placed GSK650394 first (−38.2 ± 14.6 kcal/mol) with the six TCM monomers clustering at −21.1 to −33.3 kcal/mol, driven by van der Waals and hydrophobic terms.

This work fits the Aims & Scope of *Computational Biology and Chemistry*: the journal states that protein modelling and molecular docking studies should be thoroughly validated, and that in the absence of experimental results, molecular dynamics simulations with detailed free energy calculations should be used as complementary techniques to support the major conclusions. Our study is exactly such a complementary-validation design — a purely computational assessment in which 100 ns MD and MM-PBSA serve as the validation layer for docking-derived hypotheses, with all limitations (rebuilt receptor, single-trajectory MM-PBSA without entropy correction, no experimental assay) stated openly.

We have provided the input files, analysis scripts, and output tables in a version-controlled repository (https://github.com/yyx-4113/pnd-sgk1-tcm-md, tag v1.0.1) with a MANIFEST checksum; data are publicly available, not "on request".

A large-language model was used for writing assistance and language polishing; the computational design, simulations, and all reported numerical results were produced by the author. No AI tool performed the calculations.

We confirm that this manuscript is original, not under consideration elsewhere, and that all authors have approved the submission. There is no funding to declare and no conflict of interest.

Thank you for your consideration.

Sincerely,

**Yongxin Yang, B.M.** (ORCID: 0009-0004-9698-6552)
Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China
Email: 960856791@qq.com
