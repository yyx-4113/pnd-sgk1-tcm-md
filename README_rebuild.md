# Required re-simulation on a continuous-chain SGK1 receptor

This note records the single conclusion-level correction identified by the
Round-1 independent review (2026-10-05) and the work that remains before the
manuscript's molecular-dynamics validation can be presented as definitive.

## What the review found

The first-round receptor (`input/receptor_md.pdb`) was rebuilt from PDB 7PUE
chain A with PDBFixer, preserving the three native crystallographic chain
breaks (residues 108-109, 135-148, 246-253) as four separate chains A-D.
Production MD applied no inter-chain restraints, so the fragments drifted
apart. Raw `rg.xvg` / `rmsf.xvg` show radius of gyration swinging ~1.9 -> ~4.9 nm
and 0.000 fraction of residues with RMSF < 0.30 nm in all seven systems: the
"complex" is not a folded protein-ligand system. Consequently the absolute
MM-PBSA dG and the "stays in the pocket" claim are valid only as a
pose-retention screen, not as a binding-energy validation. (The ligand
remained in close receptor contact throughout, which is itself expected for a
docked-start pose and does not by itself prove pocket occupancy.)

## What is done

- Manuscript now discloses the four-chain limitation explicitly and frames the
  MD as a pose-retention screen; the absolute dG is reported as qualitative.
- `rebuild_receptor_contiguous.py` builds a continuous-chain receptor by
  MODELLER loop refinement of the three missing loops.
- `re-dock` + `run_md.sh` are re-runnable against the continuous receptor.

## What the author must still run (requires GPU + MODELLER)

1. Fetch human SGK1 kinase-domain sequence, residues 82-376, from UniProt
   (verify the accession at time of use) into `sgk1_82_376.fasta`.
2. `python rebuild_receptor_contiguous.py` -> `input/receptor_contiguous.pdb`.
3. Re-dock the seven ligands into `receptor_contiguous.pdb` (AutoDock Vina;
   grid centred on the ATP cavity; co-crystal ligand already removed).
4. Re-run `run_md.sh` on the continuous receptor (optionally apply light
   position restraints to the structured core during NVT/NPT).
5. Recompute MM-PBSA over a well-equilibrated 200-snapshot window; recompute
   per-residue contact occupancy and pose clustering; report real convergence
   (RMSF fraction of structured core, Rg time series, block-averaged dG error).
   Independently validate the three modelled loops BEFORE trusting the rebuilt
   receptor: build a reference for the full SGK1 kinase domain with AlphaFold2 or
   ESMFold (or run a short unrestrained equilibration of the continuous receptor),
   and report loop heavy-atom RMSD of the three modelled segments (108-109,
   135-148, 246-253) against that reference. Loops that collapse or diverge from
   the reference indicate a failed rebuild and must be remodelled. Only after this
   loop validation and the continuous-chain re-simulation (steps 2-5) succeed may
   the MD layer be presented as a definitive binding validation.
6. Update Tables 1-2 with the corrected trajectories and bump the repository
   tag.

Until step 5 is completed, the manuscript's MD layer is a pose-retention
screen; the docking shortlist and the positive-control recovery (GSK650394)
remain valid.
