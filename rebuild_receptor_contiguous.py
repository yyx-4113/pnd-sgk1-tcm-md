#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
rebuild_receptor_contiguous.py
==============================
Build a SINGLE-CONTINUOUS SGK1 kinase-domain chain from the four-chain rebuilt
receptor (receptor_md.pdb) used in the first MD round, by modelling the three
crystallographically unresolved loops with MODELLER loop refinement.

WHY THIS SCRIPT EXISTS
---------------------
The first-round receptor (receptor_md.pdb) preserves the three native chain
breaks of PDB 7PUE as separate chains A-D (residues 82-107 / 110-134 /
149-245 / 254-376) across gaps 108-109, 135-148 and 246-253. With no
inter-chain restraints in production MD those fragments drift apart, so the
global fold (and therefore the absolute MM-PBSA dG) is not trustworthy. This
script closes the three loops so the rebuilt receptor is one continuous
polypeptide that can be restrained as a single rigid body during equilibration
and production.

REQUIREMENTS (run on a machine WITH GPU / MODELLER, NOT in the writing sandbox)
--------------------------------------------------------------------------------
1. MODELLER >= 10 (free academic licence; set MODELLER_KEY or MODELLER_CLEAN_UP).
2. The full SGK1 kinase-domain sequence, residues 82-376, in FASTA form
   (default: sgk1_82_376.fasta). Fetch human SGK1 from UniProt, verify the
   accession at the time of use, and extract residues 82-376. The three gap
   regions (108-109, 135-148, 246-253) must be present in this sequence.
3. input/receptor_md.pdb (the four-chain structure from the first round).

USAGE
-----
    python rebuild_receptor_contiguous.py \
        --receptor input/receptor_md.pdb \
        --fasta   sgk1_82_376.fasta \
        --out     input/receptor_contiguous.pdb

OUTPUT
------
    input/receptor_contiguous.pdb  (continuous chain A, 82-376, ~294 residues)
After this, re-dock the seven ligands into receptor_contiguous.pdb and re-run
run_md.sh (continuous chain -> inter-chain drift is no longer possible; apply
light position restraints to the structured core during NVT/NPT if desired).
"""

import argparse
import os
import sys


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def write_pir(target_id, target_seq, template_id, template_seq, path):
    """Minimal PIR alignment for MODELLER loop modelling.

    CRITICAL: the template and target records MUST be equal length and aligned
    position-by-position. The template is the full 82-376 SGK1 sequence with
    "-" gap characters at the three crystallographic loop regions that have no
    coordinates (108-109, 135-148, 246-253); the target is the full 82-376
    sequence with no gaps. MODELLER models the target residues that fall in the
    template's gap positions.
    """
    if len(template_seq) != len(target_seq):
        sys.exit("ERROR: PIR template (%d) and target (%d) length differ; "
                 "alignment is invalid." % (len(template_seq), len(target_seq)))
    with open(path, "w") as fh:
        fh.write(">P1;%s\nsequence:%s:82:376:::undefined::::\n" % (template_id, template_id))
        fh.write("%s*\n" % template_seq)
        fh.write(">P1;%s\nsequence:%s:82:376:::undefined::::\n" % (target_id, target_id))
        fh.write("%s*\n" % target_seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", default="input/receptor_md.pdb")
    ap.add_argument("--fasta", default="sgk1_82_376.fasta")
    ap.add_argument("--out", default="input/receptor_contiguous.pdb")
    ap.add_argument("--template-id", default="7PUE")
    ap.add_argument("--target-id", default="SGK1_82_376")
    args = ap.parse_args()

    if not os.path.exists(args.receptor):
        sys.exit("ERROR: receptor %s not found" % args.receptor)
    if not os.path.exists(args.fasta):
        sys.exit("ERROR: FASTA %s not found. Fetch human SGK1 (UniProt, verify "
                 "accession) residues 82-376 into this file." % args.fasta)

    target_seq = read_fasta(args.fasta)
    if not (80 <= len(target_seq) <= 320):
        sys.exit("ERROR: target sequence length %d looks wrong (expect ~294 "
                 "for residues 82-376)." % len(target_seq))

    # Template sequence = the FULL 82-376 target string, with "-" gap characters
    # inserted at the three crystallographic loop regions so that template and
    # target are equal length and MODELLER models those loops. Residues 82-376
    # map to 0-based indices 0..294; the gap index ranges (0-based) are:
    #   108-109 -> 26:28,  135-148 -> 53:67,  246-253 -> 164:172
    # (Note: residue 110 is RESOLVED and must NOT be gapped; residue 135 is the
    #  first residue of the 135-148 gap and IS gapped by design, as is 108-109
    #  and 246-253. The four originally-resolved segments are 82-107, 110-134,
    #  149-245 and 254-376.)
    template_seq = list(target_seq)
    for i in range(26, 28):      # 108-109
        template_seq[i] = "-"
    for i in range(53, 67):      # 135-148
        template_seq[i] = "-"
    for i in range(164, 172):    # 246-253
        template_seq[i] = "-"
    template_seq = "".join(template_seq)
    assert len(template_seq) == len(target_seq), \
        "PIR alignment length mismatch: template %d vs target %d" % (
            len(template_seq), len(target_seq))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    pir = os.path.join(os.path.dirname(os.path.abspath(args.out)),
                       "align_loop.pir")
    write_pir(args.target_id, target_seq, args.template_id, template_seq, pir)

    try:
        from modeller import Environ, log, refine, assess
        from modeller.automodel import loopmodel
    except ImportError:
        sys.exit("ERROR: MODELLER is not installed in this environment. Install "
                 "it (academic licence) on a GPU machine and re-run this script.")

    log.verbose()
    env = Environ()

    # Map the template PDB chains onto the template sequence.
    env.io.atom_files_directory = [os.path.dirname(os.path.abspath(args.receptor))]

    class MyLoop(loopmodel):
        def select_loop_atoms(self):
            # model the three unresolved loops
            return self.select_atoms(
                self.residue_range("108:", "109:"),
                self.residue_range("135:", "148:"),
                self.residue_range("246:", "253:"),
            )

    a = MyLoop(env, alnfile=os.path.basename(pir),
               knowns=args.template_id, sequence=args.target_id,
               loop_assess_methods=assess.DOPE)
    a.starting_model = 1
    a.ending_model = 1
    a.loop.starting_model = 1
    a.loop.ending_model = 5
    a.loop.md_level = refine.fast
    a.make()

    print("Loop modelling done. Inspect the DOPE-assessed model and copy the "
          "best to %s" % args.out)
    print("Then re-dock the seven ligands and re-run run_md.sh on the "
          "continuous receptor.")


if __name__ == "__main__":
    main()
