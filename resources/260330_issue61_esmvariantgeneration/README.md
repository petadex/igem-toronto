# ESM‑2 Mutagenesis Prototype

Prototype pipeline for generating PETase sequence variants using ESM‑2 with site‑specific masking, temperature‑scaled sampling, and model‑based scoring.

---

## Files

| File | Purpose |
|---|---|
| `esm2_mutationvariant_prototype.py` | Runs the full ESM‑2 mutagenesis workflow: loads PETase sequence + structure, identifies catalytic‑pocket interface residues, masks positions, samples variants at multiple temperatures, and outputs ranked sequences. |
| `notes.md` | Additional notes on PETase residue numbering, PDB alignment, masking logic, and test variants generated |
| `PETase.pdb` | PDB file used for generation of the ten test variants |

---

## Prerequisites

- Python 3.10+
- PyTorch (GPU recommended but optional)
- `fair-esm` for loading the ESM‑2 model
- `biopython` for PDB structure parsing
- A PETase PDB structure file (AlphaFold or crystal structure)
- PETase sequence aligned to the PDB residue numbering (sequence with N and C terminal domains removed to align with the PDB file for generating the ten test variants: `GSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLE`)

Install dependencies (`pip install torch biopython fair-esm`)

## Step 1 - Prepare PETase inputs
- A PETase PDB file
- A PETase amino‑acid sequence with N‑ and C‑terminal tags removed to match PDB numbering

NOTE: If your PDB file is stored elsewhere, update the path inside `esm2_mutationvariant_prototype.py`

## Step 2 - Run the ESM‑2 mutagenesis script (esm2_mutationvariant_prototype.py)

This script performs the following:

Loads the ESM‑2 model (esm2_t33_650M_UR50D)
1) Parses the PETase structure using Biopython
2) Identifies interface residues within 8 Å of the catalytic triad
3) Masks interface residues with probability 0.6
4) Samples replacements at temperatures 0.6, 0.8, 1.0, and 1.2
5) Applies a hydrophobic amino‑acid bias (F/W/Y/L/I/V)
6) Computes ESM‑2 log‑likelihoods for each generated variant
7) Ranks variants and prints the top 10 sequences
   Example output:
   Variant 1
   Score: `-123.45`
   Temperature: `0.8`
   Mutations: `[45, 112, 178]`
   Sequence: `GSHMRGP...`

## NOTE: These variants can be used for downstream analyses such as:

- ESMFold structure prediction
- PET docking simulations
- Stability or ΔΔG prediction
- Experimental screening

