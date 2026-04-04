#!/usr/bin/env python
# coding: utf-8

# In[169]:


## Core dependencies

get_ipython().system('pip install torch biopython fair-esm')


# In[170]:


##Imports

import torch
import random
from typing import List, Dict

# ESM (ESM-2 local model)
import esm

# Structure parsing
from Bio.PDB import PDBParser


# In[171]:


##Load ESM-2 model + alphabet

device = "cuda" if torch.cuda.is_available() else "cpu"

# Load pretrained ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.eval().to(device)

# Batch converter to handle sequences
batch_converter = alphabet.get_batch_converter()
mask_idx = alphabet.mask_idx  # for masking positions


# In[172]:


##Load PETase sequence

#note starting his-tag not included in PDB so removed from sequence: MGSSHHHHHHSSGLVPR
#note ending n-terminal sequence not included in PDB so removed from sequence: DPAANKARKEAELAAATAEQ

sequence = "GSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLE"  # replace with your PETase sequence

data = [("PETase", sequence)]
labels, strs, tokens = batch_converter(data)
tokens = tokens.to(device)


# In[173]:


##Load structure (PDB)

pdb_file = "replace with your filepath.pdb"  # AlphaFold or crystal structure

parser = PDBParser()
structure = parser.get_structure("PETase", pdb_file)

##PDB structure check (ranges and gaps)

#for checking ranges
for pdb_model in structure:
    for chain in pdb_model:
        print("Chain:", chain.id)
        residues = [res.id[1] for res in chain]
        print("Min:", min(residues), "Max:", max(residues))
        print("Unique ranges:", sorted(set(residues))[:10], "...", sorted(set(residues))[-10:])

residues = sorted([res.id[1] for res in chain])

#for checking gaps
gaps = []
for i in range(len(residues)-1):
    if residues[i+1] - residues[i] > 1:
        gaps.append((residues[i], residues[i+1]))

print("Gaps:", gaps)


# In[174]:


##Set valid PDB range based on structure check results
VALID_PDB_RANGE = (30, 292)

##Determine residue shift value (how residue number shifts based on PDB residue labelling)
shift = VALID_PDB_RANGE[0] - 1


# In[175]:


##Define catalytic residues

catalytic_residues = [160, 206, 237]  # Ser–His–Asp (example)
#catalytic residues based on removal of start and end blanked out sequences
#catalytic_residues = [131, 177, 208]


# In[176]:


##Find interface via catalytic pocket neighbors and define function

def get_interface_from_catalytic_pocket(
    structure,
    catalytic_residues,
    cutoff=8.0,
    target_chain="A"
):
    interface_residues = set()
    catalytic_atoms = []

    # --- Step 1: collect catalytic atoms ---
    for pdb_model in structure:
        for chain in pdb_model:
            if chain.id != target_chain:
                continue

            for res in chain:
                res_id = res.id[1]

                # restrict to valid PETase region
                if not (VALID_PDB_RANGE[0] <= res_id <= VALID_PDB_RANGE[1]):
                    continue

                if res_id in catalytic_residues:
                    for atom in res:
                        catalytic_atoms.append(atom)

    # --- Step 2: find neighboring residues ---
    for pdb_model in structure:
        for chain in pdb_model:
            if chain.id != target_chain:
                continue

            for res in chain:
                res_id = res.id[1]

                # restrict to valid PETase region
                if not (VALID_PDB_RANGE[0] <= res_id <= VALID_PDB_RANGE[1]):
                    continue

                for atom in res:
                    for cat_atom in catalytic_atoms:
                        if atom - cat_atom < cutoff:
                            interface_residues.add(res_id)
                            break

    return sorted(list(interface_residues))


# In[177]:


##Run interface finding function to get catalytic residues

interface_positions = get_interface_from_catalytic_pocket(
    structure,
    catalytic_residues,
    cutoff=8.0
)


# In[178]:


##Remove catalytic residues (protect them)

protected_positions = catalytic_residues

interface_positions = [
    pos for pos in interface_positions
    if VALID_PDB_RANGE[0] <= pos <= VALID_PDB_RANGE[1]
    if pos not in protected_positions
]


# In[179]:


##Masking function

def mask_interface(tokens, interface_positions, mask_prob=0.6):
    masked = tokens.clone()
    masked_positions = []

    for pos in interface_positions: # (ctrl + / for multi-line comments)
        seq_pos = pos - shift - 1        # convert once
        token_pos = seq_pos + 1          # account for <cls>
    #     if pos >= masked.shape[1]: #**masked.shape is problematic since tokenizing uses sequence index and interface_positions uses residue id from PDB
    #         continue #meant to be that if the residue is out of bounds of the token tensor (token tensor sequence length) it will skip that residue for masking (problem is that PDB residue id does not align with sequence index values which screws up what is being masked)

        if random.random() < mask_prob:
            masked[0, token_pos] = mask_idx
            masked_positions.append(seq_pos)

    return masked, masked_positions


# In[180]:


##Amino acid bias (PET-friendly)

def apply_amino_acid_bias(probs, alphabet, bias_strength=1.2):
    preferred_aas = ["F", "W", "Y", "L", "I", "V"]

    for aa in preferred_aas:
        aa_id = alphabet.get_idx(aa)
        probs[aa_id] *= bias_strength

    return probs / probs.sum()


# In[181]:


##Sampling with temperature

def sample_with_temperature(
    logits,
    masked_tokens,
    masked_positions,
    alphabet,
    temperature=1.0,
    top_k=5
):
    new_tokens = masked_tokens.clone()

    for pos in masked_positions:
        scaled_logits = logits[0, pos] / temperature
        probs = torch.softmax(scaled_logits, dim=-1)

        probs = apply_amino_acid_bias(probs, alphabet)

        topk_probs, topk_indices = torch.topk(probs, top_k)
        topk_probs = topk_probs / topk_probs.sum()

        sampled_idx = torch.multinomial(topk_probs, 1)
        new_tokens[0, pos + 1] = topk_indices[sampled_idx]

    return new_tokens


# In[182]:


##Decode tokens - get sequence

def decode_sequence(tokens, alphabet):
    return "".join([
        alphabet.get_tok(t) 
        for t in tokens[0].cpu().numpy()
        if alphabet.get_tok(t) not in ["<cls>", "<eos>", "<pad>"]
    ])


# In[183]:


##Compute log-likelihood

def compute_log_likelihood(model, tokens):
    with torch.no_grad():
        logits = model(tokens)["logits"]

    log_probs = torch.log_softmax(logits, dim=-1)

    ll = 0.0
    for i in range(tokens.shape[1]):
        aa = tokens[0, i]
        ll += log_probs[0, i, aa]

    return ll.item()


# In[184]:


##Generate variants

def generate_variants(
    model,
    alphabet,
    base_tokens,
    interface_positions,
    n_variants=100,
    temperatures=[0.6, 0.8, 1.0, 1.2]
):
    variants = []

    for temp in temperatures:
        n_per_temp = n_variants // len(temperatures)

        for _ in range(n_per_temp):

            masked_tokens, masked_positions = mask_interface(
                base_tokens,
                interface_positions
            )

            if len(masked_positions) == 0:
                continue

            with torch.no_grad():
                logits = model(masked_tokens)["logits"]

            new_tokens = sample_with_temperature(
                logits,
                masked_tokens,
                masked_positions,
                alphabet,
                temperature=temp
            )

            seq = decode_sequence(new_tokens, alphabet)
            score = compute_log_likelihood(model, new_tokens)

            variants.append({
                "sequence": seq,
                "score": score,
                "temperature": temp,
                "mutations": masked_positions
            })

    return sorted(variants, key=lambda x: x["score"], reverse=True)


# In[185]:


##Run pipeline

variants = generate_variants(
    model,
    alphabet,
    tokens,
    interface_positions,
    n_variants=100
)


# In[186]:


##Output top variants

for i, v in enumerate(variants[:10]):
    print(f"\nVariant {i+1}")
    print(f"Score: {v['score']:.2f}")
    print(f"Temperature: {v['temperature']}")
    print(f"Mutations (positions): {v['mutations']}")
    print(f"Sequence: {v['sequence']}")


# In[187]:

##View list of interface residues
print(interface_positions)