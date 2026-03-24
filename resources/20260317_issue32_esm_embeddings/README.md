# Protein Embeddings Extraction with Pretrained ESM Models

This repository contains scripts for extracting protein sequence embeddings using pretrained ESM (Evolutionary Scale Modeling) models. The scripts provide functionality to process FASTA files and save the resulting embeddings for downstream applications.

## Overview

The repository contains three main modules:

1. **extract_embeddings**  
   - **Purpose:** Extracts per-residue embeddings from a FASTA file containing one or more protein sequences.  
   - **Usage:** Suitable for batch processing and saving the embeddings in individual files (one per sequence).  
   - **Output:** Saves a PyTorch checkpoint (`.pt` file) for each sequence in the specified output directory.

2. **single_protein_embeddings**  
   - **Purpose:** Extracts per-residue embeddings from a FASTA file containing a single protein sequence.  
   - **Usage:** Designed for cases where only one protein is present. Returns a tensor of shape `(L, D)`, where `L` is the number of residues and `D` is the embedding dimension.  
   - **Output:** Prints the shape of the final embedding tensor.

3. **multiple_protein_embeddings**  
   - **Purpose:** Extracts embeddings from a FASTA file containing multiple protein sequences and saves both the sequence names and corresponding embeddings as NumPy arrays.  
   - **Usage:** Useful for downstream analysis when you need to work with all extracted embeddings together.  
   - **Output:** Creates a folder named `npy_embeddings` (if it does not exist) and saves two files inside it: `names.npy` (sequence identifiers) and `embeddings.npy` (embedding arrays).

## Dependencies

Before running these scripts, ensure you have the following installed:

- Python 3.7+
- [PyTorch](https://pytorch.org/) (with CUDA support if running on GPU)
- [ESM](https://github.com/facebookresearch/esm) package (for loading pretrained models and processing FASTA files)
- NumPy
- pathlib (comes with Python 3.4+)

## How to Use

### 1. Extract Embeddings (Batch Processing)

The `extract_embeddings` function processes a FASTA file with potentially multiple protein sequences. To run this:

- **Edit the `main()` function in the script** to update the following paths as needed:
  - `fasta_file`: Path to your input FASTA file.
  - `output_dir`: Path to the directory where you want to save the PyTorch checkpoint files.

### 2. Extract Single Protein Embeddings

For a FASTA file that contains a single protein sequence, use the `extract_per_residue_embeddings` function. Update the file path in the `main()` function and then run

The script will print the shape of the extracted embedding tensor.

### 3. Extract Multiple Protein Embeddings and Save as NumPy Files

This module extracts sequence names and embeddings from a multi-sequence FASTA file and saves them as NumPy arrays. The `main()` function sets up an output directory (`npy_embeddings`) and saves two files:
- `names.npy`: Contains sequence identifiers.
- `embeddings.npy`: Contains embedding arrays (with dtype `object` to handle variable shapes).

To run:

- Update the `fasta_file` path in the script’s `main()` function.

## Customization

- **Model Choice:**  
  The scripts use pretrained ESM models (e.g., `esm2_t33_650M_UR50D` or `esm2_t6_8M_UR50D`). You can change the `model_name` parameter in the functions to switch between models.

- **Sequence Length and Layers:**  
  Adjust the `seq_length` and `repr_layers` parameters to control the maximum sequence length (including special tokens) and the model layers from which to extract representations.

- **Output Directories:**  
  Modify the file paths in the `main()` functions to suit your file organization preferences.
