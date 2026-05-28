import sys
import polars as pl

def build_mapping_dict(master_mapping_gz, source_col=3, target_col=0, target_ids=None):
    print("[*] Scanning master mapping file and filtering columns...")

    source_name = f"source_col_{source_col}"
    target_name = f"target_col_{target_col}"

    q = (
        pl.scan_csv(
            master_mapping_gz,
            separator="\t",
            has_header=False,
            with_column_names=lambda cols: [
                target_name if i == target_col else source_name if i == source_col else f"col_{i}"
                for i in range(len(cols))
            ],
            infer_schema_length=0,
            rechunk=False,
        )
        .select([source_name, target_name])
        .filter(pl.col(source_name).is_not_null() & (pl.col(source_name) != ""))
    )

    print("[*] Exploding and normalizing accessions...")
    q = (
        q.with_columns(pl.col(source_name).str.split(";"))
        .explode(source_name)
        .with_columns(
            pl.col(source_name).str.strip_chars()
            # BUG 1 FIX: strip version suffixes so NP_001234.1 → NP_001234
            .str.replace(r"\.\d+$", "")
            .str.replace(r"_[A-Za-z0-9]+$", "")  # only applies for pdb
        )
        # BUG 2 FIX: drop empty strings produced by empty cells after split
        .filter(pl.col(source_name).is_not_null() & (pl.col(source_name) != ""))
    )

    if target_ids is not None:
        q = q.filter(pl.col(source_name).is_in(list(target_ids)))

    # BUG 3 FIX: prefer reviewed (Swiss-Prot) UniProt ACs over unreviewed (TrEMBL)
    # Swiss-Prot ACs are 6 chars; TrEMBL are 10. Sort so shorter (reviewed) sorts first,
    # then deduplicate keeping first.
    mapping_df = (
        q.sort(
            pl.col(target_name).str.len_chars(),  # reviewed ACs are shorter
            descending=False
        )
        .unique(subset=[source_name], keep="first")
        .collect(engine="streaming")
    )

    print("[*] Building lookup dictionary...")
    lookup_dict = dict(zip(mapping_df[source_name], mapping_df[target_name]))
    print(f"[+] Built {len(lookup_dict):,} unique mappings.")
    return lookup_dict

def read_target_ids(input_csv):
    print(f"[*] Reading your target accessions from {input_csv}...")

    header_labels = {"accession", "id", "from", "entry"}
    target_ids = []
    with open(input_csv, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            token = line.strip().split(',')[0].split('\t')[0]
            if token.lower() in header_labels:
                continue
            target_ids.append(token)

    return target_ids

def convert_accessions(target_ids, lookup_dict, output_tsv):
    print(f"[*] Looping through {len(target_ids):,} accessions to map them locally...")
    
    # Open output file and write your mapping pairs instantly
    with open(output_tsv, 'w') as out_f:
        out_f.write("From\tEntry\n")
        
        mapped_count = 0
        for refseq_id in target_ids:
            # Look up the ID in your O(1) dictionary. Returns None if not found.
            uniprot_ac = lookup_dict.get(refseq_id)
            
            if uniprot_ac:
                out_f.write(f"{refseq_id}\t{uniprot_ac}\n")
                mapped_count += 1
                
    print(f"[+] Finished! Successfully converted and saved {mapped_count:,} matching pairs to {output_tsv}.")

if __name__ == "__main__":
    # Change these filenames to match your exact paths
    MASTER_FILE = "idmapping_selected.tab.gz"
    
    # take the command line argument for the input accessions file
    if len(sys.argv) < 3:
        print("Usage: python local_accession_mapping.py <input_accessions.csv> <output_mappings.tsv>")
        sys.exit(1)
    
    INPUT_IDS = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]

    target_ids = read_target_ids(INPUT_IDS)
    target_set = set(target_ids)

    # Run the pipeline
    if "refseq" in INPUT_IDS.lower():
        print("[*] Detected RefSeq accessions in input. Building RefSeq to UniProt mapping...")
        refseq_to_uniprot = build_mapping_dict(MASTER_FILE, source_col=3, target_col=0, target_ids=target_set)
    elif "genbank" in INPUT_IDS.lower():
        print("[*] Detected GenBank accessions in input. Building GenBank to UniProt mapping...")
        refseq_to_uniprot = build_mapping_dict(MASTER_FILE, source_col=17, target_col=0, target_ids=target_set)
    else:
        print("[!] Warning: Could not detect accession type from input filename. Defaulting to PDB mapping.")
        refseq_to_uniprot = build_mapping_dict(MASTER_FILE, source_col=5, target_col=0, target_ids=target_set)

    convert_accessions(target_ids, refseq_to_uniprot, OUTPUT_FILE)