import os
import pandas as pd
import glob

def main():
    input_dir = 'resources/260629_issue28_automated_metadata/query results'
    output_file = 'resources/260629_issue28_automated_metadata/data/unique_taxa.csv'
    
    print(f"Scanning for CSV files in {input_dir}...")
    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    print(f"Found {len(csv_files)} CSV files.")
    
    if len(csv_files) == 0:
        print("No CSV files found. Please ensure you have downloaded the Athena outputs from S3.")
        return
        
    unique_taxa = {}
    total_processed_rows = 0
    
    for i, file_path in enumerate(csv_files):
        if os.path.getsize(file_path) < 100:
            continue
            
        print(f"[{i+1}/{len(csv_files)}] Reading {os.path.basename(file_path)}...")
        try:
            chunk_size = 100000
            for chunk in pd.read_csv(file_path, chunksize=chunk_size):
                total_processed_rows += len(chunk)
                chunk = chunk.dropna(subset=['taxon_id'])
                for _, row in chunk.iterrows():
                    tax_id = int(row['taxon_id'])
                    tax_name = str(row['taxon_name'])
                    unique_taxa[tax_id] = tax_name
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            
    print(f"Processed {total_processed_rows} total rows.")
    print(f"Found {len(unique_taxa)} unique taxonomic groups.")
    
    df_out = pd.DataFrame([
        {'taxon_id': tax_id, 'taxon_name': tax_name}
        for tax_id, tax_name in unique_taxa.items()
    ])
    df_out.to_csv(output_file, index=False)
    print(f"Successfully wrote unique taxonomic groups to {output_file}")

if __name__ == '__main__':
    main()
