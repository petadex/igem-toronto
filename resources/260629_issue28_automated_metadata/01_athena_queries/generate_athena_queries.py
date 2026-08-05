import os
import pandas as pd

def main():
    input_file = 'resources/260629_issue28_automated_metadata/data/petadex_biosamples.csv'
    output_dir = 'notebooks/metadata-automated-bacdrive/queries'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
        
    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file)
    
    # Extract 6-character prefixes (e.g., SAMN31, SAMEA1)
    df['prefix'] = df['biosample'].str[:6]
    prefixes = df['prefix'].dropna().unique()
    print(f"Found {len(prefixes)} unique prefixes.")
    
    for prefix in sorted(prefixes):
        query = f"""-- Query for BioSample prefix: {prefix}
SELECT 
    meta.biosample AS biosample,
    meta.acc AS run_id,
    tax.tax_id AS taxon_id,
    tax.name AS taxon_name,
    tax.self_count AS self_count,
    tax.total_count AS total_count
FROM 
    sra.metadata meta
JOIN 
    sra_tax_analysis_tool.tax_analysis tax ON meta.acc = tax.acc
JOIN
    default.petadex_biosamples pb ON meta.biosample = pb.biosample
WHERE 
    meta.biosample LIKE '{prefix}%';
"""
        out_path = os.path.join(output_dir, f"{prefix}.sql")
        with open(out_path, 'w') as f:
            f.write(query)
            
    print(f"Successfully generated {len(prefixes)} SQL queries in {output_dir}")

if __name__ == '__main__':
    main()
