import os
import pandas as pd
import numpy as np

# Paths
base_dir = "."
query_dir = os.path.join(base_dir, "query results")
bacdive_clean_path = os.path.join(base_dir, "data", "bacdive_clean.csv")
out_dir = os.path.join(base_dir, "data", "bacdive_only")
means_csv_path = os.path.join(out_dir, "biosample_bacdive_means.csv")

os.makedirs(out_dir, exist_ok=True)

# List of files to process (SAMN14 through SAMN18)
target_files = ["SAMN14.csv", "SAMN15.csv", "SAMN16.csv", "SAMN17.csv", "SAMN18.csv"]

print(f"Loading BacDive clean taxonomy from {bacdive_clean_path}...")
bacdive_df = pd.read_csv(bacdive_clean_path)
valid_taxa = set(bacdive_df['taxon_id'].dropna().astype(int).unique())

bacdive_subset = bacdive_df[['taxon_id', 'temp_optimum', 'ph_optimum']].drop_duplicates(subset=['taxon_id']).copy()
bacdive_subset['taxon_id'] = bacdive_subset['taxon_id'].astype(int)

all_new_means = []

chunk_size = 5_000_000

for filename in target_files:
    file_path = os.path.join(query_dir, filename)
    trimmed_out_path = os.path.join(out_dir, filename)
    
    if not os.path.exists(file_path):
        print(f"File {filename} not found in {query_dir}, skipping...")
        continue

    print(f"\n==========================================")
    print(f"Processing {filename} in chunks...")
    print(f"==========================================")
    
    # Overwrite trimmed file to clean any partial writes (e.g., failed SAMN14)
    if os.path.exists(trimmed_out_path):
        os.remove(trimmed_out_path)

    biosample_accum = {}
    first_chunk = True
    total_rows_written = 0

    for chunk_idx, chunk in enumerate(pd.read_csv(file_path, chunksize=chunk_size, low_memory=False)):
        print(f"  [{filename}] Chunk {chunk_idx + 1} ({len(chunk):,} rows)...")
        
        # Clean taxon_id
        chunk['taxon_id'] = pd.to_numeric(chunk['taxon_id'], errors='coerce')
        trimmed_chunk = chunk[chunk['taxon_id'].isin(valid_taxa)].copy()
        
        if len(trimmed_chunk) > 0:
            # Write trimmed chunk to CSV
            trimmed_chunk.to_csv(trimmed_out_path, mode='a', header=first_chunk, index=False)
            first_chunk = False
            total_rows_written += len(trimmed_chunk)
            
            # Merge with BacDive subset
            merged = trimmed_chunk.merge(bacdive_subset, on='taxon_id', how='left')
            
            # Aggregate totals per biosample across chunks
            for b_id, group in merged.groupby('biosample'):
                temp_vals = group['temp_optimum'].dropna()
                ph_vals = group['ph_optimum'].dropna()
                
                sum_t, cnt_t = temp_vals.sum(), len(temp_vals)
                sum_p, cnt_p = ph_vals.sum(), len(ph_vals)
                
                if b_id not in biosample_accum:
                    biosample_accum[b_id] = [sum_t, cnt_t, sum_p, cnt_p]
                else:
                    biosample_accum[b_id][0] += sum_t
                    biosample_accum[b_id][1] += cnt_t
                    biosample_accum[b_id][2] += sum_p
                    biosample_accum[b_id][3] += cnt_p

    print(f"Finished {filename}! Total trimmed rows saved: {total_rows_written:,}")

    # Build dataframe of biosample means for this file
    results = []
    for b_id, (sum_t, cnt_t, sum_p, cnt_p) in biosample_accum.items():
        mean_t = (sum_t / cnt_t) if cnt_t > 0 else np.nan
        mean_p = (sum_p / cnt_p) if cnt_p > 0 else np.nan
        results.append({
            'biosample': b_id,
            'mean_temp_optimum': mean_t,
            'organisms_used_for_temp': cnt_t,
            'mean_ph_optimum': mean_p,
            'organisms_used_for_ph': cnt_p
        })

    if results:
        file_means_df = pd.DataFrame(results)
        all_new_means.append(file_means_df)
        print(f"Compiled means for {len(file_means_df):,} biosamples from {filename}")

# Update biosample_bacdive_means.csv
if all_new_means:
    print("\nMerging all new biosample means into biosample_bacdive_means.csv...")
    new_means_combined = pd.concat(all_new_means, ignore_index=True)
    
    if os.path.exists(means_csv_path):
        existing_df = pd.read_csv(means_csv_path)
        print(f"Existing biosample_bacdive_means.csv count: {len(existing_df):,} biosamples")
        combined_df = pd.concat([existing_df, new_means_combined], ignore_index=True)
    else:
        combined_df = new_means_combined

    # Re-aggregate across any potential duplicates across files and deduplicate strictly
    final_df = combined_df.groupby('biosample', as_index=False).agg({
        'mean_temp_optimum': 'mean',
        'organisms_used_for_temp': 'sum',
        'mean_ph_optimum': 'mean',
        'organisms_used_for_ph': 'sum'
    })
    
    final_df = final_df.drop_duplicates(subset=['biosample'], keep='last')
    final_df.to_csv(means_csv_path, index=False)
    print(f"Successfully saved updated {means_csv_path} with total {len(final_df):,} unique biosamples.")
