import os
import pandas as pd
import numpy as np

# Paths
base_dir = "."
samn3_path = os.path.join(base_dir, "query results", "SAMN3.csv")
bacdive_clean_path = os.path.join(base_dir, "data", "bacdive_clean.csv")
out_dir = os.path.join(base_dir, "data", "bacdive_only")
trimmed_out_path = os.path.join(out_dir, "SAMN3.csv")
means_csv_path = os.path.join(out_dir, "biosample_bacdive_means.csv")

os.makedirs(out_dir, exist_ok=True)

print(f"Loading BacDive clean taxonomy from {bacdive_clean_path}...")
bacdive_df = pd.read_csv(bacdive_clean_path)
valid_taxa = set(bacdive_df['taxon_id'].dropna().astype(int).unique())

bacdive_subset = bacdive_df[['taxon_id', 'temp_optimum', 'ph_optimum']].drop_duplicates(subset=['taxon_id']).copy()
bacdive_subset['taxon_id'] = bacdive_subset['taxon_id'].astype(int)

print(f"Processing massive file in chunks: {samn3_path}...")

# Remove existing trimmed file if present
if os.path.exists(trimmed_out_path):
    os.remove(trimmed_out_path)

# Dictionary to accumulate totals per biosample across chunks
# biosample -> [sum_temp, count_temp, sum_ph, count_ph]
biosample_accum = {}

chunk_size = 5_000_000
first_chunk = True
total_rows_written = 0

for chunk_idx, chunk in enumerate(pd.read_csv(samn3_path, chunksize=chunk_size, low_memory=False)):
    print(f"  Processing chunk {chunk_idx + 1} ({len(chunk):,} rows)...")
    
    # Clean taxon_id
    chunk['taxon_id'] = pd.to_numeric(chunk['taxon_id'], errors='coerce')
    trimmed_chunk = chunk[chunk['taxon_id'].isin(valid_taxa)].copy()
    
    if len(trimmed_chunk) > 0:
        # Save trimmed chunk to CSV
        trimmed_chunk.to_csv(trimmed_out_path, mode='a', header=first_chunk, index=False)
        first_chunk = False
        total_rows_written += len(trimmed_chunk)
        
        # Merge to get optimum values
        merged = trimmed_chunk.merge(bacdive_subset, on='taxon_id', how='left')
        
        # Group by biosample and update accumulation dict
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

print(f"Finished chunk processing! Total trimmed rows written to SAMN3.csv: {total_rows_written:,}")

# Build dataframe of SAMN3 biosample means
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

samn3_means_df = pd.DataFrame(results)
print(f"Calculated means for {len(samn3_means_df):,} biosamples from SAMN3.csv")

# Merge with existing biosample_bacdive_means.csv
if os.path.exists(means_csv_path):
    print(f"Merging SAMN3 results into existing {means_csv_path}...")
    existing_df = pd.read_csv(means_csv_path)
    combined_df = pd.concat([existing_df, samn3_means_df], ignore_index=True)
    
    # Re-aggregate in case any biosample was split
    final_df = combined_df.groupby('biosample', as_index=False).agg({
        'mean_temp_optimum': 'mean',
        'organisms_used_for_temp': 'sum',
        'mean_ph_optimum': 'mean',
        'organisms_used_for_ph': 'sum'
    })
else:
    final_df = samn3_means_df

final_df.to_csv(means_csv_path, index=False)
print(f"Updated {means_csv_path} successfully with total {len(final_df):,} biosamples.")
