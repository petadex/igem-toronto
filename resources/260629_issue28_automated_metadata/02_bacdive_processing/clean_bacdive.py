import os
import re
import pandas as pd

def parse_range_or_value(val_str):
    """
    Parses a string representing a temperature or pH value/range.
    Returns (optimum, min_val, max_val).
    """
    if pd.isna(val_str) or not str(val_str).strip():
        return None, None, None
        
    val_str = str(val_str).strip()
    
    # Split by semicolon to handle multiple entries
    parts = [p.strip() for p in val_str.split(";")]
    
    single_vals = []
    ranges = []
    
    for p in parts:
        # Regex to match ranges like "10-45" or "4.5-5" or negative-looking ranges
        range_match = re.match(r"^(\d+\.?\d*)\s*-\s*(\d+\.?\d*)$", p)
        if range_match:
            try:
                low = float(range_match.group(1))
                high = float(range_match.group(2))
                ranges.append((low, high))
            except ValueError:
                pass
        else:
            # Check if it's a single value
            single_match = re.match(r"^(\d+\.?\d*)$", p)
            if single_match:
                try:
                    single_vals.append(float(single_match.group(1)))
                except ValueError:
                    pass
                    
    # Determine optimum
    optimum = None
    if single_vals:
        # Take the first single value as optimum (BacDive usually lists optimum first)
        optimum = single_vals[0]
    elif ranges:
        # If only ranges, take the midpoint of the first range
        optimum = (ranges[0][0] + ranges[0][1]) / 2.0
        
    # Determine absolute min and max
    all_nums = single_vals.copy()
    for r in ranges:
        all_nums.extend(r)
        
    min_val = min(all_nums) if all_nums else None
    max_val = max(all_nums) if all_nums else None
    
    return optimum, min_val, max_val

def main():
    input_file = "resources/260629_issue28_automated_metadata/data/bacdive_data.csv"
    output_file = "resources/260629_issue28_automated_metadata/data/bacdive_clean.csv"
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return
        
    print(f"Loading raw BacDive data from {input_file}...")
    df = pd.read_csv(input_file)
    
    clean_rows = []
    for idx, row in df.iterrows():
        taxon_id = row["taxon_id"]
        organism = row["organism"]
        found = row["bacdive_found"]
        
        if found != "Yes":
            continue
            
        temp_opt, temp_min, temp_max = parse_range_or_value(row.get("bacdive_temp_c"))
        ph_opt, ph_min, ph_max = parse_range_or_value(row.get("bacdive_ph"))
        
        clean_rows.append({
            "taxon_id": int(taxon_id),
            "organism": organism,
            "temp_optimum": temp_opt,
            "temp_min": temp_min,
            "temp_max": temp_max,
            "ph_optimum": ph_opt,
            "ph_min": ph_min,
            "ph_max": ph_max,
            "oxygen_tolerance": row.get("bacdive_oxygen"),
            "gram_stain": row.get("bacdive_gram"),
            "morphology": row.get("bacdive_morphology")
        })
        
    clean_df = pd.DataFrame(clean_rows)
    clean_df.to_csv(output_file, index=False)
    
    print(f"Saved {len(clean_df)} cleaned records to {output_file}")
    print(f"Stats:")
    print(f"  Total with Temperature Optimum: {clean_df['temp_optimum'].notna().sum()}")
    print(f"  Total with pH Optimum:          {clean_df['ph_optimum'].notna().sum()}")

if __name__ == "__main__":
    main()
