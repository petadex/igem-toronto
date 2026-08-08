# process_bacdive_means.R

library(data.table)

# Define paths
base_dir <- "resources/260629_issue28_automated_metadata"
query_results_dir <- file.path(base_dir, "query results")
bacdive_clean_path <- file.path(base_dir, "data", "bacdive_clean.csv")
out_dir <- file.path(base_dir, "data", "bacdive_only")

# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1. Load BacDive clean data
cat(sprintf("Loading %s...\n", bacdive_clean_path))
bacdive_dt <- fread(bacdive_clean_path)

# Filter out rows with missing taxon_id
bacdive_dt <- bacdive_dt[!is.na(taxon_id)]

# Keep only necessary columns and drop duplicates
bacdive_subset <- unique(bacdive_dt[, .(taxon_id, temp_optimum, ph_optimum)], by = "taxon_id")
valid_taxa <- bacdive_subset$taxon_id

# 2. Process query result CSVs
csv_files <- list.files(query_results_dir, pattern = "\\.csv$", full.names = TRUE)

all_biosample_stats <- list()

for (file_path in csv_files) {
  file_name <- basename(file_path)
  
  if (file_name == "uniqueTaxonIDs.csv") next
  
  cat(sprintf("Processing %s...\n", file_name))
  
  # Read query results
  query_dt <- tryCatch({
    fread(file_path)
  }, error = function(e) {
    cat(sprintf("Error reading %s: %s\n", file_name, e$message))
    return(NULL)
  })
  
  if (is.null(query_dt) || nrow(query_dt) == 0) next
  
  if (!("taxon_id" %in% names(query_dt)) || !("biosample" %in% names(query_dt))) {
    cat(sprintf("Skipping %s - missing required columns.\n", file_name))
    next
  }
  
  # Trim to only include taxon_ids present in bacdive_clean
  trimmed_dt <- query_dt[taxon_id %in% valid_taxa]
  
  # Save the trimmed CSV
  trimmed_out_path <- file.path(out_dir, file_name)
  fwrite(trimmed_dt, trimmed_out_path)
  cat(sprintf("  Saved trimmed CSV to %s (%d rows)\n", trimmed_out_path, nrow(trimmed_dt)))
  
  # Calculate means per biosample
  if (nrow(trimmed_dt) > 0) {
    # Merge with BacDive subset
    merged_dt <- merge(trimmed_dt, bacdive_subset, by = "taxon_id", all.x = TRUE)
    
    # Aggregate stats per biosample
    biosample_stats <- merged_dt[, .(
      mean_temp_optimum = mean(temp_optimum, na.rm = TRUE),
      organisms_used_for_temp = sum(!is.na(temp_optimum)),
      mean_ph_optimum = mean(ph_optimum, na.rm = TRUE),
      organisms_used_for_ph = sum(!is.na(ph_optimum))
    ), by = biosample]
    
    # Clean up NaNs from mean calculations (when all values are NA)
    biosample_stats[is.nan(mean_temp_optimum), mean_temp_optimum := NA]
    biosample_stats[is.nan(mean_ph_optimum), mean_ph_optimum := NA]
    
    all_biosample_stats[[file_name]] <- biosample_stats
  }
}

# 3. Combine and save all biosample stats
if (length(all_biosample_stats) > 0) {
  final_stats_dt <- rbindlist(all_biosample_stats, use.names = TRUE, fill = TRUE)
  
  # It's possible some biosamples were split across query result files (though unlikely). 
  # We can aggregate them again just to be safe.
  final_stats_aggregated <- final_stats_dt[, .(
    mean_temp_optimum = mean(mean_temp_optimum, na.rm = TRUE),
    organisms_used_for_temp = sum(organisms_used_for_temp, na.rm = TRUE),
    mean_ph_optimum = mean(mean_ph_optimum, na.rm = TRUE),
    organisms_used_for_ph = sum(organisms_used_for_ph, na.rm = TRUE)
  ), by = biosample]
  
  final_stats_aggregated[is.nan(mean_temp_optimum), mean_temp_optimum := NA]
  final_stats_aggregated[is.nan(mean_ph_optimum), mean_ph_optimum := NA]
  
  final_csv_path <- file.path(out_dir, "biosample_bacdive_means.csv")
  fwrite(final_stats_aggregated, final_csv_path)
  cat(sprintf("\nSuccessfully compiled biosample means into %s (%d biosamples)\n", final_csv_path, nrow(final_stats_aggregated)))
} else {
  cat("No data processed for biosample means.\n")
}

