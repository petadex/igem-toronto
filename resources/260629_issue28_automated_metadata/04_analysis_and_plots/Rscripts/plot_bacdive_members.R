options(scipen = 999)

library(ggplot2)
library(data.table)
library(scales)

# Define paths
data_file <- "resources/260629_issue28_automated_metadata/data/bacdive_only/biosample_bacdive_means.csv"
output_dir <- "resources/260629_issue28_automated_metadata/plots"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Load data
message("Loading data...")
cols_to_load <- c("organisms_used_for_temp", "organisms_used_for_ph")
df <- fread(data_file, select = cols_to_load)
total_biosamples <- nrow(df)

# Function to generate and save log1p-scaled histogram including 0 / Missing values
create_plot <- function(data, column_name, title_name, unit, file_name, bins = 50) {
    message(paste("Processing", column_name, "..."))
    
    val_vector <- data[[column_name]]
    
    # Identify missing/NA/NULL/empty and convert to 0
    is_na_or_null <- is.na(val_vector) | val_vector == "" | val_vector == "NULL"
    clean_vals <- ifelse(is_na_or_null, 0, as.numeric(val_vector))
    clean_vals[is.na(clean_vals)] <- 0
    
    count_zero_missing <- sum(clean_vals == 0)
    count_valid_pos <- sum(clean_vals > 0)
    
    sub_df <- data.frame(
        val = clean_vals,
        category = ifelse(clean_vals == 0, "Zero / Missing", "Positive Count")
    )
    
    # Calculate statistics for positive values (>0)
    pos_vals <- clean_vals[clean_vals > 0]
    med_val <- median(pos_vals)
    mean_val <- mean(pos_vals)
    sd_val <- sd(pos_vals)
    
    med_text <- paste0("Median (>0): ", round(med_val, 2))
    mean_text <- paste0("Mean (>0): ", round(mean_val, 2))
    
    # Subtitle text with summary stats
    sub_text <- paste0("Total Biosamples: ", comma(total_biosamples), 
                       " | Positive (>0): ", comma(count_valid_pos), 
                       " | Zero/Missing: ", comma(count_zero_missing),
                       "\nMedian (>0): ", comma(round(med_val, 2)), 
                       " | Mean (>0): ", comma(round(mean_val, 2)), 
                       " | SD (>0): ", comma(round(sd_val, 2)))
    
    p <- ggplot(sub_df, aes(x = val, fill = category)) + 
        geom_histogram(bins = bins, color = "black", alpha = 0.85) + 
        geom_vline(xintercept = med_val, color = "red", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_val, color = "darkblue", linetype = "dotted", linewidth = 1) +
        annotate(
            "label", 
            x = med_val, 
            y = Inf, 
            label = med_text, 
            vjust = 1.5, 
            hjust = -0.1, 
            color = "red", 
            fontface = "bold", 
            fill = "white",
            size = 3.5
        ) +
        annotate(
            "label", 
            x = mean_val, 
            y = Inf, 
            label = mean_text, 
            vjust = 3.5, 
            hjust = -0.1, 
            color = "darkblue", 
            fontface = "bold", 
            fill = "white",
            size = 3.5
        ) +
        scale_fill_manual(
            values = c("Zero / Missing" = "coral2", "Positive Count" = "steelblue"),
            name = ""
        ) +
        scale_x_continuous(
            trans = "log1p",
            breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
            labels = c("0 / Missing", "1", "10", "100", "1,000", "10,000", "100,000")
        ) +
        scale_y_continuous(labels = scales::comma) +
        theme_minimal() + 
        theme(legend.position = "top") +
        labs(
            title = paste("Histogram of", title_name, "(log1p Scale, Including 0/Missing)"), 
            subtitle = sub_text, 
            x = paste0(title_name, " (log1p scale, ", unit, ")"), 
            y = "Frequency (Number of Biosamples)"
        )
    
    out_path <- file.path(output_dir, file_name)
    message(paste("Saving plot to", out_path))
    ggsave(out_path, plot = p, width = 9, height = 6, dpi = 300)
    
    alt_file_name <- paste0("bacdive_", column_name, "_histogram.png")
    if (alt_file_name != file_name) {
        alt_out_path <- file.path(output_dir, alt_file_name)
        message(paste("Saving copy to", alt_out_path))
        ggsave(alt_out_path, plot = p, width = 9, height = 6, dpi = 300)
    }
}

# Generate plot for organisms_used_for_temp
create_plot(
    df, 
    "organisms_used_for_temp", 
    "Organisms Used for Temperature", 
    "organisms", 
    "bacdive_temperature_histogram.png", 
    bins = 50
)

# Generate plot for organisms_used_for_ph
create_plot(
    df, 
    "organisms_used_for_ph", 
    "Organisms Used for pH", 
    "organisms", 
    "bacdive_ph_histogram.png", 
    bins = 50
)

message("Done generating BacDive plots.")
