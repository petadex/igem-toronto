# plot_biosample_bacdive_means.R

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
df <- fread(data_file)
total_biosamples <- nrow(df)

# Function to generate and save histogram with MEAN line & label
create_plot <- function(data, column_name, organisms_col, title_name, unit, file_name, bin_width) {
    message(paste("Processing", column_name, "..."))
    
    # Calculate NULL / NA counts
    val_vector <- data[[column_name]]
    is_na <- is.na(val_vector)
    count_na <- sum(is_na)
    count_valid <- total_biosamples - count_na
    
    # Filter valid data
    sub_df <- data.frame(val = as.numeric(val_vector[!is_na]))
    
    # Calculate overall mean across all biosamples
    mean_val <- mean(sub_df$val, na.rm = TRUE)
    mean_text <- paste0("Mean of Means: ", round(mean_val, 2), " ", unit)
    
    # Subtitle text
    sub_text <- paste0("Total Biosamples: ", comma(total_biosamples), 
                       " | Valid Biosamples: ", comma(count_valid), 
                       " | Missing/NULL: ", comma(count_na),
                       "\n", mean_text)
    
    p <- ggplot(sub_df, aes(x = val)) + 
        geom_histogram(binwidth = bin_width, fill = "steelblue", color = "black") + 
        geom_vline(xintercept = mean_val, color = "blue", linetype = "dashed", linewidth = 1) +
        annotate(
            "label", 
            x = mean_val, 
            y = Inf, 
            label = mean_text, 
            vjust = 1.5, 
            hjust = -0.1, 
            color = "blue", 
            fontface = "bold", 
            fill = "white",
            size = 4
        ) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
        scale_y_continuous(labels = scales::comma) +
        theme_minimal() + 
        labs(
            title = paste("Histogram of Biosample Mean", title_name), 
            subtitle = sub_text, 
            x = paste("Mean", title_name, "(", unit, ")"), 
            y = "Frequency (Number of Biosamples)"
        )
    
    out_path <- file.path(output_dir, file_name)
    message(paste("Saving plot to", out_path))
    ggsave(out_path, plot = p, width = 9, height = 6, dpi = 300)
}

# Generate plots for the available optimum columns
create_plot(df, "mean_temp_optimum", "organisms_used_for_temp", "Temperature Optimum", "°C", "biosample_mean_temp_optimum_histogram.png", bin_width = 1)
create_plot(df, "mean_ph_optimum", "organisms_used_for_ph", "pH Optimum", "pH", "biosample_mean_ph_optimum_histogram.png", bin_width = 0.2)

message("Done generating plots for biosample means.")

