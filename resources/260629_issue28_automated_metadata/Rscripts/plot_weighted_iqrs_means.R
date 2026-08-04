# plot_weighted_iqrs_means.R

options(scipen = 999)

library(ggplot2)
library(data.table)

# Define paths
data_file <- "data/Weighted IQRs.csv"
output_dir <- "plots"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Load data for all temperature IQR metrics
message("Loading data...")
cols_to_load <- c("temp_min_weighted_iqr", "temp_optimum_weighted_iqr", "temp_max_weighted_iqr")
df <- fread(data_file, select = cols_to_load)

# Function to generate and save histogram with MEAN line & label
create_temp_histogram <- function(data, column_name, title_name, file_name) {
    message(paste("Processing", column_name, "..."))
    
    # Filter out 0s
    sub_df <- data[get(column_name) > 0, .(val = get(column_name))]
    
    # Calculate mean
    mean_val <- mean(sub_df$val, na.rm = TRUE)
    mean_text <- paste0("Mean: ", round(mean_val, 2), " °C")
    
    p <- ggplot(sub_df, aes(x = val)) + 
        geom_histogram(binwidth = 1, fill = "steelblue", color = "black") + 
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
        scale_x_continuous(
            limits = c(0, 40), 
            breaks = seq(0, 40, by = 5),
            oob = scales::oob_keep
        ) +
        scale_y_continuous(labels = scales::comma) +
        theme_minimal() + 
        labs(
            title = paste("Histogram of", title_name), 
            subtitle = paste0("Excluding values equal to 0 | ", mean_text), 
            x = paste(title_name, "(°C)"), 
            y = "Frequency"
        )
    
    out_path <- file.path(output_dir, file_name)
    message(paste("Saving plot to", out_path))
    ggsave(out_path, plot = p, width = 9, height = 6, dpi = 300)
}

# Generate plots for temp_min, temp_optimum, and temp_max with mean
create_temp_histogram(df, "temp_min_weighted_iqr", "Minimum Temperature Weighted IQR", "temp_min_weighted_iqr_mean_histogram.png")
create_temp_histogram(df, "temp_optimum_weighted_iqr", "Temperature Optimum Weighted IQR", "temp_optimum_weighted_iqr_mean_histogram.png")
create_temp_histogram(df, "temp_max_weighted_iqr", "Maximum Temperature Weighted IQR", "temp_max_weighted_iqr_mean_histogram.png")

message("Done generating all Weighted IQR mean plots.")
