# plot_bacdive_clean_means.R

options(scipen = 999)

library(ggplot2)
library(data.table)
library(scales)

# Define paths
data_file <- "data/bacdive_clean.csv"
output_dir <- "plots"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Load data
message("Loading data...")
cols_to_load <- c("taxon_id", "organism", "temp_optimum", "temp_min", "temp_max", "ph_optimum", "ph_min", "ph_max")
df <- fread(data_file, select = cols_to_load)
total_organisms <- nrow(df)

# Function to generate and save histogram with MEAN line & label and NA counts
create_plot <- function(data, column_name, title_name, unit, file_name, bin_width) {
    message(paste("Processing", column_name, "..."))
    
    # Calculate NULL / NA counts
    val_vector <- data[[column_name]]
    is_na <- is.na(val_vector) | val_vector == "" | val_vector == "NULL"
    count_na <- sum(is_na)
    count_valid <- total_organisms - count_na
    
    # Filter valid data
    sub_df <- data.frame(val = as.numeric(val_vector[!is_na]))
    
    # Exclude NAs introduced by coercion if any
    sub_df <- sub_df[!is.na(sub_df$val), , drop = FALSE]
    
    # Calculate mean
    mean_val <- mean(sub_df$val, na.rm = TRUE)
    mean_text <- paste0("Mean: ", round(mean_val, 2), " ", unit)
    
    # Subtitle text
    sub_text <- paste0("Total Organisms: ", comma(total_organisms), 
                       " | Valid Results: ", comma(count_valid), 
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
            title = paste("Histogram of", title_name), 
            subtitle = sub_text, 
            x = paste(title_name, "(", unit, ")", sep=""), 
            y = "Frequency"
        )
    
    out_path <- file.path(output_dir, file_name)
    message(paste("Saving plot to", out_path))
    ggsave(out_path, plot = p, width = 9, height = 6, dpi = 300)
}

# Generate plots for temperatures with mean annotation
create_plot(df, "temp_optimum", "Temperature Optimum", "°C", "bacdive_temp_optimum_mean_histogram.png", bin_width = 1)
create_plot(df, "temp_min", "Minimum Temperature", "°C", "bacdive_temp_min_mean_histogram.png", bin_width = 1)
create_plot(df, "temp_max", "Maximum Temperature", "°C", "bacdive_temp_max_mean_histogram.png", bin_width = 1)

# Generate plots for pH with mean annotation
create_plot(df, "ph_optimum", "pH Optimum", "pH", "bacdive_ph_optimum_mean_histogram.png", bin_width = 0.2)
create_plot(df, "ph_min", "Minimum pH", "pH", "bacdive_ph_min_mean_histogram.png", bin_width = 0.2)
create_plot(df, "ph_max", "Maximum pH", "pH", "bacdive_ph_max_mean_histogram.png", bin_width = 0.2)

message("Done generating all BacDive mean plots.")
