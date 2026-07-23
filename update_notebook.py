import json

notebook_path = 'c:/Users/denis/OneDrive - University of Toronto/Documents/iGEM/igem-repo/igem-toronto/notebooks/260629_issue28_automated_metadata.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

r_code = """%%R -w 900 -h 1200
options(scipen = 999)
library(ggplot2)
library(data.table)

# Load data for all temperature IQR metrics
file_path <- "../resources/260629_issue28_automated_metadata/data/Weighted IQRs.csv"
cols_to_load <- c("temp_min_weighted_iqr", "temp_optimum_weighted_iqr", "temp_max_weighted_iqr")
df <- fread(file_path, select = cols_to_load)

# Function to render histogram with median
create_plot <- function(data, column_name, title_name) {
    sub_df <- data[get(column_name) > 0, .(val = get(column_name))]
    med_val <- median(sub_df$val, na.rm = TRUE)
    med_text <- paste0("Median: ", round(med_val, 2), " °C")
    
    ggplot(sub_df, aes(x = val)) + 
        geom_histogram(binwidth = 1, fill = "steelblue", color = "black") + 
        geom_vline(xintercept = med_val, color = "red", linetype = "dashed", linewidth = 1) +
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
            size = 4
        ) +
        scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5), oob = scales::oob_keep) +
        scale_y_continuous(labels = scales::comma) +
        theme_minimal() + 
        labs(
            title = paste("Histogram of", title_name), 
            subtitle = paste0("Excluding values equal to 0 | ", med_text), 
            x = paste(title_name, "(°C)"), 
            y = "Frequency"
        )
}

p1 <- create_plot(df, "temp_min_weighted_iqr", "Minimum Temperature Weighted IQR")
p2 <- create_plot(df, "temp_optimum_weighted_iqr", "Temperature Optimum Weighted IQR")
p3 <- create_plot(df, "temp_max_weighted_iqr", "Maximum Temperature Weighted IQR")

print(p1)
print(p2)
print(p3)
"""

for cell in nb['cells']:
    if cell.get('id') == 'r-plot':
        cell['source'] = [line + '\n' for line in r_code.split('\n')]
        cell['source'][-1] = cell['source'][-1].rstrip('\n')

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)
