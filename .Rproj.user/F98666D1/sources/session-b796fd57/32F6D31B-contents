##############################################################
####################    Figure 5    ##########################
##############################################################
library(tidyverse)
library(readxl)
library(patchwork)

# Define colors
custom_colors <- c(
  "JRDNN-KM" = "#377EB8",   # Blue
  "SC3" = "red",            # Red
  "BLGGM" = "#4DAF4A",      # Green
  "Seurat" = "orange",      # Orange
  "GENIE3" = "#FFCC00",     # Purple
  "JGNsc" = "#80B7FF",      # Dark Orange
  "JSEM" = "#778899",       # Grey
  "locCSN" = "#8FE506",     # Light Green
  "SpQN" = "#F781BF",       # Pink
  "Normalisr" = "#984EA3",  # Violet
  "CS-CORE" = "#A65628"     # Brown
)

# Define method list
ari_methods <- c("JRDNN", "SC3", "BLGGM", "Seurat")
other_methods <- c("JRDNN", "BLGGM", "GENIE3", "JGNsc", "JSEM", "locCSN","SpQN","Normalisr","CS-CORE")

# Update paths and subfolders
# Assuming 'complex generative mechanisms' is moved to the root 'code_and_data'
base_path <- here("Simulations", "result_data","complex generative mechanisms") 
subfolders <- c("GSD","scMultiSim-T3","SERGIO-DS1")

# Modify data loading function for new folder structure
load_and_process_data <- function(base_path, subfolders) {
  all_data <- list()
  
  for (folder in subfolders) {
    files <- list.files(path = file.path(base_path, folder), pattern = "*.xlsx", full.names = TRUE)
    
    data_list <- lapply(files, function(file) {
      data <- read_excel(file)
      method <- tools::file_path_sans_ext(basename(file))
      data <- data %>% mutate(Method = method)
    })
    
    combined_data <- bind_rows(data_list)
    
    # New grouping method: use folder name as Setting directly
    long_data <- combined_data %>%
      pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value") %>%
      mutate(Setting = folder)
    
    all_data[[folder]] <- long_data
  }
  
  return(bind_rows(all_data))
}

# Load and process data
processed_data <- load_and_process_data(base_path, subfolders)

# Process ARI data
filtered_data_ARI <- processed_data %>%
  filter(Method %in% ari_methods) %>%
  filter(Metric == "ARI") %>%
  mutate(Setting = factor(Setting, levels = subfolders)) %>%
  mutate(Method = case_when(
    Method == "JRDNN" ~ "JRDNN-KM",
    TRUE ~ Method
  )) %>%
  mutate(Method = factor(Method, levels = c("JRDNN-KM", "BLGGM", "SC3", "Seurat")))

# Process F1 data
filtered_data_F1 <- processed_data %>%
  filter(Method %in% other_methods) %>%
  filter(Metric == "F1") %>%
  mutate(Setting = factor(Setting, levels = subfolders)) %>%
  mutate(Method = case_when(
    Method == "JRDNN" ~ "JRDNN-KM",
    TRUE ~ Method
  )) %>%
  mutate(Method = factor(Method, levels = c("JRDNN-KM", "BLGGM",  "JGNsc", "JSEM","SpQN","Normalisr","GENIE3", "locCSN","CS-CORE")))

# Process Recall data
filtered_data_Recall <- processed_data %>%
  filter(Method %in% other_methods) %>%
  filter(Metric == "Recall") %>%
  mutate(Setting = factor(Setting, levels = subfolders)) %>%
  mutate(Method = case_when(
    Method == "JRDNN" ~ "JRDNN-KM",
    TRUE ~ Method
  )) %>%
  mutate(Method = factor(Method, levels = c("JRDNN-KM", "BLGGM",  "JGNsc", "JSEM","SpQN","Normalisr","GENIE3", "locCSN","CS-CORE")))

# Process Precision data
filtered_data_Precision <- processed_data %>%
  filter(Method %in% other_methods) %>%
  filter(Metric == "Precision") %>%
  mutate(Setting = factor(Setting, levels = subfolders)) %>%
  mutate(Method = case_when(
    Method == "JRDNN" ~ "JRDNN-KM",
    TRUE ~ Method
  )) %>%
  mutate(Method = factor(Method, levels = c("JRDNN-KM", "BLGGM",  "JGNsc", "JSEM","SpQN","Normalisr","GENIE3", "locCSN","CS-CORE")))

# Modify plotting function
plot_metric <- function(data, y_label, y_breaks, y_limits, show_legend = FALSE) {
  p <- ggplot(data, aes(x = Setting, y = Value, fill = Method, color = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.size = 0.5) +
    labs(x = "Data Setting", y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
      axis.text.y = element_text(size = 12),                        
      legend.title = element_blank(),
      legend.position = ifelse(show_legend, "left", "none"),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 10)
    ) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits)  
  
  return(p)
}

# Generate plots for each metric
plot_ARI <- plot_metric(
  filtered_data_ARI, 
  y_label = "ARI Value", 
  y_breaks = seq(0.5, 1, 0.1),  
  y_limits = c(0.5, 1), 
  show_legend = TRUE
)

plot_F1 <- plot_metric(
  filtered_data_F1, 
  y_label = "F1 Score", 
  y_breaks = seq(0.25, 0.75, 0.1),  
  y_limits = c(0.25, 0.75)
)

plot_Recall <- plot_metric(
  filtered_data_Recall, 
  y_label = "Recall", 
  y_breaks = seq(0.25, 0.85, 0.1),  
  y_limits = c(0.25, 0.85), 
  show_legend = TRUE
)

plot_Precision <- plot_metric(
  filtered_data_Precision, 
  y_label = "Precision", 
  y_breaks = seq(0.2, 0.7, 0.1),  
  y_limits = c(0.2, 0.7)
)

# Combine plots
combined_plot <- (plot_ARI | plot_F1) / (plot_Recall | plot_Precision) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") &
  theme(legend.position = 'left')

# Display and save plot
print(combined_plot)

