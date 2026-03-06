library(tidyverse)

read_data <- function(data_path, type) {
  # Convert to long format for plotting
  read_csv(data_path, show_col_types = FALSE) %>%
    pivot_longer(
      cols = c("total_barcodes", "unique_barcodes"), # Removed 'data' argument
      names_to = "Status",
      values_to = "Count"
    ) %>%
    group_by(Status) %>%
    mutate(median_count = median(Count)) %>%
    ungroup() %>%
    filter(case_when(
      Status == "unique_barcodes" ~ Count <= 125,
      Status == "total_barcodes" ~ Count <= 1000,
      TRUE ~ TRUE
    )) %>%
    mutate(type = type)
}

# Load data
df_pre <- read_data("barcode_counts.csv", "Pre-Correction")
df_post <- read_data("barcode_counts_after_correction.csv", "Post-Correction")
data <- bind_rows(
  df_pre,
  df_post
)

# Set the order of the Status for consistent coloring
data$Status <- factor(
  data$Status,
  levels = c("total_barcodes", "unique_barcodes"),
  labels = c("Total Barcodes", "Unique Barcodes")
)

# Set the order of the type for consistent faceting
data$type <- factor(
  data$type,
  levels = c("Pre-Correction", "Post-Correction")
)

# Create histogram
p <- ggplot(data, aes(x = Count)) +
  geom_histogram(
    aes(fill = Status),
    alpha = 0.7,
    position = "identity",
    bins = 100
  ) +
  facet_grid(type ~ Status, scales = 'free') +
  theme_cowplot(12) +
  labs(
    x = "Barcode count",
    y = "Count",
    fill = "Status"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Total Barcodes" = "steelblue",
      "Unique Barcodes" = "salmon"
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0.15)
  ) +
  geom_vline(
    aes(xintercept = median_count),
    color = "black",
    linetype = "dashed"
  ) +
  geom_text(
    aes(x = median_count, label = paste0("Median: ", median_count)),
    y = Inf,
    x = 10,
    vjust = 2,
    hjust = -0.1,
    size = 3.5,
  )

# Save plot
ggsave(
  "06_barcode_counts_histogram.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)
