library(tidyverse)

# Load CSV file
data <- read_csv("barcode_counts.csv")

# Convert to long format for plotting
data_long <- pivot_longer(
  data,
  cols = c("total_barcodes", "unique_barcodes"),
  names_to = "Status",
  values_to = "Count"
) %>%
  group_by(Status) %>%
  mutate(median_count = median(Count)) %>%
  ungroup() %>%
  # remove barcodes with extremely high counts
  # for better visualization
  filter(ifelse(Status == "unique_barcodes", Count <= 125, Count <= 1000))

# Set the order of the Status for consistent coloring
data_long$Status <- factor(
  data_long$Status,
  levels = c("total_barcodes", "unique_barcodes"),
  labels = c("Total Barcodes", "Unique Barcodes")
)

# Create histogram
p <- ggplot(data_long, aes(x = Count)) +
  geom_histogram(
    aes(fill = Status),
    alpha = 0.7,
    position = "identity",
    bins = 100
  ) +
  facet_wrap(~Status, ncol = 2, scales = 'free') +
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
  "05_barcode_counts_histogram.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)
