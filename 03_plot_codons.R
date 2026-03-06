library(tidyverse)
library(cowplot)

# Read the CSV file into R
df <- read.csv(
  "/mnt/4TB_SSD/scripts/gps-library-ont-sequencing/stop_codon_counts.csv"
)

# Reshape the data frame for plotting
df_long <- pivot_longer(
  df,
  cols = c("TAA", "TAG", "TGA", "Other"),
  names_to = "StopCodon",
  values_to = "Count"
)

# Set the order of the stop codons for consistent coloring
df_long$StopCodon <- factor(
  df_long$StopCodon,
  levels = c("TAA", "TGA", "TAG", "Other")
)

# Convert barcode_length to a factor for better control over the x-axis
df_long$barcode_length <- factor(df_long$barcode_length)

# Create the stacked bar plot
p <- ggplot(df_long, aes(x = barcode_length, y = Count, fill = StopCodon)) +
  geom_bar(stat = "identity") +
  theme_cowplot(12) +
  labs(
    x = "Barcode Length",
    y = "Count",
    fill = "Stop Codon"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = c(
      "TAA" = "steelblue",
      "TAG" = "salmon",
      "TGA" = "lightgreen",
      "Other" = "gray"
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
  )

# Save the plot as a PNG file
ggsave(
  "/mnt/4TB_SSD/scripts/gps-library-ont-sequencing/stop_codon_counts.png",
  plot = p,
  width = 3,
  height = 4,
  dpi = 300,
  bg = "white"
)
