library(tidyverse)
library(cowplot)

# Load data
data_24 <- read_table(
  "barcode_24_length_distribution.txt",
  col_names = c("Count", "Length")
) %>%
  mutate(barcode_length = "24 bp")

data_30 <- read_table(
  "barcode_30_length_distribution.txt",
  col_names = c("Count", "Length")
) %>%
  mutate(barcode_length = "30 bp")

data <- bind_rows(data_24, data_30) %>%
  # Filter out lengths < 20 and > 32
  # There are very few barcodes with length < 20 or > 32
  filter(Length >= 20, Length <= 32)

# Plot histogram
p <- ggplot(data, aes(x = Length, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_cowplot(12) +
  facet_wrap(~barcode_length, ncol = 2) +
  labs(
    x = "Barcode Length (bp)",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_continuous(
    breaks = seq(20, 32, by = 2)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max(data$Count) * 1.1)
  )

# Save plot
ggsave(
  "barcode_length_distribution.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)
