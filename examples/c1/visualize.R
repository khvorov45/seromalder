library(tidyverse)

input <- read_csv("examples/c1/input.csv")

input_extra <- input %>%
  mutate(titre = 2^log2titre)

input_plot <- input_extra %>%
  ggplot(aes(time, titre)) +
  theme_bw() +
  scale_y_log10("Titre", breaks = 5 * 2^(0:15)) +
  geom_point(alpha = 0.3, shape = 18)

ggsave("examples/c1/input_plot.pdf", input_plot)
