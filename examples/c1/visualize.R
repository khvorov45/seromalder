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

output <- read_csv("examples/c1/output.csv")

output_long <- output %>%
  pivot_longer(-iteration, names_to = "parameter", values_to = "value")

trace_plots <- output_long %>%
  ggplot(aes(iteration, value)) +
  theme_bw() +
  facet_wrap(~parameter, scales = "free_y") +
  geom_line()

ggsave("examples/c1/trace_plots.pdf", trace_plots)
