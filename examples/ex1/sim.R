# Single sample per individual, single virus, single infection.
# Everyone infected, all infection times known.
# Only 4 parameters in the model.

library(tidyverse)

# NOTE(sen) Log2-difference induced by infection. Long-term because it doesn't
# wane.
long_term_boost <- 4

# NOTE(sen) Log2-difference induced by infection. Short-term because it wanes
# over time. Adds to long-term boost
short_term_boost <- 2

# NOTE(sen) Time units to titre peak (from infection)
time_to_peak <- 14

# NOTE(sen) Time units for the short-term titre to wane completely (from peak)
time_to_wane <- 50

n_individuals <- 1000

censor_titres <- function(titres) {
  cuts <- cut(titres, c(-Inf, 5 * 2^(1:10), Inf)) %>% as.integer()
  5 * 2^(cuts - 1)
}

# NOTE(sen) Time is units from reference date
sim_data <- tibble(
  pid = 1:n_individuals,
  virus = 1,
  time_infection = as.integer(runif(n_individuals, 0, 75)),
  time_sample = as.integer(runif(n_individuals, 35, 100)),
  time_peak = time_infection + time_to_peak,
  time_wane = time_peak + time_to_wane,
  log2titre_expected = case_when(
    time_sample < time_infection ~ log2(5),
    time_sample < time_peak ~ log2(5) +
      (long_term_boost + short_term_boost) * (time_sample - time_infection) / time_to_peak,
    time_sample < time_wane ~ log2(5) + long_term_boost +
      short_term_boost * (1 - (time_sample - time_peak) / time_to_wane),
    TRUE ~ log2(5) + long_term_boost
  ),
  log2titre_actual = rnorm(n_individuals, log2titre_expected, 0.5),
  titre_actual = 2^log2titre_actual,
  titre_censored = censor_titres(titre_actual),
  log2titre_censored = log2(titre_censored),
  time_from_infection = time_sample - time_infection
)

write_csv(sim_data, "examples/ex1/sim-data.csv")

titre_plot <- sim_data %>%
  select(pid, time_from_infection, titre_actual, titre_censored) %>%
  pivot_longer(contains("titre"), names_to = "type", values_to = "titre") %>%
  ggplot(aes(time_from_infection, titre)) +
  theme_bw() +
  theme(panel.spacing = unit(0, "null"), strip.background = element_blank()) +
  facet_wrap(~type, labeller = as_labeller(~ recode(.x, "titre_actual" = "Actual", "titre_censored" = "Censored"))) +
  scale_x_continuous("Time from infection", breaks = c(-50, 0, time_to_peak, time_to_peak + time_to_wane, 100)) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:15)) +
  geom_point(alpha = 0.5, shape = 18) +
  geom_segment(aes(x = -50, y = 5, xend = 0, yend = 5), color = "red") +
  geom_segment(aes(
    x = 0, y = 5, xend = time_to_peak, yend = 5 * 2^long_term_boost * 2^short_term_boost
  ), color = "red") +
  geom_segment(aes(
    x = time_to_peak, y = 5 * 2^long_term_boost * 2^short_term_boost,
    xend = time_to_peak + time_to_wane, yend = 5 * 2^long_term_boost
  ), color = "red") +
  geom_segment(aes(
    x = time_to_peak + time_to_wane, y = 5 * 2^long_term_boost,
    xend = 100, yend = 5 * 2^long_term_boost
  ), color = "red")

ggsave("examples/ex1/plot.pdf", titre_plot, width = 15, height = 10, units = "cm")

dyn.load("build/seromalder.so")
.C(
  "sml_r_mcmc",
  n_individuals = 2L,
  logtitres = c(4, 6), times_sample = c(50, 60), times_infection = c(70, 80),
  start_short_term_boost = 1, start_long_term_boost = 1,
  start_time_to_peak = 1, start_time_to_wane = 1,
  iterations = 10L,
  out_long_term_boost = rep(0, 10), out_short_term_boost = rep(0, 10),
  out_time_to_peak = rep(0, 10), out_time_to_wane = rep(0, 10)
)
