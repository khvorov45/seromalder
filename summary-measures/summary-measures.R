library(tidyverse)

cdc_viruses <- read_csv("data/cdc-virus.csv")
cdc_vaccine <- read_csv("data/cdc-vaccine.csv")
cdc_obj1_participants <- read_csv("data/cdc-obj1-participant.csv")

cdc_isolate_names <- read_csv("virus-distances//cdc-isolate-names.csv")

cdc_virus_dist_json <- jsonlite::read_json("virus-distances/cdc-fasta-distances.json")

extract_isolate_name <- \(meta) str_replace(meta, "^.*\\|.*\\|(.*)$", "\\1")

cdc_virus_dist <- imap_dfr(cdc_virus_dist_json$nodes, function(node_value, node_name) {
  imap_dfr(node_value$distance, function(distance, node2_name) {
    tibble(
      vaccine_isolate_name = extract_isolate_name(node_name),
      isolate_name = extract_isolate_name(node2_name),
      distance
    )
  })
})

cdc_virus_dist_extra <- cdc_virus_dist %>%
  inner_join(
    cdc_isolate_names %>%
      select(vaccine_isolate_name = fasta_virus_name, vaccine_virus = virus_full),
    c("vaccine_isolate_name")
  ) %>%
  inner_join(
    cdc_isolate_names %>%
      select(isolate_name = fasta_virus_name, virus_full),
    c("isolate_name")
  ) %>%
  inner_join(cdc_vaccine, c("vaccine_virus" = "virus_full"))

cdc_obj1_hi <- read_csv("data/cdc-obj1-hi.csv") %>%
  inner_join(cdc_obj1_participants, "pid") %>%
  inner_join(cdc_viruses, "virus_full") %>%
  left_join(cdc_vaccine %>% mutate(vaccine_strain = TRUE), c("study_year", "virus_full")) %>%
  mutate(vaccine_strain = replace_na(vaccine_strain, FALSE)) %>%
  group_by(pid, timepoint) %>%
  mutate(distance_from_vaccine = abs(virus_year - virus_year[vaccine_strain])) %>%
  ungroup() %>%
  left_join(cdc_virus_dist_extra, c("virus_full", "study_year"))

cdc_obj1_ratios <- cdc_obj1_hi %>%
  filter(timepoint %in% c("Pre-vax", "Post-vax")) %>%
  select(-bleed_date) %>%
  pivot_wider(names_from = "timepoint", values_from = titre) %>%
  mutate(post_pre_ratio = `Post-vax` / `Pre-vax`)

cdc_obj1_ratios %>% arrange(desc(post_pre_ratio))

cdc_obj1_ratios_plot <- cdc_obj1_ratios %>%
  filter(!is.na(post_pre_ratio), !is.na(distance)) %>%
  mutate(
    y_pos = exp(rnorm(n(), log(post_pre_ratio), 0.1)),
    x_pos = rnorm(n(), distance, 0.5)
  ) %>%
  ggplot(aes(x_pos, y_pos)) +
  theme_bw() +
  theme(panel.spacing = unit(0, "null")) +
  scale_y_log10() +
  facet_wrap(~prior_vacs, ncol = 1, strip.position = "right") +
  geom_line(aes(group = pid), alpha = 0.2) +
  geom_point(alpha = 0.2, shape = 18)

ggsave(
  "summary-measures/cdc-obj1-ratios.pdf", cdc_obj1_ratios_plot,
  width = 15, height = 25, units = "cm"
)

sim_data <- tibble(pid = 1:100) %>%
  slice(rep(1:n(), each = 10)) %>%
  mutate(
    virus = rep(1:10, 100),
    virus_distance_from_ref = virus - 1,
    expected_ratio = 1 + (3 - 1) * exp(-0.05 * virus_distance_from_ref),
    log_expected_ratio = log(expected_ratio),
    random_effect = rnorm(n(), 0.05),
    log_ratio = rnorm(n(), log_expected_ratio + random_effect, 0.1),
    ratio = exp(log_ratio)
  )

sim_data %>%
  ggplot(aes(virus_distance_from_ref, expected_ratio)) +
  scale_y_log10() +
  geom_point()
