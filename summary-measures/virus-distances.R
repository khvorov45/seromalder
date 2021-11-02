library(tidyverse)

cdc_viruses <- read_csv("summary-measures/cdc-virus.csv")

system("gisaid-data/combine-fasta.sh")

read_fasta_tsv <- function(path) {
  read_tsv(path, col_names = FALSE) %>%
    select(X1, X2) %>%
    rename(meta = X1, seq = X2) %>%
    mutate(isolate_name = str_replace(meta, "^.*\\|.*\\|(.*)$", "\\1"))
}

combined_fasta <- read_fasta_tsv("gisaid-data/combined-fasta.tsv")

cdc_viruses_extra <- cdc_viruses %>%
  mutate(
    to_search = virus_full %>%
      str_replace("/\\d{2}e?$", "/") %>%
      str_replace_all(" ", "_") %>%
      recode("a/singapore/16-0019/" = "a/singapore/16")
  )

cdc_viruses_extra %>%
  arrange(to_search) %>%
  print(n = 100)

fasta_virus_names <- unique(combined_fasta$isolate_name)

cdc_isolate_names <- map_dfr(
  cdc_viruses_extra$to_search,
  ~ tibble(to_search = .x, fasta_virus_name = fasta_virus_names[str_detect(tolower(fasta_virus_names), .x)])
)

cdc_isolate_names_filtered <- cdc_isolate_names %>%
  inner_join(cdc_viruses_extra, "to_search") %>%
  mutate(
    fasta_virus_year = str_replace(fasta_virus_name, "^.*/(\\d{4})$", "\\1") %>%
      str_replace("^.*/(\\d{4})_.*$", "\\1") %>%
      as.integer()
  ) %>%
  filter(fasta_virus_year == virus_year) %>%
  group_by(virus_full) %>%
  filter(str_length(fasta_virus_name) == min(str_length(fasta_virus_name))) %>%
  filter(row_number() == 1) %>%
  ungroup()

cdc_isolate_names_filtered %>%
  arrange(virus_full) %>%
  print(n = 100)

cdc_fasta <- combined_fasta %>%
  filter(isolate_name %in% cdc_isolate_names_filtered$fasta_virus_name) %>%
  group_by(isolate_name) %>%
  filter(row_number() == 1) %>%
  ungroup()

write_tsv(cdc_fasta %>% select(meta, seq), "summary-measures/cdc-fasta.tsv", col_names = FALSE)

system("summary-measures/tsvtofasta.sh")

system("summary-measures/align.sh")

cdc_fasta_aligned <- read_fasta_tsv("summary-measures/cdc-fasta-aligned.tsv")

cdc_vaccine_fasta_aligned <- cdc_fasta_aligned %>%
  filter(isolate_name %in% c("A/Singapore/16-0059/2016", "A/Hong_Kong/4801/2014")) %>%
  group_by(isolate_name) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(
    virus_full = recode(
      isolate_name,
      "A/Hong_Kong/4801/2014" = "a/hong kong/4801/14e",
      "A/Singapore/16-0059/2016" = "a/singapore/16-0019/16e"
    )
  ) %>%
  inner_join(cdc_vaccine, "virus_full") %>%
  select(seq, virus_full, study_year) %>%
  rename(vaccine_seq = seq, vaccine_virus_full = virus_full)

cdc_vaccine_diffs <- cdc_fasta_aligned %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(study_year = rep(1:3, length.out = n())) %>%
  inner_join(cdc_vaccine_fasta_aligned, "study_year") %>%
  mutate(
    vaccine_diff_count = map2_dbl(str_split(seq, ""), str_split(vaccine_seq, ""), ~ sum(.x != .y)),
  )

setdiff(cdc_vaccine_diffs$isolate_name, cdc_isolate_names_filtered$fasta_virus_name)
setdiff(cdc_isolate_names_filtered$fasta_virus_name, cdc_vaccine_diffs$isolate_name)

cdc_vaccine_diffs_extra <- cdc_vaccine_diffs %>%
  inner_join(cdc_isolate_names_filtered, c("isolate_name" = "fasta_virus_name"))

write_csv(
  cdc_vaccine_diffs_extra %>% select(virus_full, study_year, vaccine_diff_count),
  "summary-measures/cdc-virus-distances.csv"
)
