library(tidyverse)

cdc_viruses <- read_csv("data/cdc-virus.csv", col_types = cols()) %>%
  mutate(
    to_search = virus_full %>%
      str_replace("/\\d{2}e?$", "/") %>%
      str_replace_all(" ", "_") %>%
      recode("a/singapore/16-0019/" = "a/singapore/[^/]*16-0019")
  )

combined_fasta <- read_tsv(
  "data/combined-gisaid-sequences.tsv",
  col_names = FALSE, col_types = cols()
) %>%
  select(X1, X2) %>%
  rename(meta = X1, seq = X2) %>%
  mutate(isolate_name = str_replace(meta, "^.*\\|.*\\|(.*)$", "\\1"))

fasta_virus_names <- unique(combined_fasta$isolate_name)

cdc_isolate_names <- map_dfr(
  cdc_viruses$to_search,
  ~ tibble(
    to_search = .x,
    fasta_virus_name = fasta_virus_names[str_detect(tolower(fasta_virus_names), .x)]
  )
)

cdc_isolate_names_filtered <- cdc_isolate_names %>%
  inner_join(cdc_viruses, "to_search") %>%
  mutate(
    fasta_virus_year = str_replace(fasta_virus_name, "^.*/(\\d{4})$", "\\1") %>%
      str_replace("^.*/(\\d{4})_.*$", "\\1") %>%
      as.integer()
  ) %>%
  filter(fasta_virus_year == virus_year) %>%
  group_by(virus_full) %>%
  # TODO(sen) Replace with something better (passage filter?)
  filter(str_length(fasta_virus_name) == min(str_length(fasta_virus_name))) %>%
  filter(row_number() == 1) %>%
  ungroup()

cdc_isolate_names_filtered %>%
  arrange(virus_full) %>%
  print(n = 100)

write_csv(cdc_isolate_names_filtered, "virus-distances/cdc-isolate-names.csv")

cdc_fasta <- combined_fasta %>%
  filter(isolate_name %in% cdc_isolate_names_filtered$fasta_virus_name) %>%
  group_by(isolate_name) %>%
  filter(row_number() == 1) %>%
  ungroup()

write_tsv(cdc_fasta %>% select(meta, seq), "virus-distances/cdc-fasta.tsv", col_names = FALSE)

system("virus-distances/tsvtofasta.sh")

system("virus-distances/align.sh")
