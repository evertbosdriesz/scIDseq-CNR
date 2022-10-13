## Clean up antibody and treatment names, and perform TMM normalization

library(tidyverse)
library(here)


# load data -------------------------------------------------------------------

dat <- read_tsv(here("data", "raw", "REVI182_raw_Evert.txt"))
# This contains the cells and antibodies that are kept after QC performend by
# van Eijl after initial data-analysis with subsampling
keep_cells <- read_tsv(here("data", "raw", "scIDseqData_vanEijl.txt"))
keep_cells <- unique(keep_cells$sampleID)

# Clean up the raw read counts -----------------------------------------------
dat_raw <- dat %>%
  filter(sample_folder_well_name %in% keep_cells) %>%
  filter(sample_folder_well_name %in% keep_cells) %>%
  select(-starts_with('X')) %>% # Remove last 4 spurious columns
  # Rename treatment
  mutate(
    treatment = case_when(
      treatment == "+EGF_ip70S6K" ~ "ip70S6K_EGF",
      treatment == "-EGF_vehicle" ~ "No_EGF",
      treatment == "+EGF_iRSK1/2/3/4" ~ "iRSK_EGF",
      treatment == "+EGF_vehicle" ~ "EGF",
      TRUE ~ NA_character_
    )
  ) %>%
  # Some antibody name cleaning
  mutate(
    Ab_name = str_to_upper(Ab_name) %>%
      str_replace_all("-", "_")
  ) %>%
  mutate(Ab_name = str_replace(Ab_name, "AKT_P_4060", "AKT1_P")) %>%
  mutate(Ab_name = str_replace(Ab_name, "EGFR_P_AF1095", "EGFR_P_Y1173")) %>%
  mutate(Ab_name = str_replace(Ab_name, "EGFR_P$", "EGFR_P_Y1045")) %>%
  mutate(Ab_name = str_replace(Ab_name, "GDK3B_P", "GSK3B_P")) %>%
  mutate(Ab_name = str_replace(Ab_name, "IKB_A", "IKBA")) %>%
  mutate(Ab_name = str_replace(Ab_name, "NOTCH_1_CLEAVED", "NICD")) %>%
  mutate(Ab_name = str_replace(Ab_name, "RSK_P$", "RSK1_P_S380")) %>%
  mutate(Ab_name = str_replace(Ab_name, "RSK_P90_P", "RSK1_P_T359")) %>%
  mutate(Ab_name = str_replace(Ab_name, "SMAD2_3", "SMAD2_3_P")) %>%
  mutate(Ab_name = str_replace(Ab_name, "STAT1_P_P91", "STAT1")) %>%
  mutate(Ab_type = if_else(str_detect(Ab_name, '_P'), "phospho", "total")) %>%
  # Select and rename columns
  select(
    sample_id = sample_folder_well_name,
    plate_number = sample_folder,
    treatment,
    ab_name = Ab_name,
    ab_type = Ab_type,
    ab_count_raw = antibody_count,
    treatment
    ) %>%
  # Remove low qualitiy cells and Abs
  filter(sample_id %in% keep_cells) %>%
  filter(ab_name != "CDC2_P")


# Perform TMM normalisation ----------------------------------------------------

# Get the normalisation factors using TMM
mat <-
  dat_raw %>%
  select(sample_id, ab_name, ab_count_raw) %>%
  pivot_wider(names_from = ab_name, values_from = ab_count_raw ) %>%
  filter(sample_id %in% keep_cells) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

normfactors <- edgeR::calcNormFactors(t(mat), sumTrim = 0.05, logratioTrim = 0) %>%
  enframe("sample_id", "normfactor")

# Extract sample annotations
sample_annot <-
  dat_raw %>%
  group_by(sample_id) %>%
  mutate(lib_size = sum(ab_count_raw)) %>%
  select(sample_id, plate_number, treatment, lib_size) %>%
  ungroup() %>%
  distinct()


# Scale counts to TMM-corrected library sizes
dat_tmm <- dat_raw %>%
  left_join(normfactors, by = "sample_id") %>%
  left_join(sample_annot) %>%
  mutate(ab_count_tmm_batch = ab_count_raw/(lib_size*normfactor))

# Remove the plate batch effects
batch_corrected <- select(dat_tmm, sample_id, ab_name, ab_count_tmm_batch) %>%
  pivot_wider(names_from = ab_name, values_from = ab_count_tmm_batch) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

batch = column_to_rownames(sample_annot, var = "sample_id")[
  rownames(batch_corrected),
  "plate_number"
]

batch_corrected <-
  limma::removeBatchEffect(t(batch_corrected), batch = batch) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "ab_name", values_to = "ab_count_tmm")

dat_tmm <- dat_tmm %>%
  left_join(batch_corrected)

# Write the results to file ----------------------------------------------------
write_tsv(dat_tmm, here("data", "processed", "scIDseq-vanEijl-tmm.tsv"))
write_tsv(dat_raw, here("data", "processed", "scIDseq-venEijl-raw-counts.tsv"))
