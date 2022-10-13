# Map kinase substrate relations using the phosphosite-plus database.
# Kinase_Substrate_Dataset.gz is downloaded from https://www.phosphosite.org/ on 2021-02-22

library(tidyverse)
library(here)


df_psp <- read_tsv(here("data/external/Kinase_Substrate_Dataset.gz"), skip = 3) %>%
  filter(SUB_ORGANISM == "human" & KIN_ORGANISM == "human")


df_reg <- read_tsv(here("data/external/Regulatory_sites.gz"), skip = 3) %>%
  filter(ORGANISM == "human")

df_panel <- read_tsv(here("data/annotations/node2gene.txt"))

select(df_panel, KIN_NODE = NODE_NAME, KIN_GENE = GENE)

genes_in_panel <- unique(filter(df_panel, !is.na(GENE))$GENE)

df <-
  df_panel %>%
  select(-SUBSTRATE) %>%
  rename(SUB_NODE_NAME = NODE_NAME, SUB_GENE = GENE, SUB_MOD_RSD = PSITE) %>%
  left_join(
    select(df_psp, KIN_GENE = GENE, SUB_GENE, SUB_MOD_RSD, KINASE, SUBSTRATE)
   ) %>%
  filter(SUB_GENE %in% genes_in_panel & KIN_GENE %in% genes_in_panel) %>%
  left_join(
    select(df_panel, KIN_NODE_NAME = NODE_NAME, KIN_GENE = GENE), by = "KIN_GENE"
   ) %>%
  filter(SUB_NODE_NAME != KIN_NODE_NAME) %>% # Remove auto-phosphorylation
  select(SUB_NODE_NAME, KIN_NODE_NAME, SUB_MOD_RSD,
         SUBSTRATE, KINASE,
         SUB_GENE, KIN_GENE
         ) %>%
  distinct()

write_tsv(df, here("data", "annotations", "substrate-kinase-psp.tsv"))
