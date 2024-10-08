library(here)
library(tidyverse)
source(here("src", "graphics-options.R"))



dat_tmm <- read_tsv(here("data", "processed", "scIDseq-vanEijl-tmm.tsv")) %>%
  select(sample_id, treatment, ab_name, ab_type, ab_count_tmm) %>%
  filter(treatment != "No_EGF") %>%
  mutate(treatment = factor(treatment, levels = c("EGF", "ip70S6K_EGF", "iRSK_EGF")))


dat_zscore <- read_csv(
  here("data", "processed" ,"TMM_normalised_z_transformed_concencus_clusters.csv")
) %>%
  select(-...1, -plate_number, -treatment) %>%
  mutate(cluster = as_factor(cluster)) %>%
  pivot_longer(!c(cluster, sample_id), names_to = "ab_name", values_to = "zscore")


dat <- left_join(dat_tmm, dat_zscore)

plot_correlation <- function(df, ab1, ab2, c, count_col = 'zscore'){
  df %>%
    select(sample_id, ab_name, treatment, cluster, !!sym(count_col)) %>%
    pivot_wider(
      names_from = ab_name,
      values_from = !!sym(count_col)
    ) %>%
    mutate(in_cluster = if_else(cluster == c, TRUE, FALSE)) %>%
    ggplot(aes(x = !!sym(ab1), y = !!sym(ab2))) +
    geom_point(aes(col = treatment), size = 0.5) +
    facet_wrap(~in_cluster) +
    labs(x = ab1, y = ab2, title = str_c("Cluster ", c), subtitle = count_col) +
    geom_smooth(method = "lm") +
    presentation_theme +
    theme_light() +
    scale_color_manual(values = col_lst)
}

plt <- plot_correlation(dat,"MAPK_APK2_P",  "CYCLIN_B1", "5") +
  labs(x = "APK2 (z-score normalized expression)", y = "Cyclin B1  (z-score normalized expression)",
       subtitle = "") +
  theme(legend.position = "none") +
  theme_light() +
  manuscript_theme
  scale_color_manual(values = col_lst, name = "Treatment")
plt
ggsave(
    here("figures","APK2-CYCLINB1-correlation.pdf"),
    plt,
    width = 4, height = 2
  )

