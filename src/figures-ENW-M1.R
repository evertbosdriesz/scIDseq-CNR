library(tidyverse)
library(tidymodels)
library(here)

source(here("src", "graphics-options.R"))

dat_tmm <- read_tsv(here("data", "processed", "scIDseq-vanEijl-tmm.tsv")) %>%
  select(sample_id, treatment, ab_name, ab_type, ab_count_tmm) %>%
  filter(treatment != "No_EGF") %>%
  mutate(treatment = factor(treatment, levels = c("EGF", "ip70S6K_EGF", "iRSK_EGF")))


dat_zscore <- read_csv(
  here("data", "processed" ,"TMM_normalised_z_transformed_concencus_clusters.csv")) %>%
  select(-`...1`, -plate_number, -treatment) %>%
  mutate(cluster = as_factor(cluster)) %>%
  pivot_longer(!c(cluster, sample_id), names_to = "ab_name", values_to = "zscore")


dat <- left_join(dat_tmm, dat_zscore) %>%
  mutate(C5 = cluster == "5")


plot_ab <- function(ab) {
  ggplot(
    filter(dat, ab_name == ab),
    aes(x = C5, y = zscore, fill = treatment)) +    #expand_limits(y = 0) +
    my_theme +
    scale_fill_manual(values = col_lst) +
    labs(title = ab, x = "Cell state cluster", y = "z-score") +
    scale_x_discrete(limits = c(TRUE, FALSE), labels = c("G2M", "Others")) +
    geom_boxplot(lwd=0.25, outlier.size = 0.25) #+
    #theme(legend.position = "bottom")
  }

plot_ab("CYCLIN_B1") + labs(title = "Cyclin B1 expression")
# ggsave(
#   "~/My Drive/Proposals/2021-4-ENW-M1/cyclin-b1.pdf",
#   plot_ab("CYCLIN_B1"),
#   width = 2.5, height = 3
# )

plot_correlation <- function(df, ab1, ab2, c, count_col = 'zscore'){
  df %>%
    select(sample_id, ab_name, treatment, cluster, !!sym(count_col)) %>%
    pivot_wider(
      names_from = ab_name,
      values_from = !!sym(count_col)
    ) %>%
    mutate(in_cluster = if_else(cluster == c, TRUE, FALSE)) %>%
    ggplot(aes(x = !!sym(ab1), y = !!sym(ab2), col = in_cluster)) +
    geom_point() +
    # facet_wrap(~treatment) +
    labs(x = ab1, y = ab2, title = str_c("Cluster ", c)) +#, subtitle = count_col) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = c("gray", "red"), name = str_c("In cluster ", c)) +
    my_theme
}

plot_correlation(dat,"MAPK_APK2_P",  "CYCLIN_B1", "5") +
  labs(title = "Cyclin B1 is upregulated through MK2 in G2M cells")


plot_correlation(dat,"H2A_P",  "RB_P", "5")
# ggsave(
#   "~/My Drive/Proposals/2021-4-ENW-M1/MK2-vs-CylinB1.pdf",
#   plot_correlation(dat,"MAPK_APK2_P",  "CYCLIN_B1", "5") +
#     labs(title = "Cyclin B1 is upregulated through MK22\nin G2M cells",
#           x = "MK2", y = "Cyclin B1") +
#     theme(legend.position = "none"),
#   width = 2.5, height = 2
# )
#
ggsave(
  here("figures", "presentations", "MK2-vs-CylinB1.pdf"),
  plot_correlation(dat,"MAPK_APK2_P",  "CYCLIN_B1", "5") +
    labs(title = "", x = "MK2", y = "Cyclin B1") +
    presentation_theme +
    annotate(geom="text", x=3, y=5.5, label="Cell state 5", color="red") +
    annotate(geom="text", x=3, y=.8, label="All other cells", color="dark gray") +
    theme(legend.position = "none"),
  width = 3, height = 3
)

ab_list <- sort(unique(dat_tmm$ab_name))
sorted(ab_list)
