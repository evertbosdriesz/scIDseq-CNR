library(tidyverse)
library(tidymodels)
library(here)
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

# All antibodies
ab_list <- sort(unique(dat$ab_name))


## Plot all antibodies.

plot_ab <- function(ab) {
  ggplot(
    filter(dat, ab_name == ab),
    aes(x = treatment, y = zscore, fill = treatment)
  ) +
    #expand_limits(y = 0) +
    my_theme +
    scale_fill_manual(values = col_lst) +
    labs(title = ab, x = "Treatment", y = "z-score") +
    geom_boxplot(lwd=0.25, outlier.size = 0.25) +
    theme(legend.position = "none")
}

save_plot <- function(ab){
  ggsave(
    here("ab-profiles", "population", str_c(ab, "_expression.pdf")),
    plot_ab(ab),
    width = 2, height = 1.75
  )
}

purrr::map(ab_list, save_plot)

## Plot all antibodies, slit out by cell-state-cluster

background_col <- "gray90"

plot_ab_cs <- function(ab) {
  ggplot(
    filter(dat, ab_name == ab),
    aes(x = cluster, y = zscore, fill = treatment)
  ) +
    #expand_limits(y = 0) +
    my_theme +
    scale_fill_manual(values = col_lst) +
    labs(title = ab, x = "Cell state cluster", y = "z-score") +
    geom_rect(xmin=0.5, xmax=1.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
    geom_rect(xmin=2.5, xmax=3.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
    geom_rect(xmin=4.5, xmax=5.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
    geom_rect(xmin=6.5, xmax=7.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
    geom_rect(xmin=8.5, xmax=9.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
    geom_boxplot(lwd=0.25, outlier.size = 0.25) +
    theme(legend.position = "none")
}



save_plot_cs <- function(ab){
  ggsave(
    here("ab-profiles", "per-cell-state", str_c(ab, "_expression.pdf")),
    plot_ab_cs(ab),
    width = 2, height = 1.75
  )
}


purrr::map(ab_list, save_plot_cs)



# Plots for manuscript
#
# ## ERK1/2
# ab <- "ERK1_2_P"
# ggsave(
#   here("figures", "scIDseq", str_c(ab, "_expression.pdf")),
#   plot_ab(ab) + labs(title = "phospho-ERK1/2"),
#   width = 2, height = 1.75
# )
#
# ## Cyclin B1
# ab <- "CYCLIN_B1"
# ggsave(
#   here("figures", "scIDseq", str_c(ab, "_expression.pdf")),
#   plot_ab(ab) + labs(title = "Cyclin B1"),
#   width = 2, height = 1.75
# )
#
# ## MK2
# ab <- "MAPK_APK2_P"
# ggsave(
#   here("figures", "scIDseq", str_c(ab, "_expression.pdf")),
#   plot_ab(ab) + labs(title = "phospho-MK2"),
#   width = 2, height = 1.75
# )
#
# ## p-H3
# ab <- "H3_P"
# ggsave(
#   here("figures", "scIDseq", str_c(ab, "_expression.pdf")),
#   plot_ab(ab) + labs(title = "phospho-H3"),
#   width = 2, height = 1.75
# )
#
#
# ## p-H3
# ab <- "RB_P"
# ggsave(
#   here("figures", "scIDseq", str_c(ab, "_expression.pdf")),
#   plot_ab(ab) + labs(title = "phospho-RB"),
#   width = 2, height = 1.75
# )
#
# plot_ab_lst <- function(ab_lst, labels) {
#   ggplot(
#     mutate(filter(dat, ab_name %in% ab_lst), ab_name = factor(ab_name, levels = levels(ab_lst))),
#     aes(x = cluster, y = zscore, fill = treatment)) +
#     #expand_limits(y = 0) +
#     my_theme +
#     scale_fill_manual(values = col_lst) +
#     labs(title = "Antibody expression", x = "Cell state cluster", y = "z-score") +
#     geom_rect(xmin=0.5, xmax=1.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
#     geom_rect(xmin=2.5, xmax=3.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
#     geom_rect(xmin=4.5, xmax=5.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
#     geom_rect(xmin=6.5, xmax=7.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
#     geom_rect(xmin=8.5, xmax=9.5, ymin=-Inf, ymax=+Inf, fill = background_col, color = background_col) +
#     geom_boxplot(lwd=0.25, outlier.size = 0.25) +
#     theme(legend.position = "none") +
#     facet_wrap(
#       ~ab_name, scales = "free",
#       labeller = as_labeller(labels), nrow = 2)
# }
#
# labels <- c(
#   ERK1_2_P = "phospho-ERK1/2",
#   MAPK_APK2_P = "phospho-MK2",
#   CYCLIN_B1 = "Cyclin B1",
#   H3_P = "phospho-H3"
# )
#
# ab_facts <- factor(names(labels), levels = c("ERK1_2_P", "MAPK_APK2_P",  "CYCLIN_B1", "H3_P"))
# plt <- plot_ab_lst(ab_facts, labels)
#
# ggsave(
#   here("figures","AB_expression.pdf"),
#   plt,
#   width = 4, height = 3
# )
