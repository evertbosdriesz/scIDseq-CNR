## CODE

# load dataframe

TMM_normalised <- read_tsv(file = input_file)

#create list of antibody names
total_abs <- colnames(TMM_normalised[-(1:3)])

# alphebatises list of antibody names
total_abs <- sort(total_abs)

#create seperate data frames for each treatment
table_treatment_1 <- TMM_normalised %>% 
  filter(treatment == treatment_1)

table_treatment_2 <- TMM_normalised %>% 
  filter(treatment == treatment_2)

# Create dataframe with p_values and Fold changes
p_value <- vector()     #creates emtpy vector to later add the p_vals to
antibody <- vector()        #creates emtpy vector to later add the antibody to
fold_change <- vector() #creates emtpy vector to later add the fold change to

for (antibody_of_interest in total_abs) {       #runs through list of abs
  val_treatment_1 <- table_treatment_1 %>% 
    pull(antibody_of_interest)                    #makes vector of values for this treatment with a antibody
  
  val_treatment_2 <- table_treatment_2 %>% 
    pull(antibody_of_interest)                    #makes vector of values for this treatment with a antibody
  
  antibody <- append(antibody, antibody_of_interest)      #adds antibody name to vector antibody
  
  x <- ks.test(val_treatment_2, val_treatment_1)$p.value   #runs statistical test
  p_value <- append(p_value, x)                   #adds p_val to vector of p_value
  
  fold_change <- append(fold_change, mean(val_treatment_1) / mean(val_treatment_2))
}

table_fold_changes <- data.frame(antibody, fold_change, p_value)

table_fold_changes<- table_fold_changes %>% 
  mutate(log2_fold_change = log2(fold_change)) %>% 
  mutate(minus_log10_p_value = -log10(p_value))

# convert infinite values in table to 20, since the -10log conversion 
# transforms very low p-vals to Inf since there's too little memory
# for the computer. This value is chosen for clean plotting.

table_fold_changes <- table_fold_changes %>% 
  mutate(minus_log10_p_value = recode(minus_log10_p_value, "Inf" = 20))

#create tables with only values that interest us

p_val_lim = 5
pos_fold_change_lim = 0
neg_fold_change_lim = -pos_fold_change_lim

#create table with negative fold changes
table_negative_fold_changes <- table_fold_changes %>% 
  filter(log2_fold_change <= neg_fold_change_lim) %>% 
  filter(minus_log10_p_value >= 5)

#create table with positive fold changes
table_positive_fold_changes <- table_fold_changes %>% 
  filter(log2_fold_change >= pos_fold_change_lim) %>% 
  filter(minus_log10_p_value >= 5)

#create table with values that don't change
table_no_fold_changes <- table_fold_changes %>%
  filter(log2_fold_change <= pos_fold_change_lim) %>%
  filter(log2_fold_change >= neg_fold_change_lim) %>%
  filter(minus_log10_p_value <= 5)

## PLOTTING

pdf(file = paste(figure_dir, '/vulcano_plot_', treatment_2, "_vs_", treatment_1, ".pdf", sep = ""),
    width = 4,
    height = 4)

# create vulcano plot
ggplot() +
  geom_point(data=table_fold_changes, 
             aes(x = log2_fold_change, 
                 y=minus_log10_p_value),
             size = 2,
             colour = "darkgreen") +
  geom_point(data=table_fold_changes, 
             aes(x = log2_fold_change, 
                 y=minus_log10_p_value),
             size = 2,
             shape = 21,
             colour = "black") +
  geom_hline(yintercept= p_val_lim, linetype = "dashed", color = "gray") +
  coord_cartesian(xlim = c(-1,1), ylim=c(0, 22)) +
  geom_point(data= table_positive_fold_changes,
             aes(x = log2_fold_change, y=minus_log10_p_value),
             colour = "red", show.legend = FALSE,
             size = 2) +
  geom_point(data= table_positive_fold_changes,
             aes(x = log2_fold_change, y=minus_log10_p_value),
             colour = "darkred", show.legend = FALSE,
             size = 2,
             shape = 21) +
  geom_text_repel(data = table_positive_fold_changes,
                  aes(x = log2_fold_change, y=minus_log10_p_value,
                      label = antibody),
                  xlim = c(pos_fold_change_lim, NA), 
                  size = 2.5, show.legend = FALSE) +
  geom_point(data= table_negative_fold_changes,
             aes(x = log2_fold_change, y=minus_log10_p_value),
             colour = "blue", 
             show.legend = FALSE,
             size = 2) +
  geom_point(data= table_negative_fold_changes,
             aes(x = log2_fold_change, y=minus_log10_p_value),
             colour = "darkblue", 
             show.legend = FALSE,
             size = 2,
             shape = 21) +
  geom_text_repel(data = table_negative_fold_changes,
                  aes(x = log2_fold_change, y=minus_log10_p_value,
                      label = antibody),
                  xlim = c(NA, neg_fold_change_lim),
                  size = 2.5, show.legend = FALSE) +
  ggtitle(paste(treatment_2, "vs.", treatment_1)) +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  theme(panel.background = element_rect(fill = "white", colour = "black", line = "solid"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA)) +
  theme(text = element_text(size = 9)) +
  theme(plot.title = element_text(hjust=0.5, size = 10))

dev.off()
