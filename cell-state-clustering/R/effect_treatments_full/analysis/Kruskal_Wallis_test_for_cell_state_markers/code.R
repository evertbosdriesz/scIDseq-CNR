## CODE

# load dataframe

TMM_normalised <- read_tsv(file = input_file)

#create list of antibody names
total_abs <- colnames(TMM_normalised[-(1:3)])

# alphebatises list of antibody names
total_abs <- as.factor(sort(total_abs))

TMM_normalised <- TMM_normalised %>% 
  filter(treatment == treatment_1 | treatment == treatment_2 | treatment == treatment_3)

# Create dataframe with p_values and Fold changes
p_value <- vector()     #creates emtpy vector to later add the p_vals to
antibody <- vector()        #creates emtpy vector to later add the antibody to

for (antibody_of_interest in total_abs) {       #runs through list of abs
    TMM_temp <- vector("list", 3)
    
    EGF <- TMM_normalised %>% 
                filter(treatment == treatment_1) %>% 
                pull(antibody_of_interest)
    
    TMM_temp[[1]] <- EGF
    
    ip70S6K <- TMM_normalised %>% 
      filter(treatment == treatment_2) %>% 
      pull(antibody_of_interest)
    
    TMM_temp[[2]] <- ip70S6K
    
    iRSK <- TMM_normalised %>% 
      filter(treatment == treatment_3) %>% 
      pull(antibody_of_interest)
  
    TMM_temp[[3]] <- iRSK
    
  antibody <- append(antibody, antibody_of_interest)      #adds antibody name to vector antibody
  
  
  x <- kruskal.test(TMM_temp)
  
  p_value <- append(p_value, x[["p.value"]])                   #adds p_val to vector of p_value

}

# make a data frame out of the antibody and p-values
table_effects <- data.frame(antibody, p_value)

# transform the p-value in the -log10 p-value
table_effects <- table_effects %>% 
  mutate(minus_log10_p_value = -log10(p_value))

# save table as csv for later analysis
write.csv(table_effects, paste(data_dir, "/table_significance_all_abs_kruskal_test.csv", sep = ""))

#create tables with only values that interest us

p_val_lim = 5

#create table with no significant effects
table_no_effects <- table_effects %>% 
  filter(minus_log10_p_value <= p_val_lim)

write.csv(table_no_effects, paste(data_dir, "/table_cell_state_markers_kruskal_test.csv", sep = ""))

