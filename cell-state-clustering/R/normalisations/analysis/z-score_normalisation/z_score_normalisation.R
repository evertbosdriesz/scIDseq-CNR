## LIBRARIES

library(tidyverse)



## VARIABLES

# set the input file, this needs to be a matrix with antibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 
input_file <- 'R/01_data_in_use/TMM_normalised_matrix.tsv'

treatment_1 <- "ip70S6K_EGF"
treatment_2 <- "EGF"
treatment_3 <- "iRSK_EGF"

output_file <- "TMM_normalised_z_transformed_matrix.csv"

save_table <- "yes"

data_dir <- paste(paste(getwd(),'/R/normalisations/normalised_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}


## CODE

# load data frame
TMM_normalised <- read_tsv(file = input_file)

TMM_normalised <- TMM_normalised %>% 
  filter(treatment == treatment_1| treatment == treatment_2 | treatment == treatment_3)

#create alphabetised list of antibody names (exclude meta data)
total_abs <- colnames(TMM_normalised[-(1:3)]) %>% 
  sort()

# turn list of antibody names into factor in order for looping to go right
total_abs <- as_factor(total_abs)

#create a new data frame to store the z-transformed values
TMM_z_score <- TMM_normalised

## TRANSFORM THIS PIECE OF CODE TO A LOOP SO THE ENTIRE DATAFRAME 
## IS NORMALISED PER COLUMN IN ORDER FOR EACH COLUMN TO GET THE SAME WEIGHT

for (antibody in total_abs){
  TMM_z_score[[antibody]] <- scale(TMM_z_score[[antibody]], center = TRUE, scale = TRUE)
}

# save new data frame in the normalised data_folder
if(save_table == "yes") {
  write.csv(TMM_z_score, file = paste(data_dir, "/", output_file, sep = ""))
  print("Table has been saved in your normalised data folder")
} else {
  print("Table has not been saved")
}

# get antibody list
antibodies <- colnames(TMM_normalised[4:length(TMM_normalised)])

# transform into tibble

antibodies <- tibble(antibodies)

# save new data frame in the normalised data_folder
if(save_table == "yes") {
  write.csv(antibodies, file = paste(data_dir, "/antibodies.csv", sep = ""))
  print("Table has been saved in your normalised data folder")
} else {
  print("Table has not been saved")
}

