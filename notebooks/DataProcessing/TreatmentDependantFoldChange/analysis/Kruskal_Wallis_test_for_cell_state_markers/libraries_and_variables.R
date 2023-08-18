library(tidyverse)

## VARIABLES

# set the input file, this needs to be a matrix with z-transformed
# antibody signals where the first three columns are the meta data and 
# called 'sample_id', 'plate_number' and 'treatment'. 

input_file <- '../../../../../data/processed/TMM_normalised_matrix.tsv'

#set the treatments you compare
treatment_1 <- "EGF"    
treatment_2 <- "iRSK_EGF"
treatment_3 <- "ip70S6K_EGF"

save_table = "yes"


## CREATE OUTPUT FOLDERS IF NEED BE

data_dir <- paste(paste(getwd(),'/../../clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}
