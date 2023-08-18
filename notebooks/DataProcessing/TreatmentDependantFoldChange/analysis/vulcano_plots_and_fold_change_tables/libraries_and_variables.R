# LIBRARIES
library(tidyverse)
library(ggrepel)


## VARIABLES

# set the input file, this needs to be a matrix with z-transformed
# antibody signals where the first three columns are the meta data and 
# called 'sample_id', 'plate_number' and 'treatment'. 

input_file <- '../../../../../data/processed/TMM_normalised_matrix.tsv'

#set the treatments you compare (choose from iRSK_EGF, ip70S6K_EGF and EGF)
treatment_1 <- "iRSK_EGF"    # treatment values
treatment_2 <- "EGF" # control values

figure_dir <- paste(paste(getwd(),'/../../figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'/../../clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}