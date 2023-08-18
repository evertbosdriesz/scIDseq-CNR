# LIBRARIES

library(tidyverse)
library(ggridges) # for plotting the graph type with shifted y-positions
library(yarrr)


palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")


## VARIABLES

# set the input file, this needs to be a matrix with z-transformed
# antibody signals where the first three columns are the meta data and 
# called 'sample_id', 'plate_number' and 'treatment'. 

input_file <- '../../../../data/TMM_normalised_z_transformed_matrix.csv'

# set whether you want to save the output 

save_figures = "yes" #choose "yes" if you want to save the figures

## CREATE OUTPUT FOLDERS IF NEED BE

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}
