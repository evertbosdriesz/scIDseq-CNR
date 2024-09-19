## Libraries

library(tidyverse)
library(yarrr)            # for using the right colours



## VARIABLES

# set the input file, this needs to be a matrix with abtibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 

# set the treatments to filter on

treatment_1 <- "ip70S6K_EGF"
treatment_2 <- "EGF"
treatment_3 <- "iRSK_EGF"


# want to save?

save_figure = "yes"


## CREATE OUTPUT FOLDERS IF NEED BE

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'/R/umap_distributions/figures/', sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'/R/umap_distributions/clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}
