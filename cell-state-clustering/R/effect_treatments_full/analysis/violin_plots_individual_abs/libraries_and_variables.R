## libraries

library(tidyverse)
library(yarrr)            # for using the right colours
library(ggpubr)           # for creating nice plots
library(cowplot)

## VARIABLES

# set the input file, this needs to be a matrix with abtibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 

input_file <- 'R/01_data_in_use/TMM_normalised_z_transformed_matrix.csv'

# set the treatments you want to compare (choose from ip70S6K_EGF, iRSK_EGF and EGF)

treatment_1 <- "ip70S6K_EGF" # the test treatment
treatment_2 <- "EGF"      # control

# automatically sets the right colour for later analysis
if (treatment_1 == "iRSK_EGF") {
  correct_colour = 3
} 
if (treatment_1 == "ip70S6K_EGF") {
  correct_colour = 2
}


# want to save?

save_figure = "yes"


## CREATE OUTPUT FOLDERS IF NEED BE

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'/R/effect_treatments_full/figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'/R/effect_treatments_full/clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}

## themes and palettes

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

counter <- 1
