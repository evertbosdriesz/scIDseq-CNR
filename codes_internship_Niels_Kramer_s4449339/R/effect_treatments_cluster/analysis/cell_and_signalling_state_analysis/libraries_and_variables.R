## libraries

library(tidyverse)
library(yarrr)            # for using the right colours
library(ggpubr)
library(cowplot)


## VARIABLES

# set the input file, this needs to be a matrix with abtibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 

input_file <- 'R/01_data_in_use/TMM_normalised_z_transformed_concencus_clusters.csv'

# set the file to take the categories of signals from

model <- "R/01_data_in_use/model_cell_states.csv"
cycling_markers <- "R/01_data_in_use/model_cell_cycle.csv"
umap_coord <- "R/01_data_in_use/TMM_normalised_z_transformed_clusters_with_umap.tsv"

# want to save?

save_figure = "yes"
save_table = "yes"


## CREATE OUTPUT FOLDERS IF NEED BE

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'/R/effect_treatments_cluster/figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'/R/effect_treatments_cluster/clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}

## themes and palettes

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")
