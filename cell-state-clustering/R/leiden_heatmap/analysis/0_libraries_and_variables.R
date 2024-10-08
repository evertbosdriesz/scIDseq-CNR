## libraries

library(tidyverse)
library(igraph)           # for creating the network and adjacency matrices 
library(WGCNA)            # for calculating adjacency from distance
library(RANN)             # for calculating nearest neighbours
library(leiden)           # for clustering with the Leiden algorithm
library(yarrr)            # for using the right colours
library(ggplot2)
library(ComplexHeatmap)
library(umap)
library(kmed)
library(ggpubr)
library(cowplot)


## variables

## VARIABLES

# set the input file, this needs to be a matrix with abtibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 

input_file <- 'R/01_data_in_use/TMM_normalised_z_transformed_matrix.csv'

# set the treatments to filter on

treatment_1 <- "ip70S6K_EGF"
treatment_2 <- "EGF"
treatment_3 <- "iRSK_EGF"

# sets the name of the data frame to take the non-changing antibodies from.
# Needs to have the columns 'antibody' and 'sum_ranks'.

abs_no_fold_change <- "R/01_data_in_use/non_significant_abs_kruskal_test.csv"

# set the settings for leiden

no_of_nearest_neighbours <- 15 # set number of nearest neighbours for SSN calculations
res = 1.2 # set resolution for leiden clustering
no_iter <- -2L #set the number of iterations the algorithm is run, standard is 2
no_repeats <- 10


# want to save?

save_figure = "yes"
save_table = "yes"


## CREATE OUTPUT FOLDERS IF NEED BE

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'/R/leiden_heatmap/figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'/R/leiden_heatmap/clean_data/',sep=""), Sys.Date(), sep="")
if (file.exists(data_dir)){print(paste(data_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(data_dir))
  print('generating output directory')
}

## themes and palettes

palette_plots <- c(as.character(piratepal(palette = "basel")), "8A2BE2")

my_theme <- theme(plot.title = element_text(hjust = 0.5, size = 25),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15),
                  axis.text.y = element_text(size = 15),
                  legend.text = element_text(size = 15),
                  legend.title = element_text(size = 15),
                  legend.key.size = unit(0.5, "cm"),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  strip.text.x = element_text(size = 12, color = "black"),
                  strip.background = element_blank(),
                  legend.background = element_rect())

# set counter for suffix after plots and tables

counter = 1
