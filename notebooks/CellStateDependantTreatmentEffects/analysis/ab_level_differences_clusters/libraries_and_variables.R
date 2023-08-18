## libraries

library(tidyverse)
library(yarrr)            # for using the right colours
library(ggplot2)
library(ComplexHeatmap)
library(umap)
library(kmed)
library(ggpubr)
library(cowplot)
library(multcompView)     # for creating the significance labels
# library(ppcor)            # for doing partial correlations


## variables

## VARIABLES

# set the input file, this needs to be a matrix with abtibody signals
# where the first three columns are the meta data and called 'sample_id', 
# 'plate_number' and 'treatment'. 

input_file <- '../../../../data/processed/TMM_normalised_z_transformed_concencus_clusters.csv'

# set the treatments to filter on

treatment_1 <- "ip70S6K_EGF"
treatment_2 <- "EGF"
treatment_3 <- "iRSK_EGF"

# sets the name of the data frame to take the non-changing antibodies from.
# Needs to have the columns 'antibody' and 'sum_ranks'.

abs_no_fold_change <- "../../../../data/annotations/non_significant_abs_kruskal_test.csv"


# want to save?

save_figure = "yes"
save_table = "yes"


## CREATE OUTPUT FOLDERS IF NEED BE

#set figure directory and creates it if it doesn't exist
figure_dir <- paste(paste(getwd(),'/../../figures/',sep=""), Sys.Date(), sep="")
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating output directory')
}

data_dir <- paste(paste(getwd(),'/../../clean_data/',sep=""), Sys.Date(), sep="")
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

## Created functions

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}