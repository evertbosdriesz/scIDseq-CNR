my_theme <- theme(
  # Overall
  text = element_text(colour = "black", size = 6, family = 'Helvetica'),
  line = element_line(colour = "black", size = .3),
  title = element_text(colour = "black", size = 6, family = 'Helvetica'),
  axis.title = element_text(colour = "black"),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  # Axis
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(colour = "black", size = 6, family = 'Helvetica'),
  # Legend
  legend.text = element_text(colour = "black", size = 6, family = 'Helvetica'),
  legend.key = element_blank(),
  # Facet labels
  strip.text = element_text(colour = "black", size = 6, family = 'Helvetica',
                            hjust = 0),
  strip.background = element_rect(fill = NA)
)

presentation_theme <- theme(
  # Overall
  text = element_text(colour = "black", size = 10, family = 'Helvetica'),
  #line = element_line(colour = "black", size = 1),
  title = element_text(colour = "black", size = 10, family = 'Helvetica'),
  axis.title = element_text(colour = "black"),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  # Axis
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(colour = "black", size = 10, family = 'Helvetica'),
  # Legend
  legend.text = element_text(colour = "black", size = 8, family = 'Helvetica'),
  legend.key = element_blank(),
  # Facet labels
  strip.text = element_text(colour = "black", size = 8, family = 'Helvetica',
                            hjust = 0),
  strip.background = element_rect(fill = NA)
)

col_lst <- c('EGF' = 'forestgreen', # "#228B22"
             'ip70S6K_EGF' = 'lightblue', # HEX code "#ADD8E6"
             'iRSK_EGF' = 'firebrick') # Hex code"#B22222"

tx_labs <- c('EGF' = 'Control',
             'ip70S6K_EGF' = 'p70-S6Ki',
             'iRSK_EGF' = 'RSKi')
