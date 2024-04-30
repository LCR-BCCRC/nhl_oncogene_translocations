library(tidyverse)
library(readxl)
library(circlize)
library(data.table)
library(ComplexHeatmap)

source("src/circos/prepare_circos.R")
source("src/circos/plot_circos.R")


circos_data <- readRDS("src/circos/circos_data.rds")

circos <- prepare_circos(
    breakpoints = circos_data$myc_zoomed,
    genes = circos_data$genes_ig_myc_zoomed, # Coordinates of genes you want plotted with the colour of the rectangle in the `colour` column. All columns are required.
    region_order = circos_data$region_order, # Order of the regions as you want them plotted.
    partner_upstream_buffer = circos_data$buffer_zoomed, # How much space there should be up- and downstream of the genes at the ends of the plotting regions.
    partner_downstream_buffer = circos_data$buffer_zoomed,
    target_downstream_buffer = circos_data$buffer_pvt1,
    maf = maf_motif
)
pdf("test_circos.pdf")
plot_circos(circos,
    prop_chroms = circos_data$proportions_zoomed, # Proportions of the regions
    angle = 178,
    nest_colours = circos_data$nest_colours, # Named vector of colours for the shape that connects the region to the ideogram
    link_colours = circos_data$link_colours, # Named vector of colours for the links between regions. There should be a colour for every unique item in the column specified in the next argument.
    colour_links_by = "IG_feature", # The name of the column in the breakpoint table that should be used to assign colours.
    title = "IGH-MYC Breakpoints"
)
dev.off()
