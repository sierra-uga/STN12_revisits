## PACKAGE INSTALL
# List of required packages
packages <- c("dplyr", "tidyr", "phyloseq", "ggplot2", "plyr", "fantaxtic", "ggpubr", "RColorBrewer")
# Check if packages are not already installed and install them
install_packages <- packages[!packages %in% installed.packages()]
if (length(install_packages) > 0) {
  install.packages(install_packages)
}

# load libraries
# use code above if you have not run this code before!
library("dplyr") # part of tidyverse - data manipulation
library("tidyr") # part of tidyverse - data manipulation
library("phyloseq") # phylogenic analysis package
library("ggplot2") # for plotting
library("plyr") # part of tidyverse - data manipulation
library("fantaxtic") # for arranging plots together
library("ggpubr") # for arranging plots together
library("RColorBrewer") # for getting colors for plots

##################
#  data culling  #
##################

# select only bacteria, remove chloroplasts
ps_sub <- ps_d %>% # the %>% is part of the tidyverse (tidyr/dyplr), it acts like a transition, using the same object without having to rename it
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

wanted_stations <- c("STN012", "STN115", "STN12.3") # make a vector of desired stations
ps_sub <- ps_sub %>% subset_samples(Station %in% wanted_stations) # filter phyloseq object by the vector

# free-living phyloseq, filter by pore size, then remove the taxa that have 0
ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq, filter by pore size, then remove the taxa that have 0
ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 


# Create a data frame for freeliving, agglomerate by Order, transform to rel.abundance
data_free <- ps_free %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at Order level, can change to different taxonomic level!
  transform_sample_counts(function(x) {x/sum(x)})# Transform to rel. abundance (normalize data

top_free <- top_taxa(data_free, # selecting the top 12 taxa from the data_free phyloseq object
                      n_taxa = 12, # if you change this number, you will have to add/remove colors from the myColors (see below)
                      include_na_taxa = T)

data_top_free <- top_free$ps_obj %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>%  # Filter out low abundance taxa
  arrange(Order)# arrange by Order

# particle-associated
data_part <- ps_part %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at Order level
  transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance (normalize data)

top_part <- top_taxa(data_part, # selecting the top 12 taxa from the data_free phyloseq object
                     n_taxa = 12, # if you change this number, you will have to add/remove colors from the myColors (see below)
                     include_na_taxa = T)

data_top_part <- top_part$ps_obj %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>%   # Filter out low abundance taxa
  arrange(Order)  # arrange by Order

###################
# Plot variables! #
###################

plot_labels <- unique(sample_data(ps_free)$Station) # plot labels for graph, each Station (3 total)
plot_breaks <- unique(sample_data(ps_free)$Station) # plot breaks for graph, each Station (3 total)
myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkgoldenrod", "lightblue2", "darkblue", "dodgerblue", "seagreen", "purple", "black") # this must equal the levels of the Order
data_top_free$Order <- as.factor(data_top_free$Order) # setting the Order columns to factor
data_top_part$Order <- as.factor(data_top_part$Order) # setting the Order columns to factor
names(myColors) <- levels(c(data_top_free$Order, data_top_part$Order)) # setting the names of the colors to coordinate with the Order columns of each dataframe

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

barplot_free <- ggplot(data_top_free, aes(x = Station, y = Abundance, fill = Order, group = Order)) + facet_grid(~factor(watertype, levels=c("AASW", "AASW-WW", "WW", "WW-CDW", "CDW"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = myColors, drop = FALSE) + # set the colors with custom colors (myColors)
  scale_x_discrete(
    breaks = plot_breaks, # setting breaks
    labels = plot_labels, # settting levels
    drop = FALSE
  ) +
  #theme(plot.title = element_text(hjust = 0.5, size=17)) + # remove # if you want title
  theme(axis.title.x = element_blank()) + # removing x-axis title
  theme(axis.text.x = element_text(size=9)) + # setting x-axis title
  theme(axis.title.y = element_blank()) + # removing y-axis title
  theme(legend.position = "none") + # remove legend, delete line if you want a legend
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")
ggsave("graphics/free_living_barplot.pdf", width = 8, height = 6, dpi = 150)

###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

# the following plot is basically the same as above, look at annotation for free-living barplot if confused about what each line does!
barplot_part <- ggplot(data_top_part, aes(x = Station, y = Abundance, fill = Order, group = Order)) + facet_grid(~factor(watertype, levels=c("AASW", "AASW-WW", "WW", "WW-CDW", "CDW"))~.) +
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop = FALSE) +
  scale_x_discrete(
    breaks = plot_breaks,
    labels = plot_labels,
    drop = FALSE
  ) +
  #theme(plot.title = element_text(hjust = 0.5, size=17)) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(size=9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size=9 , color="white")) + # makes color of y-axis text white so its even when combining plots, remove entire color if want back to normal.
  theme(legend.position = "none") +
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ggtitle("Particle-associated (>3 µm)")
ggsave("graphics/part_associated_barplot.pdf", width = 8, height = 6, dpi = 150)

###################### 
#  stacked barplots  #
######################
#  BOTH COMMUNITIES  # 
######################

total <- rbind(data_top_part, data_top_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Station, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) +
  guides(fill = guide_legend(override.aes = list(color = "black", size = 1))) # adds black outline around legend

legend_combined <- get_legend(legend_plot)

ps_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right", legend.grob = legend_combined
)
# to remove the white space between the two community plots, you'd have to play with the plot.margins of each plot individually!
# something like: plot.margin=unit(c(1,1,-0.5,1), "cm")), where the margins follow the following structure:
# unit(c(top, right, bottom, left), units).

annotate_figure(ps_combined, top = text_grob("Station 12 Time Series", 
                                                    color = "black", face = "bold", size = 18))

ggsave("graphics/combined_barplot.pdf", width = 13, height = 7, dpi = 150)