## PACKAGE INSTALL
# List of required packages
packages <- c("dplyr", "tidyr", "phyloseq", "qiime2R", "ggplot2", "vegan", "fantaxtic", "ggpubr", "decontam")

# Check if packages are not already installed and install them
install_packages <- packages[!packages %in% installed.packages()]
if (length(install_packages) > 0) {
  install.packages(install_packages)
}

# Additional installation for 'decontam' package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("decontam", quietly = TRUE)) {
  BiocManager::install("decontam")
}

# libraries
# use code above if you have not run this code before!
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
library("vegan")
library("fantaxtic")
library("ggpubr")
library("decontam")

##################
# phyloseq setup #
##################

# feature table
ASV <- qza_to_phyloseq(features="required_files/table.qza") # set ASV (feature) table for phyloseq
# use for later
metatable <- read.delim("required_files/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names="sample_name") 
metatable$is.neg <- metatable$Sample.Control == "Control.Sample" # making a boolean column to coordinate with controls (TRUE/FALSE) (for decontam package)
metatable$Final_Qubit <- as.numeric(metatable$Final_Qubit) # make the column numeric (for decontam package)
#metatable <- filter(metatable, Sample.Control == "True.Sample")# filter by transect
META <- sample_data(metatable) # set metatable for phyloseq

# importing taxonomy
taxonomy <- read.delim("required_files/taxonomy.tsv", sep="\t", header=TRUE) # set taxonomy for phyloseq
names(taxonomy) <- c("row", "tax", "Confidence") # rename columns # CONFIDENCE LEVEL .7 in the QIIME2 script
row.names(taxonomy) <- taxonomy[[1]] # making the row names the tax name [[1]] are row.names
taxonomy <- taxonomy[,(-1)] # removing the duplicate column
# SILVA taxonomy is in one column, separate to be able to work with different taxonomic levels:
taxonomy <-  separate(taxonomy, tax, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", 
                                       "D7", "D8", "D9", "D10", "D11", 
                                       "D12", "D13", "D14"), sep = ";", fill = "right") # creates columns for each taxonomy
taxonomy <- taxonomy[,c(1:7)] # select only first 1-7 columns (to genus)
nm1 <- colnames(taxonomy) # setting a vector column names (sample names)
taxonomy[nm1] <- lapply(taxonomy[nm1], gsub, pattern = "D_.__", replacement = "") # cycling through the taxonomy to remove D__2_ before taxonomy names.
taxonomy <- taxonomy[-1, ] # removes first row of taxonomy table
names <- row.names(taxonomy)
row.names(ASV) <- names # matching the names from taxonomy to otu table
taxmat <- as.matrix(taxonomy) # changing to matrix for phyloseq setup
TAX <- tax_table(taxmat) # set tax table for phyloseq

# import rooted tree
TREE <- qza_to_phyloseq(tree="required_files/rooted-tree.qza") # setting tree for phyloseq object
TREE[["tip.label"]] <- names # changing the labels to the same taxonomy names (vector from before)

# merge all imported objects into phyloseq
ps <- merge_phyloseq(ASV, TAX, META, TREE) # merge all imported objects into semi-final phyloseq
ps <- subset_samples(ps, sample.illumina != "078_1040_PRE") # removed bc qubit value is 0
# and we are done!!!!!! but! we need to look through to make sure we dont have any contamination.
# there is the package "decontam", which looks at the DNA concentration before PCR and the frequency+prevelance of taxa in the kitblanks

ps_decontam <- isContaminant(ps, conc="Final_Qubit", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined") #decontam, using phyloseq object (ps), setting method to "combined", that uses frequency/prev, with 0.5 threshold
table(ps_decontam$contaminant) # it identified 133 potential contaminants
ps_decontam_list <- tibble::rownames_to_column(ps_decontam, var = "ASV") #using tibble package, might have to install?
#write.csv(ps_decontam_list, "ps_decontam_list.csv") # can export as csv if needed

ps_d <- prune_taxa(!ps_decontam$contaminant, ps) # taking out taxa based if the contaminant is TRUE, thats the ! before the dataframe
ps_d #  4559 taxa and 279 samples (originally was 4692 taxa)


