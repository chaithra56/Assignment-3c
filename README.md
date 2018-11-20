# Assignment-3c
# PART C0: Biological Question + Hypothesis ----

#Do subfamilies (Viperinae and Crotalinae ) within the Viperidae family live in the same geographic region or different geographic regions?

#Viperidae is a large family of venomous snakes with 3 main sub-families (albeit 1 is very small). Within those sub-families is Viperinae and Crotalinae that are sister groups and are distinguished by the lack of (Viperinae) or presence of (Crotalinae) a heat-sensing pit organ. Given the distinguished morphological trait, I believe the sister groups would not share geographical regions as by evolutionary theory, Crotalinae developed the need to have a heat-sensing pit organ for a reason.

# PART C1: LOAD PACKAGES REQUIRED FOR PART C----

install.packages("rworldmap")
library(rworldmap)
library(plyr)
library(seqinr)
#Install packages from Bioconductor
source("https://bioconductor.org/biocLite.R")
library(Biostrings)
library(stringr)
library(tidyverse)
library(ape)
library(seqinr)
library(muscle)
library(msa)
library(DECIPHER)
library(RSQLite)
library(kmer)
library(ggdendro)
library(ggplot2)
library(cluster)
library(dendextend)
library(dplyr)

# PART C2: OBTAIN DATA FROM BOLD TO BE USED FOR PHYLOGENETIC ANALYSIS ----

Viperidae <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Viperidae&format=tsv")

write_tsv(Viperidae, "Viperidae_BOLD_data.tsv")

viper <- read_tsv("Viperidae_BOLD_data.tsv")

rm(Viperidae)
#the data of Viperidae and viper are same. Viper data is used to perform other functions in the program so the data of Viperidae can be eliminated to save the storage space.


#Read data from BOLD and created a file in case a duplicate of the original file is required. 
# PART C3: FILTER AND TIDY UP DATA TO BE USED FOR PHYLOGENETIC ANALYSIS ----

Sequence.Viperidae <- viper %>%
  select(subfamily_name, country, nucleotides) %>%
  drop_na()
#Created a dataframe that has the sub-family names, country, and nucleotides.

Sequence.Viperinae <- Sequence.Viperidae %>%
  filter(str_detect(subfamily_name, "Viperinae")) %>%
  filter(str_detect(nucleotides, "AGTC"))
#Created a dataframe that has the Viperinae sequences only.

Sequence.Crotalinae <- Sequence.Viperidae %>%
  filter(str_detect(subfamily_name, "Crotalinae")) %>%
  filter(str_detect(nucleotides, "AGTC"))
#Created a dataframe that has the Crotalinae sequences only.

# to perform some functions like alignments or to plot a data the data shouldn't be null. In order to eliminate the empty data the function na.omit can be used. For example :
length(viper$nucleotides)
Fulltest <- na.omit(viper$nucleotides)
length(Fulltest)



# PART C4: RANDOMLY SAMPLE SEQUENCES PER COUNTRY ----

RS.Viperinae <- Sequence.Viperinae %>%
  group_by(country) %>%
  sample_n(1)
#Randomly sample one sequence per country 

RS.Crotalinae <- Sequence.Crotalinae %>%
  group_by(country) %>%
  sample_n(1)
#Randomly sample one sequence per country 

# PART C5: CONDUCT A MULTIPLE SEQUENCE ALIGNMENT ON EACH FAMILY %>% ----

class(RS.Viperinae$nucleotides)
class(RS.Crotalinae$nucleotides)
#Convert 'character' into DNAStringSet.

RS.Viperinae$nucleotides <- DNAStringSet(RS.Viperinae$nucleotides)
RS.Crotalinae$nucleotides <- DNAStringSet(RS.Crotalinae$nucleotides)
#Converted into a DNAStringSet object that can be used for a MUSCLE analysis.

Viperinae.Sequence <- RS.Viperinae$nucleotides
Crotalinae.Sequence <- RS.Crotalinae$nucleotides
#Create a DNAStringSet object that contains the nucleotie sequences to be aligned.

Viperinae.MUSCLE <- DNAStringSet(muscle::muscle(Viperinae.Sequence, maxiters = 2, diags = TRUE))
Crotalinae.MUSCLE <- DNAStringSet(muscle::muscle(Crotalinae.Sequence, maxiters = 2, diags = TRUE))
#Conduct a Multiple Sequence Alignment (MSA) of the randomly sampled sequences from each country for each subfamily using MUSCLE

# PART C6: CONDUCT A CLUSTER ANALYSIS FOR EACH SUB-FAMILY TO CREATE A DENDROGRAM----

class(Viperinae.MUSCLE)
#Need to convert DNAStringSet into DNAbin

dnaBin.Viperinae <- as.DNAbin(Viperinae.MUSCLE)
dnaBin.Crotalinae <- as.DNAbin(Crotalinae.MUSCLE)
#Converted DNAStringSet into DNAbin.

Viperinae.distMatrix <- dist.dna(dnaBin.Viperinae, model = "TN93", as.matrix = TRUE, 
                                 pairwise.deletion = TRUE)
Crotalinae.distMatrix <- dist.dna(dnaBin.Crotalinae, model = "TN93", as.matrix = TRUE, 
                                  pairwise.deletion = TRUE)
#Create a distance matrix using the alignment from muscle.

Viperinae.Clusters <- IdClusters(Viperinae.distMatrix,
                                 method = "single",
                                 cutoff= 0.02,
                                 showPlot = TRUE,
                                 type = "both",
                                 verbose = TRUE)
Crotalinae.Clusters <- IdClusters(Crotalinae.distMatrix,
                                  method = "single",
                                  cutoff= 0.02,
                                  showPlot = TRUE,
                                  type = "both",
                                  verbose = TRUE)
#Created clusters of both sub-families.

Viperinae.Dendrogram <- IdClusters(Viperinae.distMatrix,
                                   method = "single",
                                   cutoff= 0.02,
                                   showPlot = TRUE,
                                   type = "dendrogram",
                                   verbose = TRUE)
Crotalinae.Dendrogram <- IdClusters(Crotalinae.distMatrix,
                                    method = "single",
                                    cutoff= 0.02,
                                    showPlot = TRUE,
                                    type = "dendrogram",
                                    verbose = TRUE)
#Created dendrograms for both sub-families of each country. 

# PART C7: CREATING THE GEOGRAPHICAL MAP OF SPECIES----

Test <- viper %>%
  select(subfamily_name, country) %>%
  drop_na()


#Creates a dataframe that includes only the sub-family names and the corresponding country.

match <- joinCountryData2Map(Test, joinCode="NAME", nameJoinColumn="country")
#This will match the countries that are in the Test dataframe to be used on the world map.

mapParams <- mapCountryData(match, nameColumnToPlot = "subfamily_name", mapRegion = "world", catMethod = "categorical",colourPalette = "heat", addLegend = FALSE, borderCol = "grey",mapTitle = "Geographical Distribution of Crotalinae and Viperinae", oceanCol = 'lightblue', aspect = 1,missingCountryCol = "white", add = FALSE, nameColumnToHatch = "", lwd = 0.5)
#Generates a world map that includes the distribution of the Viperinae and Crotalinae sub-families across the world. 

do.call( addMapLegendBoxes, c(mapParams,x='bottom',horiz=TRUE,title="Sub-Families"))
#Customize the legend title so that everything looks seamless

#Viperidae has sub-families Viperinae and Crotalinae that are sister groups and are distinguished by the lack of (Viperinae) or presence of (Crotalinae) a heat-sensing pit organ. These subfamilies are represented in world map all over the countries.The map represents the subfamily Crotoalinae is mostly present in North America and approximately one-third of the South america and Eastern part of Asia and northern part of Asia. Viperinae subfamily mostly present in India, Western part of Europe and few countries in Africa continent.As per the hypothesis both the subfamilies doesn't share same region and live totally different climatic conditions.
