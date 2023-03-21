## Script to access and process the data ##

## Load libraries
library(phyloseq)   # Manipulation of metabarcoding data  
library(tidyverse)  # Plotting and data manipulation 
library(stringi)

## Files available
list.files(path = "data/")
#data_filtered_Fung02.zip: Filtered sequence data from Fung02 marker in compressed format.
#data_filtered_Sper01.zip: Filtered sequence data from Sper01 marker in compressed format.
#env.myco.09.02.23.csv: Environmental data containing plot information, time since glacier retreat, soil chemistry (pH, C, N, P), temperature (meanT), productivity (NDVI), wetness (WTI) and plot coordinates. 
#fung.phy.relax.rds: Phyloseq object of the fungal data (can be created using this script).
#sper.phy.relax.rds: Phyloseq object of the plant data (can be created using this script).
#full.table.09.02.23.csv: Table containing calculated diversity values with corresponding and environmental values (can be created using this script).

## Load data
## Sequences (from zipped files)
fung <- read.csv(unz("data/data_filtered_Fung02.zip","data_filtered_Fung02.csv"))
sper <- read.csv(unz("data/data_filtered_Sper01.zip","data_filtered_Sper01.csv"))

## Sequence data ####
## Fungi
## Select taxonomy table
fung.tax <- fung[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
fung.tax <-  fung.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
fung.motu <- fung[,-c(2:25)]
## Manipulation of OTU table
fung.motu <- fung.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
fung.motu.relax <- fung.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
fung.phy.relax = phyloseq(otu_table(fung.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(fung.tax),
                          sample_data(fung.motu.relax[,1:3]))
fung.phy.relax <- prune_taxa(taxa_sums(fung.phy.relax) > 0, fung.phy.relax) 
#saveRDS(fung.phy.relax, file = "data/fung.phy.relax.rds")

## Plants
## Select taxonomy table
sper.tax <- sper[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
sper.tax <-  sper.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
sper.motu <- sper[,-c(2:25)]
## Manipulation of OTU table
sper.motu <- sper.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
sper.motu.relax <- sper.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
sper.phy.relax = phyloseq(otu_table(sper.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(sper.tax),
                          sample_data(sper.motu.relax[,1:3]))
sper.phy.relax <- prune_taxa(taxa_sums(sper.phy.relax) > 0, sper.phy.relax) 
#saveRDS(sper.phy.relax, file = "data/sper.phy.relax.rds")

## Calculate diversity ####
## fung
fung <- readRDS("data/fung.phy.relax.rds")
sample_names(fung) <- paste(sample_data(fung)$Glacier, sample_data(fung)$Year, sample_data(fung)$Spot, sep = "_")
## Summarize sample that have the same 12 first characters (for example, amola_1850_a1 and amola_1850_a2)
sample_data(fung)$uniqPlot <- substr(sample_names(fung), start = 1, stop = 12)
fung = merge_samples(fung, "uniqPlot")
## Diversity values
fung.div <- estimate_richness(fung, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%
  mutate(q1 = exp(Shannon))

## sper
sper <- readRDS("data/sper.phy.relax.rds") # At the plot-level
sample_names(sper) <- paste(sample_data(sper)$Glacier, sample_data(sper)$Year, sample_data(sper)$Spot, sep = "_")
## Summarize sample that have the same 12 first characters (for example, amola_1850_a1 and amola_1850_a2)
sample_data(sper)$uniqPlot <- substr(sample_names(sper), start = 1, stop = 12)
sper = merge_samples(sper, "uniqPlot")
## Diversity values
sper.div <- estimate_richness(sper, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%
  mutate(q1 = exp(Shannon))

## Diversity table
div.table <- data.frame(Glacier = sper.div$Glacier, Year = as.numeric(sper.div$Year), Plot = sper.div$Plot, 
                        sper.q1 = sper.div$q1, fung.q1 = fung.div$q1,
                        sper.q0 = sper.div$Observed, fung.q0 = fung.div$Observed)

## Number of glaciers
length(unique(div.table$Glacier))
## Number of plots
phyloseq::nsamples(fung)
## Number of taxa
ntaxa(fung)
ntaxa(sper)
## Number of MOTUs
sum(sample_sums(fung))
ntaxa(sper)

## Select EcM fungal taxon
fung <- readRDS("data/fung.phy.relax.rds")
fung.ecm = subset_taxa(fung, genus_name == "Inocybe"|genus_name == "Austropaxillus"|genus_name =="Cantharellus"|   genus_name =="Cenococcum"|genus_name =="Clavulina"|family_name =="Cortinariaceae"|family_name =="Gomphidiaceae"|genus_name =="Helvella"|genus_name =="Lactarius"|genus_name =="Leucophleps"|genus_name =="Rhizopogon"|genus_name =="Russula"|family_name =="Sebacinaceae"|genus_name =="Suillus"|family_name =="Tuberaceae")
## Calculate diversity
sample_names(fung.ecm) <- paste(sample_data(fung.ecm)$Glacier, sample_data(fung.ecm)$Year, sample_data(fung.ecm)$Spot, sep = "_")
sample_data(fung.ecm)$uniqPlot <- substr(sample_names(fung.ecm), start = 1, stop = 12)
fung.ecm = merge_samples(fung.ecm, "uniqPlot")
fung.ecm.div <- estimate_richness(fung.ecm, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## Select AM fungal taxon
fung <- readRDS("data/fung.phy.relax.rds")
fung.am = subset_taxa(fung, family_name == "Acaulosporaceae"|family_name == "Archaeosporaceae"|order_name == "Archaeosporales"|family_name == "Diversisporaceae"|order_name == "Diversisporales"|family_name == "Glomeraceae"|order_name == "Glomerales"|family_name == "Paraglomeraceae")
## Calculate diversity
sample_names(fung.am) <- paste(sample_data(fung.am)$Glacier, sample_data(fung.am)$Year, sample_data(fung.am)$Spot, sep = "_")
sample_data(fung.am)$uniqPlot <- substr(sample_names(fung.am), start = 1, stop = 12)
fung.am = merge_samples(fung.am, "uniqPlot")
fung.am.div <- estimate_richness(fung.am, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%
  mutate(q1 = exp(Shannon))

div.myco.table <- data.frame(Glacier = fung.ecm.div$Glacier, Year = as.numeric(fung.ecm.div$Year), Plot = fung.ecm.div$Plot, 
                             ecm.q0 = fung.ecm.div$Observed, ecm.q1 = fung.ecm.div$q1,
                             am.q0 = fung.am.div$Observed, am.q1 = fung.am.div$q1)

## Environmental data ####
env.table <- read.csv("data/env.myco.09.02.23.csv") #PS: based on full div.clim.chem.21.03.2023.csv

## Join environmental and diversity data
full.table <- env.table %>% 
  left_join(div.myco.table, by = c("Glacier", "Year", "Plot")) %>%
  left_join(div.table, by = c("Glacier", "Year", "Plot"))
#write.csv(full.table, "data/full.table.09.02.23.csv", row.names = F)

# Transformations
full.table$sper.q0.s <- scale(full.table$sper.q0)
full.table$sper.q1.s <- scale(full.table$sper.q1)
full.table$sper.q1.l <- log(full.table$sper.q1)
full.table$fung.q1.l <- log(full.table$fung.q1)
full.table$am.q1.l <- log(full.table$am.q1)
full.table$ecm.q1.l <- log(full.table$ecm.q1)
full.table$t.s <- scale(full.table$meanT)
full.table$twi.s <- scale(log(full.table$twi))
full.table$ndvi.s <- scale(log(full.table$ndvi))
full.table$ph2 <- full.table$ph.s^2

