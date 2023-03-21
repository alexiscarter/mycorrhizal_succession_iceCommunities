## Script of all the analyses

## Load libraries
library(tidyverse)
library(phyloseq)
library(psadd)
library(brms)
library(performance)
library(gdm)
library(ncf)

## Load data
full.table <- read.csv("data/full.table.09.02.23.csv")

# Transformations
full.table$sper.q0.s <- scale(full.table$sper.q0)
full.table$sper.q1.l <- log(full.table$sper.q1)
full.table$sper.q1.s <- scale(full.table$sper.q1)
full.table$am.q1.l <- log(full.table$am.q1)
full.table$ecm.q1.l <- log(full.table$ecm.q1)
full.table$t.s <- scale(full.table$meanT)
full.table$twi.l <- log(full.table$twi)
full.table$twi.s <- scale(full.table$twi.l)
full.table$ndvi.l <- log(full.table$ndvi)
full.table$ndvi.s <- scale(full.table$ndvi.l)

## AM-EcM differences
full.table$diff.sh <- full.table$am.q1.l-full.table$ecm.q1.l

## Modeling differences ####
diff.sh.g <- brm(formula = diff.sh ~ time.log.sc + (1|Glacier/Year), 
                 data = full.table, family=gaussian(), warmup = 1000, iter = 10000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Check for spatial autocorrelation
cres = spline.correlog(x = full.table$lat, y = full.table$lon, z = resid(diff.sh.g), resamp = 100)
plot(cres, ylim = c(-1,1)) # No sign spatial autocorrelation

## GDM ####
## Only account for plots where soil data is available
env.plot <- full.table %>%
  group_by(uniqPlot) %>%
  dplyr::select(uniqPlot, lon, lat, time.log.sc, t.s, twi.s, ndvi.s, ph.s, n.s, p.s) %>% 
  summarise_all(mean) %>% 
  mutate(Glacier = substr(uniqPlot, start = 1, stop = 5)) %>% 
  drop_na()

## Plant
sper <- readRDS("data/sper.phy.relax.rds")
sample_names(sper) <- paste(sample_data(sper)$Glacier, sample_data(sper)$Year, sample_data(sper)$Spot, sep = "_")
sample_data(sper)$uniqPlot <- substr(sample_names(sper), start = 1, stop = 12) # ## Summarize sample that have the same 12 first characters (for example, amola_1850_a1 and amola_1850_a2)
sper = merge_samples(sper, "uniqPlot")
sper <- prune_samples(sample_sums(sper) > 0, sper)
sper <- subset_samples(sper, sample_names(sper) %in% env.plot$uniqPlot) # 581 samples

## Plant community
sper.pcoa <- ordinate(sper, method="PCoA", distance = "jaccard")
sper.pcoa.load <- data.frame(compo.axis1 = sper.pcoa$vectors[,1],
                             compo.axis2 = sper.pcoa$vectors[,2])
sper.pcoa.load$uniqPlot <- rownames(sper.pcoa.load)

## Add plant community composition and exclude plots where no soil variables were not measured and where no MOTU found in sper.
env.plant.plot <- sper.pcoa.load %>% 
  left_join(env.plot, by = "uniqPlot")

## EcM at the plot level
## Remove samples with no taxa
fung.ecm <- prune_samples(sample_sums(fung.ecm) > 0, fung.ecm) # 429 samples

## Only keep plots with measured environmental variables
fung.ecm <- subset_samples(fung.ecm, sample_names(fung.ecm) %in% env.plant.plot$uniqPlot)

## Extract OTU table
ecm.mat = as.data.frame(as(otu_table(fung.ecm), "matrix"))
ecm.mat$uniqPlot <- rownames(ecm.mat)

## Format for GDM
gdmTab.ecm <- formatsitepair(bioData=ecm.mat, bioFormat=1, XColumn="lon", YColumn="lat",
                             predData=env.plant.plot, siteColumn="uniqPlot", dist = "jaccard")

## Remove comparison across glaciers
gdmTab.ecm <- gdmTab.ecm[-which(gdmTab.ecm$s1.Glacier != gdmTab.ecm$s2.Glacier),]

## Remove glacier name (for calculation)
gdmTab.ecm = gdmTab.ecm[,!(names(gdmTab.ecm) %in% c("s1.Glacier","s2.Glacier"))]

## GDM
gdm.ecm <- gdm(data=gdmTab.ecm, geo=T)
summary(gdm.ecm) # 15.248 deviance explained

## Significance testing with all variables (permutation, influence only pvalue estimation)
gdm.ecm.signif <- gdm.varImp(spTable = gdmTab.ecm, geo=T, fullModelOnly = TRUE, nPerm = 1000, parallel = TRUE)

## Extract deviance and pvalues
imp.ecm <- data.frame(deviance = gdm.ecm.signif[[2]][,1], pvalue = gdm.ecm.signif[[3]][,1], 
                      marker = rep('EcM',length(gdm.ecm.signif[[2]][,1])))
imp.ecm$variables <- rownames(imp.ecm)

## AM at the plot level
## Remove samples with no taxa
fung.am <- prune_samples(sample_sums(fung.am) > 0, fung.am) # 523 samples

## Only keep plots with measured environmental variables
fung.am <- subset_samples(fung.am, sample_names(fung.am) %in% env.plant.plot$uniqPlot) 

## Extract OTU table
am.mat = as.data.frame(as(otu_table(fung.am), "matrix"))
am.mat$uniqPlot <- rownames(am.mat)

## Format for GDM
gdmTab.am <- formatsitepair(bioData=am.mat, bioFormat=1, XColumn="lon", YColumn="lat",
                            predData=env.plant.plot, siteColumn="uniqPlot", dist = "jaccard")

## Remove comparison across glaciers
gdmTab.am <- gdmTab.am[-which(gdmTab.am$s1.Glacier != gdmTab.am$s2.Glacier),]

## Remove glacier name (for calculation)
gdmTab.am = gdmTab.am[,!(names(gdmTab.am) %in% c("s1.Glacier","s2.Glacier"))]

## GDM
gdm.am <- gdm(data=gdmTab.am, geo=T)
summary(gdm.am) # 17.66% deviance explained

## Significance testing with all variables (permutation, influence only pvalue estimation)
gdm.am.signif <- gdm.varImp(spTable = gdmTab.am, geo=T, fullModelOnly = TRUE, nPerm = 1000, parallel = TRUE)

## Extract deviance and pvalues
imp.am <- data.frame(deviance = gdm.am.signif[[2]][,1], pvalue = gdm.am.signif[[3]][,1], 
                     marker = rep('AM',length(gdm.am.signif[[2]][,1])))
imp.am$variables <- rownames(imp.am)

## Make global table
imp.all <- rbind(imp.ecm, imp.am)

imp.all %>%
  group_by(marker) %>% 
  summarize(d=sum(deviance))

## Random forest models ####
library(rfPermute)
library(randomForest)

## Only account for plots where soil data is available
full.r <- full.table %>%
  group_by(uniqPlot) %>%
  dplyr::select(uniqPlot, am.q1.l, ecm.q1.l, sper.q1.l,ndvi.l,time.log,meanT,twi.l,lg_n,ph,lg_p) %>% 
  summarise_all(mean) %>% 
  mutate(Glacier = substr(uniqPlot, start = 1, stop = 5)) %>% 
  drop_na()

## Final models
## RF AM
rf_am=randomForest(am.q1.l~sper.q1.l+ndvi.l+time.log+meanT+twi.l+lg_n+ph+lg_p+Glacier, data=full.r, 
                   importance=T, 
                   localImp=T,
                   mtry = 2,
                   ntree = 600)
summary(rf_am) #0.237346 46.99

rp_am <- rfPermute(am.q1.l~sper.q1.l+ndvi.l+time.log+meanT+twi.l+lg_n+ph+lg_p+Glacier, data=full.r, 
                   ntree = 600, 
                   num.rep = 5000, 
                   mtry = 2, 
                   num.cores = NULL)
imp.am <- importance(rp_am)
imp.am <- data.frame(variables = row.names(imp.am), imp.am, myco = rep('AM', 9))
imp.am <- imp.am[order(imp.am$variables),]

## RF EcM
rf_ecm=randomForest(ecm.q1.l~sper.q1.l+ndvi.l+time.log+meanT+twi.l+lg_n+ph+lg_p+Glacier, data=full.r, 
                    importance=T, 
                    localImp=T,
                    mtry = 2,
                    ntree = 600)
summary(rf_ecm) #0.1093399 48.96

rp_ecm <- rfPermute(ecm.q1.l~sper.q1.l+ndvi.l+time.log+meanT+twi.l+lg_n+ph+lg_p+Glacier, data=full.r, 
                    ntree = 600, 
                    num.rep = 5000, 
                    mtry = 2, 
                    num.cores = NULL)
imp.ecm <- importance(rp_ecm)
imp.ecm <- data.frame(variables = row.names(imp.ecm), imp.ecm, myco = rep('EcM', 9))
imp.ecm <- imp.ecm[order(imp.ecm$variables),]

## Merge results
rf <- rbind(imp.am, imp.ecm)
