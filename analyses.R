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

# Transformation
full.table$sper.q0.s <- scale(full.table$sper.q0)
full.table$sper.q1.s <- scale(full.table$sper.q1)
full.table$am.q1.l <- log(full.table$am.q1)
full.table$ecm.q1.l <- log(full.table$ecm.q1)
full.table$t.s <- scale(full.table$meanT)
full.table$twi.s <- scale(log(full.table$twi))
full.table$ndvi.s <- scale(log(full.table$ndvi))
full.table$ph2 <- full.table$ph.s^2

## AM-EcM differences
full.table$diff.sh <- full.table$am.q1.l-full.table$ecm.q1.l

## Modeling differences ####
diff.sh.g <- brm(formula = diff.sh ~ time.log.sc + (1|Glacier/Year), 
                 data = full.table, family=gaussian(), warmup = 1000, iter = 10000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Check for spatial autocorrelation
cres = spline.correlog(x = full.table$lat, y = full.table$lon, z = resid(diff.sh.g), resamp = 100)
plot(cres, ylim = c(-1,1)) # No sign spatial autocorrelation

## PCoA

## GDM

## Random forest model

