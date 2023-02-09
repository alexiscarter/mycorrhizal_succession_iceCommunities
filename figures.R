## Script to reproduce all the figures ##
## Before plotting, need to run import.R and/or analyses.R

## Load libraries
library(tidyverse)
library(phyloseq)

# colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

## Plot diversity
## All fungi
ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=fung.q0))+ # +0.0000001 is to avoid small misplacement of the bins
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "Fungal richness") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 483), limits = c(0.95,NA)) +
  theme_bw()

ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=fung.q1.l))+
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "Fungal diversity (q=1)") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 483), limits = c(0.95,NA)) +
  theme_bw()

## Plants
ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=sper.q0))+
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "Vascular plant diversity (q=0)") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 483), limits = c(0.95,NA)) +
  theme_bw()

ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=sper.q1.l)) +
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "Vascular plant diversity (q=1)") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 483), limits = c(0.95,NA)) +
  theme_bw()

## EcM diversity
ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=ecm.q1.l))+
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "EcM fungal diversity (q=1)") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 500), limits = c(0.95,NA)) +
  scale_y_continuous(limits = c(-.01,2.92))+
  theme_bw()

ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=ecm.q0))+
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "EcM fungal richness") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 483), limits = c(0.95,NA)) +
  theme_bw()

## AM diveristy
ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=am.q1.l))+ 
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "AM fungal diversity (q=1)") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 500), limits = c(0.95,NA)) +
  scale_y_continuous(limits = c(-.01,2.92))+
  theme_bw()

ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=am.q0))+
  geom_hex(bins = 12) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  labs(title="", x = "Time after glacier retreat (years)", y = "AM fungal richness") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 483), limits = c(0.95,NA)) +
  theme_bw()

## Plot EcM MOTUs - family
ecm.melt <- psmelt(fung.ecm)
ecm.melt$Abundance[ecm.melt$Abundance > 0] <- 1 

ecm.melt.div <- ecm.melt %>% 
  separate(col = Sample, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%
  mutate(Year = as.numeric(Year)) %>% 
  left_join(full.table[,c(1,2,3,6)], by = c("Glacier", "Year", "Plot")) %>% ## PLOT LEVEL
  group_by(Glacier, time, Plot, family_name) %>%
  summarise(q0 = sum(Abundance))

ecm.melt.div$class <- '1-16'
ecm.melt.div$class[ecm.melt.div$time>16]<-'17-36'
ecm.melt.div$class[ecm.melt.div$time>36]<-'37-74'
ecm.melt.div$class[ecm.melt.div$time>74]<-'75-144'
ecm.melt.div$class[ecm.melt.div$time>144]<-'145-483'
ecm.melt.div$class <- factor(ecm.melt.div$class, levels = c("1-16", "17-36", "37-74", "75-144", "145-483"))

ggplot(ecm.melt.div, aes(x=class, y=q0, fill=family_name)) +
  geom_bar(stat = "summary") +
  scale_fill_manual(name = "EcM families", values = col_vector) +
  labs(x = 'Time after glacier retreat (years)', y = 'EcM fungal richness (q=0)') +
  theme_bw()

## AM MOTUs - family
am.melt <- psmelt(fung.am)
am.melt$Abundance[am.melt$Abundance > 0] <- 1 

am.melt.div <- am.melt %>% 
  separate(col = Sample, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%
  mutate(Year = as.numeric(Year)) %>% 
  left_join(full.table[,c(1,2,3,6)], by = c("Glacier", "Year", "Plot")) %>% 
  group_by(Glacier, time, Plot, family_name) %>% 
  summarise(q0 = sum(Abundance))

am.melt.div$class <- '1-16'
am.melt.div$class[am.melt.div$time>16]<-'17-36'
am.melt.div$class[am.melt.div$time>36]<-'37-74'
am.melt.div$class[am.melt.div$time>74]<-'75-144'
am.melt.div$class[am.melt.div$time>144]<-'145-483'
am.melt.div$class <- factor(am.melt.div$class, levels = c("1-16", "17-36", "37-74", "75-144", "145-483"))

ggplot(am.melt.div, aes(x=class, y=q0, fill=family_name)) +
  geom_bar(stat = "summary") +
  scale_fill_manual(name = "AM families", values = col_vector) +
  labs(x = 'Time after glacier retreat (years)', y = 'AM fungal richness (q=0)') +
  theme_bw()

##Plot AM-EcM differences
full.table$time.s <- scale(full.table$time)
full.table$time.log.sc <- scale(log(full.table$time)) 
marg <- conditional_effects(diff.sh.g)
marg.dat <- marg$time.log.sc
ggplot(data = full.table, aes(x=exp(time.log)+0.0000001, y=diff.sh))+ # +0.0000001 is to avoid small misplacement of the bins
  geom_hex(bins = 15) +
  scale_fill_viridis_c(name = "Number\nof plots", option = "F", direction = -1, trans="log10", n.breaks = 6, limits = c(1,140), alpha = .5) +
  geom_point(alpha=0.2, shape = 16) +
  geom_line(data = marg.dat, aes(x=exp(time.log.sc*attr(full.table$time.log.sc, 'scaled:scale') + attr(full.table$time.log.sc, 'scaled:center')), y=estimate__) ,size=1, alpha = 1) +
  geom_ribbon(data = marg.dat, aes(x=exp(time.log.sc*attr(full.table$time.log.sc, 'scaled:scale') + attr(full.table$time.log.sc, 'scaled:center')), y=estimate__, ymax=upper__, ymin=lower__), fill="lightblue", alpha=0.5) +
  labs(title="", x = "Time after glacier retreat (years)", y = "EcM diversity - AM diversity") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 500), limits = c(0.95,NA)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(-3,3)) + #, trans = "log10", limits = c(.9,48)
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black",size = .8), axis.ticks = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text = element_text(family = 'Helvetica', colour = "black"), text = element_text(family = 'Helvetica', colour = "black"), plot.title = element_text(hjust=0.5))


## Random Forest ####
# package not working on my machine yet, ran on PC
rf$variables <- ordered(rf$variables, levels = c("Glacier", "Time", "NDVI", "Plant diversity", "Nitrogen", "pH", "Phosphorus", "Temperature", "Wetness"))
rf %>% 
  filter(myco == "AM") %>% 
  ggplot(aes(x=variables, y = pourc_IncMSE, fill = variables))+
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Variables", values = c('#DDDDDD', '#FFAABB', '#BBCC33',  '#AAAA00', '#00a878', '#EE8866', '#EEDD88', '#99DDFF', '#77AADD'),
                    labels = c("Glacier", "Time", "Productivity", "Plant diversity", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  scale_x_discrete(labels = c("Glacier", "Time", "Productivity", "Plant diversity", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  labs(x = 'Drivers of alpha-diversity', y = 'Variable importance (% IncMSE)', title = 'AM fungi') +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,60), n.breaks = 6) + #, trans = "log10", limits = c(0,60)
  theme_bw() + theme(axis.text.x = element_text(angle = -70, vjust = 0.5, hjust=0), plot.title = element_text(hjust = 0.5))

rf %>% 
  filter(myco == "EcM") %>% 
  ggplot(aes(x=variables, y = pourc_IncMSE, fill = variables))+
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Variables", values = c('#DDDDDD', '#FFAABB', '#BBCC33',  '#AAAA00', '#00a878', '#EE8866', '#EEDD88', '#99DDFF', '#77AADD'),
                    labels = c("Glacier", "Time", "Productivity", "Plant diversity", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  scale_x_discrete(labels = c("Glacier", "Time", "Productivity", "Plant diversity", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  labs(x = 'Drivers of alpha-diversity', y = 'Variable importance (% IncMSE)', title = 'EcM fungi') +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,60), n.breaks = 6) + #, trans = "log10", limits = c(0,60)
  theme_bw() + theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0) plot.title = element_text(hjust = 0.5))

## PCoA plots
