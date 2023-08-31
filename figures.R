## Script to reproduce all the figures ##
## Before plotting, need to run import.R and/or analyses.R

## Load libraries
library(tidyverse)
library(phyloseq)

# colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

## Plot diversity ####
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

## Plot AM-EcM differences
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

## Plot probabilities of being present ####
## AM
marg <- conditional_effects(am.bern)
marg.dat <- marg$time.log.sc
am.bern.plot <- ggplot(data = full.table.pres, aes(x=exp(time.log)+0.0000001, y=pres.am))+ # +0.0000001 is to avoid small misplacement of the bins
  geom_line(data = marg.dat, aes(x=exp(time.log.sc*attr(full.table$time.log.sc, 'scaled:scale') + attr(full.table$time.log.sc, 'scaled:center')), y=estimate__) ,size=1, alpha = 1) +
  geom_ribbon(data = marg.dat, aes(x=exp(time.log.sc*attr(full.table$time.log.sc, 'scaled:scale') + attr(full.table$time.log.sc, 'scaled:center')), y=estimate__, ymax=upper__, ymin=lower__), fill="lightblue", alpha=0.5) +
  labs(title="", x = "Time after glacier retreat (years)", y = "Probability of AM fungi to be present\nin the overall fungal community") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 500), limits = c(0.95,NA)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,1)) + #, trans = "log10", limits = c(.9,48)
  theme_bw()

## EcM
marg <- conditional_effects(ecm.bern)
marg.dat <- marg$time.log.sc
ecm.bern.plot <- ggplot(data = full.table.pres, aes(x=exp(time.log)+0.0000001, y=pres.ecm))+ # +0.0000001 is to avoid small misplacement of the bins
  geom_line(data = marg.dat, aes(x=exp(time.log.sc*attr(full.table$time.log.sc, 'scaled:scale') + attr(full.table$time.log.sc, 'scaled:center')), y=estimate__) ,size=1, alpha = 1) +
  geom_ribbon(data = marg.dat, aes(x=exp(time.log.sc*attr(full.table$time.log.sc, 'scaled:scale') + attr(full.table$time.log.sc, 'scaled:center')), y=estimate__, ymax=upper__, ymin=lower__), fill="lightblue", alpha=0.5) +
  labs(title="", x = "Time after glacier retreat (years)", y = "Probability of EcM fungi to be present\nin the overall fungal community") +
  scale_x_continuous(expand = c(0.01, 0.01), trans = "log10", breaks=c(1, 10, 100, 500), limits = c(0.95,NA)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,1)) + #, trans = "log10", limits = c(.9,48)
  theme_bw()

## GDM plot ####
## remove variables with non-significant deviance
imp.all$sig <- imp.all$deviance
imp.all$sig[imp.all$pvalue > 0.05] <- 0 

## order variables
imp.all$variables <- ordered(imp.all$variables, levels = c("Geographic", "time.log.sc", "ndvi.s", "compo.axis1", "compo.axis2", "n.s", "ph.s", "p.s", "t.s", "twi.s"))

## Keep only first axis of plant PCoA (the second is never significant)
imp.all.plot <- imp.all %>% filter(!variables == 'compo.axis2')

## EcM
imp.all.plot %>% 
  filter(marker == "EcM") %>% 
  ggplot(aes(x=variables, y = deviance, fill = variables))+
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Variables", values = c('#DDDDDD', '#FFAABB', '#BBCC33',  '#AAAA00', '#00a878', '#EE8866', '#EEDD88', '#99DDFF', '#77AADD'),
                    labels = c("Geographic", "Time", "Productivity", "Plant community", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  scale_x_discrete(labels = c("Geographic", "Time", "Productivity", "Plant community", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  labs(x = 'Drivers of beta-diversity', y = 'Deviance explained (%)', title = 'EcM fungi') +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,25), n.breaks = 6) + #, trans = "log10", limits = c(0,60)
  theme_bw() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = -35, vjust = 0.5, hjust=0), plot.title = element_text(hjust = 0.5))

## AM
imp.all.plot %>% 
  filter(marker == "AM") %>% 
  ggplot(aes(x=variables, y = deviance, fill = variables))+
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Variables", values = c('#DDDDDD', '#FFAABB', '#BBCC33',  '#AAAA00', '#00a878', '#EE8866', '#EEDD88', '#99DDFF', '#77AADD'),
                    labels = c("Geographic", "Time", "Productivity", "Plant community", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  scale_x_discrete(labels = c("Geographic", "Time", "Productivity", "Plant community", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  labs(x = 'Drivers of beta-diversity', y = 'Deviance explained (%)', title = 'AM fungi') +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,25), n.breaks = 6) + #, trans = "log10", limits = c(0,60)
  theme_bw() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = -35, vjust = 0.5, hjust=0), plot.title = element_text(hjust = 0.5))

## Random Forest models ####
rf$variables <- ordered(rf$variables, levels = c("Glacier", "time.log", "ndvi.l", "sper.q1.l", "AM.reg", "EcM.reg", "lg_n", "ph", "lg_p", "meanT", "twi.l"))

rf %>% 
  filter(myco == "AM") %>% 
  ggplot(aes(x=variables, y = X.IncMSE, fill = variables))+
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Variables", values = c('#DDDDDD', '#FFAABB', '#BBCC33',  '#AAAA00', '#c8faa7', '#00a878', '#EE8866', '#EEDD88', '#99DDFF', '#77AADD'),
                    labels = c("Glacier", "Time", "Productivity", "Local plant diversity", "Regional AM dominance", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  scale_x_discrete(labels = c("Glacier", "Time", "Productivity", "Local plant diversity", "Regional AM dominance", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  labs(x = 'Drivers of alpha-diversity', y = 'Variable importance (% IncMSE)', title = 'AM fungi') +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,50), n.breaks = 6) +
  theme_bw() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0), plot.title = element_text(hjust = 0.5))

rf %>% 
  filter(myco == "EcM") %>% 
  ggplot(aes(x=variables, y = X.IncMSE, fill = variables))+
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Variables", values = c('#DDDDDD', '#FFAABB', '#BBCC33', '#AAAA00', '#c8faa7', '#00a878', '#EE8866', '#EEDD88', '#99DDFF', '#77AADD'),
                    labels = c("Glacier", "Time", "Productivity", "Local plant diversity", "Regional EcM dominance", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  scale_x_discrete(labels = c("Glacier", "Time", "Productivity", "Local plant diversity", "Regional EcM dominance", "Nitrogen", "pH", "Phosphorus", "Temperature", "Moisture")) +
  labs(x = 'Drivers of alpha-diversity', y = 'Variable importance (% IncMSE)', title = 'EcM fungi') +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,50), n.breaks = 6) + #, trans = "log10", limits = c(0,60)
  theme_bw() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0), plot.title = element_text(hjust = 0.5))

## Maps ####
library(rnaturalearth)
#library(ggspatial)
library(ggrepel)

## All forelands
sites_map <- full.table %>%
  dplyr::filter(!is.na(lon)) %>%
  dplyr::filter(!glacier %in% c("Pre de Bar", "Miage", "Lewis", "Kazbegi")) %>%
  group_by(glacier) %>% 
  summarise_all(mean)
length(unique(sites_map$glacier))#46

## Alps inlet
sites_alps_map <- full.table %>%
  select(glacier, lon, lat) %>%
  dplyr::filter(!glacier %in% c("Pre de Bar", "Miage", "Lewis", "Kazbegi")) %>% 
  dplyr::filter(lon > 0 & lon < 20 & lat > 44 & lat < 50) %>% 
  dplyr::filter(!is.na(lon)) %>% 
  group_by(glacier) %>% 
  summarise_all(mean)

sites_map_rest <- sites_map %>%
  dplyr::filter(!glacier %in% sites_alps_map$glacier)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_rest <- ggplot(data = world) +
  geom_sf(color = "darkgrey", fill = "white") +
  geom_point(data = sites_map_rest, mapping = aes(x = lon, y = lat), shape = 21, size = 2, alpha = .7, color = "black", fill = "red") +
  geom_point(data = sites_alps_map, mapping = aes(x = lon, y = lat), shape = 21, size = 1, alpha = .3, color = "black", fill = "red") +
  geom_rect(data=sites_map_rest, mapping=aes(xmin=5, xmax=15, ymin=44, ymax=48), color="#474747", alpha=0, size = .3) +
  geom_segment(aes(x = 5, y = 44, xend = -19, yend = 0), color="#474747", size = .2) +
  geom_segment(aes(x = 15, y = 44, xend = 120, yend = 0), color="#474747", size = .2) +
  geom_label_repel(data = sites_map_rest, mapping = aes(x = lon, y = lat, label = glacier), seed = 3, 
                   size = 3, label.padding = 0.20, force_pull = 0, force = 10,  
                   segment.size = 0.3, max.overlaps = Inf, max.iter = 10000000, max.time = 30) +
  coord_sf(ylim = c(-85, 87), expand = FALSE) + # Zoom in or out
  labs(x = 'Longitude (°)', y = 'Latitude (°)') + # Label x and y axis
  theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.background = element_rect(fill = 'aliceblue'), text = element_text(size = 14))

## Alps with glacier names simple for inlet
alps_names <- ggplot(data = world) +
  geom_sf(color = "darkgrey", fill = "white") +
  geom_point(data = sites_alps_map, mapping = aes(x = lon, y = lat), shape = 21, size = 2, alpha = .7, color = "black", fill = "red") +
  geom_label_repel(data = sites_alps_map, mapping = aes(x = lon, y = lat, label = glacier), 
                   label.padding = 0.08,  size = 2.8, seed = 10, force = 30, force_pull = 0, segment.size = 0.3,
                   max.overlaps = Inf, max.iter = 10000000, max.time = 30) +
  coord_sf(xlim = c(0, 20), ylim = c(42, 50), expand = FALSE, clip = "on") + # Zoom in or out
  theme(plot.margin = unit(c(0, 0, 0, 0), "null"), 
        axis.title=element_blank(), axis.line = element_blank(), axis.ticks=element_blank(),  axis.text=element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.background = element_rect(fill = 'aliceblue'), text = element_text(size = 6))

## World + Alps inlet
map_rest + annotation_custom(ggplotGrob(alps_names), xmin = -20, xmax = 120, 
                             ymin = -100, ymax = 20)
