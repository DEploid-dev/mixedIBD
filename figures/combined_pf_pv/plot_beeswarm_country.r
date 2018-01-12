rm(list=ls())
library(dplyr)
library(beeswarm)

load("pv-snap.Rdata")

pv.trim.idx = which(targeted.countries %in% c("Thailand", "Cambodia", "Indonesia"))
pv.adjusted.k.all.ibd_ld  =  adjusted.k.all.ibd_ld[pv.trim.idx]
targeted.countries = paste("pv.", targeted.countries[pv.trim.idx], sep="")
country.order = c("pv.Thailand", "pv.Cambodia", "pv.Indonesia", "Thailand", "Cambodia", "Vietnam", "Laos", "Myanmar", "Bangladesh", "Senegal", "The Gambia", "Malawi", "DR of the Congo", "Mali", "Nigeria", "Guinea", "Ghana")

pv.order.countries = factor(targeted.countries, levels = country.order)

ibd.prob.2 = ibd.prob.2strain.allIBD
ibd.prob.3 = 0.5 + 0.5* ibd.prob.3strain.allIBD - 0.5*ibd.prob.3strain.nonIBD
ibd.prob.4 = 0.5 + 0.5* ibd.prob.4strain.allIBD - 0.5*ibd.prob.4strain.nonIBD
ibd.prob.all.case = ifelse(is.na(ibd.prob.2), ibd.prob.3, ibd.prob.2)
ibd.prob.all.case = ifelse(is.na(ibd.prob.all.case), ibd.prob.4, ibd.prob.all.case)

pv.ibd.prob.all.case = ibd.prob.all.case[pv.trim.idx]

load("snap.Rdata")

pf.adjusted.k.all.ibd_ld  =  adjusted.k.all.ibd_ld
pf.order.countries = factor(targeted.countries, levels = country.order)

adjusted.k.all.ibd_ld = c(pv.adjusted.k.all.ibd_ld,
                          pf.adjusted.k.all.ibd_ld)

order.countries = factor(c(as.character(pv.order.countries),
                           as.character(pf.order.countries)), levels = country.order)

continent_color = c(rep("purple", 3), rep("coral", 6), rep("aquamarine", 8))

png("beeswarm_effect_k.png", width = 2000, height = 650)
par(mar = c(5,5,3,3))
par(mfrow = c(2,1))
beeswarm(adjusted.k.all.ibd_ld ~ order.countries, col = continent_color, cex.axis = 1.2, xlab= "", ylab = "Effective k", cex.lab = 1.5)

boxplot(adjusted.k.all.ibd_ld ~ order.countries, col = continent_color, cex.axis = 1.2, ylab = "Effective k", cex.lab = 1.5)
dev.off()





ibd.prob.2 = ibd.prob.2strain.allIBD
ibd.prob.3 = 0.5 + 0.5* ibd.prob.3strain.allIBD - 0.5*ibd.prob.3strain.nonIBD
ibd.prob.4 = 0.5 + 0.5* ibd.prob.4strain.allIBD - 0.5*ibd.prob.4strain.nonIBD
ibd.prob.all.case = ifelse(is.na(ibd.prob.2), ibd.prob.3, ibd.prob.2)
ibd.prob.all.case = ifelse(is.na(ibd.prob.all.case), ibd.prob.4, ibd.prob.all.case)

pf.ibd.prob.all.case = ibd.prob.all.case

ibd.prob.all.case = c(pv.ibd.prob.all.case, pf.ibd.prob.all.case)

png("beeswarm_ibd.png", width = 1500, height = 1000)

par(mfrow = c(2,1))
beeswarm(ibd.prob.all.case ~ order.countries)

boxplot(ibd.prob.all.case ~ order.countries)

dev.off()
