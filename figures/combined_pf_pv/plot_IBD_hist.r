rm(list=ls())
library(dplyr)
library(beeswarm)

load("pv-snap.Rdata")

pv.trim.idx = which(targeted.countries %in% c("Thailand", "Cambodia", "Indonesia"))
ibd.prob.2 = ibd.prob.2strain.allIBD[pv.trim.idx]
ibd.prob.2.non.idx = which(!is.na(ibd.prob.2))
ibd.prob.2.non.na = ibd.prob.2[ibd.prob.2.non.idx]
pv.ibd = ibd.prob.2.non.na[(.05 < ibd.prob.2.non.na) & (.95 > ibd.prob.2.non.na)]

load("snap.Rdata")

ibd.prob.2 = ibd.prob.2strain.allIBD
ibd.prob.2.non.idx = which(!is.na(ibd.prob.2))
ibd.prob.2.non.na = ibd.prob.2[ibd.prob.2.non.idx]

ibd.prob.2.samples = samples[ibd.prob.2.non.idx]

tmp.countries = meta$country[lapply(ibd.prob.2.samples, function(x){return(which(x == meta$sample))}) %>% unlist]

continent = ifelse(tmp.countries %in% c("Thailand", "Cambodia", "Vietnam", "Laos", "Bangladesh", "Myanmar"), "Asia", "Africa")

highlight.idx = (.05 < ibd.prob.2.non.na) & (.95 > ibd.prob.2.non.na)

png("ibd_frac_hist.png", width = 1000, height = 450)
layout(matrix(c(1,2,3,1,2,3,4,5,6), nrow = 3, ncol = 3, byrow = T))
#tmp.data = c(ibd.prob.2.non.na[highlight.idx])
hist(pv.ibd, main = paste("Vivax samples of 2 strains IBD fraction > 0.05, < 0.95"), probability=T, xlab="fraction", xlim = c(0,1), cex.lab = 1.5, cex.main = 1.5,col="purple")
lines(density(pv.ibd))

tmp.data = c(ibd.prob.2.non.na[highlight.idx & continent=="Africa"])
hist(tmp.data, main = paste("African samples of 2 strains IBD fraction > 0.05, < 0.95"), probability=T, xlab="fraction", xlim = c(0,1),cex.lab = 1.5, cex.main = 1.5,col="aquamarine")
lines(density(tmp.data))

tmp.data = c(ibd.prob.2.non.na[highlight.idx & continent=="Asia"])
hist(tmp.data, main = paste("Asian samples of 2 strains IBD fraction > 0.05, < 0.95"), probability=T, xlab="fraction", xlim = c(0,1),cex.lab = 1.5, cex.main = 1.5, col="coral")
lines(density(tmp.data))

#tmp.data = c(ibd.prob.2.non.na[highlight.idx])
beeswarm(pv.ibd, method="center", vertical=F, xlim = c(0,1),col="purple")

tmp.data = c(ibd.prob.2.non.na[highlight.idx & continent=="Africa"])
beeswarm(tmp.data, method="center", vertical=F, xlim = c(0,1),col="aquamarine")

tmp.data = c(ibd.prob.2.non.na[highlight.idx & continent=="Asia"])
beeswarm(tmp.data, method="center", vertical=F, xlim = c(0,1), col="coral")
dev.off()
