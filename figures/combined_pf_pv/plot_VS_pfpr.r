rm(list=ls())
library(dplyr)
library(ggplot2)
load("snap.Rdata")
pfpr = read.table("pfpr_from_map.txt", header=F, stringsAsFactors=F)
meta_with_yr = read.csv("pf3k_release_5_metadata_20170728_cleaned.csv", header=T, stringsAsFactors=F)
mean_eff_k = c()

ibd.prob.2 = ibd.prob.2strain.allIBD
ibd.prob.3 = 0.5 + 0.5* ibd.prob.3strain.allIBD - 0.5*ibd.prob.3strain.nonIBD
ibd.prob.4 = 0.5 + 0.5* ibd.prob.4strain.allIBD - 0.5*ibd.prob.4strain.nonIBD
ibd.prob.all.case = ifelse(is.na(ibd.prob.2), ibd.prob.3, ibd.prob.2)
ibd.prob.all.case = ifelse(is.na(ibd.prob.all.case), ibd.prob.4, ibd.prob.all.case)

mean_ibd = c()
mean_ibd.2 = c()
n.sample = c()
for ( i in 1:dim(pfpr)[1] ){

    country = pfpr$V1[i] %>% gsub(".2.*$","",.)
#    print(country)
    year = pfpr$V3[i]
    tmpSamples = meta_with_yr$sample_name[meta_with_yr$country == country & meta_with_yr$year == year]
    sampleIdx = which(samples %in% tmpSamples)
    mean_eff_k = c(mean_eff_k, adjusted.k.all.ibd_ld[sampleIdx] %>% mean(.))

    mean_ibd = c(mean_ibd, ibd.prob.all.case[sampleIdx] %>% mean(., na.rm = T))
    mean_ibd.2 = c(mean_ibd.2, ibd.prob.2[sampleIdx] %>% mean(., na.rm = T))
    n.sample = c(n.sample, length(sampleIdx))

}

predict_countries = c("Thailand","Cambodia","Vietnam","Laos", "Bangladesh", "Myanmar")
predict_use_eff_k = c()
predict_use_ibd = c()
for (country_i in predict_countries){
    tmpSamples = meta$sample[meta$country == country_i]
    sampleIdx = which(samples %in% tmpSamples)
    predict_use_eff_k = c(predict_use_eff_k, adjusted.k.all.ibd_ld[sampleIdx] %>% mean(.))
    predict_use_ibd = c(predict_use_ibd, ibd.prob.all.case[sampleIdx] %>% mean(.,na.rm = T))

}

load("pv-snap.Rdata")
pv.predict_countries = c("Thailand","Cambodia","Indonesia")
pv.predict_use_eff_k = c()
pv.predict_use_ibd = c()
for (country_i in pv.predict_countries){
    tmpSamples = meta$Sample.ID[meta$Country == country_i]
    sampleIdx = which(samples %in% tmpSamples)
    pv.predict_use_eff_k = c(pv.predict_use_eff_k, adjusted.k.all.ibd_ld[sampleIdx] %>% mean(.))
    pv.predict_use_ibd = c(pv.predict_use_ibd, ibd.prob.all.case[sampleIdx] %>% mean(.,na.rm = T))

}


pdf("prevelance.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
par(mar =c(5,5,5,2))
idx = which(!is.na(mean_eff_k) & n.sample > 5)
#idx = which(!is.na(mean_eff_k) )
plot(mean_eff_k, pfpr$V2, type = "n", main = paste("Prevalence vs effective number of strains, correlation: ", round(cor(mean_eff_k[idx], pfpr$V2[idx]), digits = 2)),  xlim=c(1,1.8), ylim = c(0,0.5), ylab="Pf prevalence", xlab = "Effective number of strains", cex.lab = 1.5)
x = mean_eff_k[idx]
lm1 = lm(pfpr$V2[idx]~x)
newx <- seq(1, 1.8, length.out=100)
preds <- predict(lm1, newdata = data.frame(x=newx),
                 interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = adjustcolor('grey80', alpha.f = 0.4), border = NA)

#text(mean_eff_k, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""), col = "aquamarine")
text(mean_eff_k, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
abline(lm1)

new_pfpr = predict(lm1, newdata=data.frame(x = predict_use_eff_k))
text(predict_use_eff_k, new_pfpr, labels = predict_countries, col="coral")

pv.new_pfpr = predict(lm1, newdata=data.frame(x = pv.predict_use_eff_k))

text(pv.predict_use_eff_k, pv.new_pfpr, labels = paste("pv.", pv.predict_countries, sep=""), col="purple")


idx = which(!is.na(mean_ibd) & n.sample > 5)
plot(mean_ibd, pfpr$V2, type = "n", main = paste("Prevalence vs IBD fraction, correlation: ", round(cor(mean_ibd[idx], pfpr$V2[idx]), digits = 2 )),  xlim=c(0.2,0.7), ylim = c(0,0.5), ylab="Pf prevalence", xlab = "IBD fraction", cex.lab = 1.5)
#text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
#abline(lm(pfpr$V2~mean_ibd))
x = mean_ibd[idx]
lm2 = lm(pfpr$V2[idx]~x)
newx <- seq(0.1, 0.8, length.out=100)
preds <- predict(lm2, newdata = data.frame(x=newx),
                 interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = adjustcolor('grey80', alpha.f = 0.4), border = NA)

#text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""), col = "aquamarine")
text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
abline(lm2)

new_pfpr = predict(lm2, newdata=data.frame(x = predict_use_ibd))
text(predict_use_ibd, new_pfpr, labels = predict_countries, col="coral")

pv.new_pfpr = predict(lm2, newdata=data.frame(x = pv.predict_use_ibd))
text(pv.predict_use_ibd, pv.new_pfpr, labels = paste("pv.", pv.predict_countries, sep=""), col="purple")


dev.off()
