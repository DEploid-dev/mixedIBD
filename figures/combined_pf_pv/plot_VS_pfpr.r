rm(list=ls())
library(dplyr)
library(ggplot2)
load("snap.Rdata")
pfpr = read.csv("latest_pfpr_summary.txt", header=F, stringsAsFactors=F)
meta_with_yr = read.table("pf3k_release_5_metadata_20170804.txt", header=T, stringsAsFactors=F, sep = "\t")
mean_eff_k = c()

ibd.prob.2 = ibd.prob.2strain.allIBD
ibd.prob.3 = 0.5 + 0.5* ibd.prob.3strain.allIBD - 0.5*ibd.prob.3strain.nonIBD
ibd.prob.4 = 0.5 + 0.5* ibd.prob.4strain.allIBD - 0.5*ibd.prob.4strain.nonIBD
ibd.prob.all.case = ifelse(is.na(ibd.prob.2), ibd.prob.3, ibd.prob.2)
ibd.prob.all.case = ifelse(is.na(ibd.prob.all.case), ibd.prob.4, ibd.prob.all.case)

mean_ibd = c()
mean_ibd.2 = c()
n.sample = c()
sib_frac = c()
unrelated_frac = c()

for ( i in 1:dim(pfpr)[1] ){

    country = pfpr$V1[i] %>% gsub(".2.*$","",.)
#    print(country)
    year = pfpr$V2[i]
    tmpSamples = meta_with_yr$sample[meta_with_yr$country == country & meta_with_yr$collection_year == year]
    sampleIdx = which(samples %in% tmpSamples)
    mean_eff_k = c(mean_eff_k, adjusted.k.all.ibd_ld[sampleIdx] %>% mean(.))

    mean_ibd = c(mean_ibd, ibd.prob.all.case[sampleIdx] %>% mean(., na.rm = T))
    mean_ibd.2 = c(mean_ibd.2, ibd.prob.2[sampleIdx] %>% mean(., na.rm = T))
    n.sample = c(n.sample, length(sampleIdx))
    num_unrelated = sum(ibd.prob.2strain.allIBD[sampleIdx] < .1, na.rm=T) +
                    sum(ibd.prob.3strain.nonIBD[sampleIdx] > 0.9, na.rm=T)
    unrelated_frac = c(unrelated_frac, num_unrelated/length(sampleIdx))

    num_sib = sum(ibd.prob.2strain.allIBD[sampleIdx] < .7 & ibd.prob.2strain.allIBD[sampleIdx] > .3, na.rm=T) +
          sum((ibd.prob.3strain.someIBD[sampleIdx] + ibd.prob.3strain.allIBD[sampleIdx])> 0.9, na.rm=T)
    sib_frac = c(sib_frac, num_sib/length(sampleIdx))

}
country.list = pfpr$V1 %>% gsub(".2.*$","",.)
asia.idx.bool = country.list %in% c("Thailand", "Cambodia", "Bangladesh", "Vietnam", "Myanmar", "Laos")
africa.idx.bool = (!country.list %in% c("Thailand", "Cambodia", "Bangladesh", "Vietnam", "Myanmar", "Laos"))

pdf("prevelance.pdf", width = 16, height = 8)
par(mfrow = c(1,2))
par(mar =c(5,5,5,2))
idx = which(!is.na(mean_eff_k) & n.sample > 15)
idx.bool = (!is.na(mean_eff_k) & n.sample > 15)
#idx = which(!is.na(mean_eff_k) )
plot(mean_eff_k, pfpr$V4, type = "n", main = paste("Prevalence vs effective number of strains, correlation: ", round(cor(mean_eff_k[idx], pfpr$V4[idx]), digits = 2)),  xlim=c(1,1.8), ylim = c(0.0001,0.5), ylab="Pf Parasite Rate", xlab = "Effective number of strains", cex.lab = 1.5)
x = mean_eff_k[idx]
lm1 = lm(pfpr$V4[idx]~x)
newx <- seq(1, 1.8, length.out=100)
preds <- predict(lm1, newdata = data.frame(x=newx),
                 interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = adjustcolor('grey80', alpha.f = 0.4), border = NA)

#text(mean_eff_k, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""), col = "aquamarine")
#text(mean_eff_k, pfpr$V4, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
abline(lm1)

p.color = rep("grey80", length(mean_eff_k))
#p.color[idx] = "black"
p.color[idx.bool & asia.idx.bool] = "blue"
p.color[idx.bool & africa.idx.bool] = "orange"
country.factor = as.factor(pfpr$V1 %>% gsub(".2.*$","",.))
p.pch = country.factor %>% as.numeric
points(mean_eff_k, pfpr$V4, col = p.color, pch = p.pch, cex = 1.5)

#for (i in 1:length(mean_eff_k)){
#    lines(rep(mean_eff_k[i],2), c(pfpr$V3[i], pfpr$V5[i]), col = p.color[i])
#}
legend("bottomright", legend = levels(country.factor), pch = 1:14)
legend("topleft",legend=c("Asia", "Africa"), bty="n", border=NA, fill=c("blue", "orange"), cex = 1.5);


idx = which(!is.na(mean_ibd) & n.sample > 15)
idx.bool = (!is.na(mean_ibd) & n.sample > 15)
plot(mean_ibd, pfpr$V4, type = "n", main = paste("Prevalence vs IBD fraction, correlation: ", round(cor(mean_ibd[idx], pfpr$V4[idx]), digits = 2 )),  xlim=c(0.2,0.7), ylim = c(0,0.5), ylab="Pf Parasite Rate", xlab = "IBD fraction", cex.lab = 1.5)
#text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
#abline(lm(pfpr$V2~mean_ibd))
x = mean_ibd[idx]
lm2 = lm(pfpr$V4[idx]~x)
newx <- seq(0.1, 0.8, length.out=100)
preds <- predict(lm2, newdata = data.frame(x=newx),
                 interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = adjustcolor('grey80', alpha.f = 0.4), border = NA)

#text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""), col = "aquamarine")
#text(mean_ibd, pfpr$V4, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
abline(lm2)


#new_pfpr = predict(lm2, newdata=data.frame(x = predict_use_ibd))
#text(predict_use_ibd, new_pfpr, labels = predict_countries, col="coral")

#pv.new_pfpr = predict(lm2, newdata=data.frame(x = pv.predict_use_ibd))
#text(pv.predict_use_ibd, pv.new_pfpr, labels = paste("pv.", pv.predict_countries, sep=""), col="purple")

p.color = rep("grey80", length(mean_ibd))
#p.color[idx] = "black"
p.color[idx.bool & asia.idx.bool] = "blue"
p.color[idx.bool & africa.idx.bool] = "orange"
points(mean_ibd, pfpr$V4, col = p.color, pch = p.pch, cex = 1.5)
legend("topright", legend = levels(country.factor), pch = 1:14)
legend("bottomleft",legend=c("Asia", "Africa"), bty="n", border=NA, fill=c("blue", "orange"), cex = 1.5);

dev.off()



tick_size = 1.5
label_size = 2
pdf("prevelance_single.pdf", width = 8, height = 8)
#par(mfrow = c(1,2))
par(mar =c(5,5,5,2))
idx = which(!is.na(mean_eff_k) & n.sample > 15)
idx.bool = (!is.na(mean_eff_k) & n.sample > 15)
#idx = which(!is.na(mean_eff_k) )
plot(mean_eff_k, pfpr$V4, type = "n", main = paste("Correlation: ", round(cor(mean_eff_k[idx], pfpr$V4[idx]), digits = 2)),  xlim=c(1,1.8), ylim = c(0.0001,0.5), ylab="Pf Parasite Rate", xlab = "Effective number of strains", cex.lab = label_size, cex.main = label_size)
x = mean_eff_k[idx]
lm1 = lm(pfpr$V4[idx]~x)
newx <- seq(1, 1.8, length.out=100)
preds <- predict(lm1, newdata = data.frame(x=newx),
                 interval = 'confidence')
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = adjustcolor('grey80', alpha.f = 0.4), border = NA)

#text(mean_eff_k, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""), col = "aquamarine")
#text(mean_eff_k, pfpr$V4, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
abline(lm1)

p.color = rep("grey80", length(mean_eff_k))
#p.color[idx] = "black"
p.color[idx.bool & asia.idx.bool] = "blue"
p.color[idx.bool & africa.idx.bool] = "orange"
country.factor = as.factor(pfpr$V1 %>% gsub(".2.*$","",.))
p.pch = country.factor %>% as.numeric

my.lwd = rep(0, length(mean_eff_k))
#p.color[idx] = "black"
my.lwd[idx.bool & asia.idx.bool] = 7
my.lwd[idx.bool & africa.idx.bool] = 2

grid = c(.375, .525, .675, 1)

col.idx = lapply(mean_ibd, function(x){which( grid > x )[1]}) %>% unlist
col.array = c("red", "orange", "yellow", "green")[col.idx]


points(mean_eff_k[idx.bool], pfpr$V4[idx.bool], col = col.array[idx.bool], pch = p.pch[idx.bool], cex = 2, lwd = my.lwd[idx.bool])

#for (i in 1:length(mean_eff_k)){
#    lines(rep(mean_eff_k[i],2), c(pfpr$V3[i], pfpr$V5[i]), col = p.color[i])
#}
legend("bottomright", legend = levels(country.factor), pch = 1:14, bty = "n")
legend("topleft",legend=c("Relatedness", "[0, 0.375)", "[0.375, 0.525)", "[0.525, 0.675)", "[0.675, 1)"), bty="n", border=NA, fill=c(NA, "red", "orange", "yellow", "green"), cex = 1.5);
legend("left", legend = c("Asia", "Africa"), lty = 1, lwd = c(7,2), cex = tick_size, bty="n")


new.data = data.frame(population = pfpr$V1)
new.data$year = pfpr$V2
new.data$pfpr = pfpr$V4
new.data$idx.bool = idx.bool
new.data$mean_eff_k = mean_eff_k
new.data$mean_relatedness = mean_ibd
new.data$n.sample = n.sample
new.data$unrelated_frac = unrelated_frac
new.data$sib_frac = sib_frac
save(new.data, file = "data_for_pfpr_MAP.Rdata")

idx = which(!is.na(mean_ibd) & n.sample > 15)
idx.bool = (!is.na(mean_ibd) & n.sample > 15)



#plot(mean_ibd, pfpr$V4, type = "n", main = paste("Prevalence vs IBD fraction, correlation: ", round(cor(mean_ibd[idx], pfpr$V4[idx]), digits = 2 )),  xlim=c(0.2,0.7), ylim = c(0,0.5), ylab="Pf Parasite Rate", xlab = "IBD fraction", cex.lab = 1.5)
##text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
##abline(lm(pfpr$V2~mean_ibd))
#x = mean_ibd[idx]
#lm2 = lm(pfpr$V4[idx]~x)
#newx <- seq(0.1, 0.8, length.out=100)
#preds <- predict(lm2, newdata = data.frame(x=newx),
#                 interval = 'confidence')
#polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = adjustcolor('grey80', alpha.f = 0.4), border = NA)

##text(mean_ibd, pfpr$V2, labels = paste(pfpr$V1, "(",n.sample,")", sep=""), col = "aquamarine")
##text(mean_ibd, pfpr$V4, labels = paste(pfpr$V1, "(",n.sample,")", sep=""))
#abline(lm2)


##new_pfpr = predict(lm2, newdata=data.frame(x = predict_use_ibd))
##text(predict_use_ibd, new_pfpr, labels = predict_countries, col="coral")

##pv.new_pfpr = predict(lm2, newdata=data.frame(x = pv.predict_use_ibd))
##text(pv.predict_use_ibd, pv.new_pfpr, labels = paste("pv.", pv.predict_countries, sep=""), col="purple")

#p.color = rep("grey80", length(mean_ibd))
##p.color[idx] = "black"
#p.color[idx.bool & asia.idx.bool] = "blue"
#p.color[idx.bool & africa.idx.bool] = "orange"
#points(mean_ibd, pfpr$V4, col = p.color, pch = p.pch, cex = 1.5)
#legend("topright", legend = levels(country.factor), pch = 1:14)
#legend("bottomleft",legend=c("Asia", "Africa"), bty="n", border=NA, fill=c("blue", "orange"), cex = 1.5);

dev.off()

