rm( list= ls() )
library(maps)
library(mapdata)
library(mapproj)
library(dplyr)
library("colorRamps")


#rm(list=ls())
#library(dplyr)
#library(ggplot2)
load("snap.Rdata")
pfpr = read.csv("latest_pfpr_summary.txt", header=F, stringsAsFactors=F)

countries = unique(meta$country)[1:14]
#a = map('world', countries)

need_countries = read.table("world_country_list.csv", header=T, sep =",", quote="", stringsAsFactors=F)$Country


pfpr$country  = pfpr$V1 %>% gsub(".2.*$","",.)
average.pfpr = aggregate(pfpr$V4, by = list(pfpr$country), FUN= mean)


mean_eff_k = c()
ibd.prob.2 = ibd.prob.2strain.allIBD
ibd.prob.3 = 0.5 + 0.5* ibd.prob.3strain.allIBD - 0.5*ibd.prob.3strain.nonIBD
ibd.prob.4 = 0.5 + 0.5* ibd.prob.4strain.allIBD - 0.5*ibd.prob.4strain.nonIBD
ibd.prob.all.case = ifelse(is.na(ibd.prob.2), ibd.prob.3, ibd.prob.2)
ibd.prob.all.case = ifelse(is.na(ibd.prob.all.case), ibd.prob.4, ibd.prob.all.case)
mean_ibd = c()
#mean_ibd.2 = c()

n.sample = c()
for ( i in 1:dim(average.pfpr)[1] ){

    country = gsub("_", " ", average.pfpr$Group.1[i])
#    print(country)
#    year = pfpr$V2[i]
    tmpSamples = meta$sample[meta$country == country]
    sampleIdx = which(samples %in% tmpSamples)
    mean_eff_k = c(mean_eff_k, adjusted.k.all.ibd_ld[sampleIdx] %>% mean(.))

    mean_ibd = c(mean_ibd, ibd.prob.all.case[sampleIdx] %>% mean(., na.rm = T))
#    mean_ibd.2 = c(mean_ibd.2, ibd.prob.2[sampleIdx] %>% mean(., na.rm = T))
    n.sample = c(n.sample, length(sampleIdx))
}


coord = rbind(c(23.6850, 90.3563),
              c(12.5657, 104.9910),
              c(4.0383-8, 21.7587-2),
              c(7.9465, 1.0232-2),
              c(9.9456-1, 9.6966-2),
              c(19.8563, 102.4955),
              c(-13.2543, 34.3015),
              c(17.5707-1, 3.9962-2),
              c(21.9162, 95.9560),
              c(9.0820, 8.6753),
              c(14.4974, -14.4524),
              c(15.8700, 100.9925),
              c(13.4432, -15.3101),
              c(14.0583, 108.2772))

off_coord = rbind(c(23.6850, 90.3563),
              c(12.5657-10, 104.9910-10),
              c(4.0383-6, 21.7587+2),
              c(7.9465-8, 1.0232-12),
              c(9.9456, 9.6966),
              c(19.8563, 102.4955),
              c(-13.2543+2, 34.3015+10),
              c(17.5707, 3.9962),
              c(21.9162, 95.9560),
              c(9.0820, 8.6753),
              c(14.4974, -14.4524-5),
              c(15.8700, 100.9925),
              c(13.4432-2, -15.3101-3),
              c(14.0583, 108.2772))


pfprGrid = seq(0, 0.44, by = 0.01)
pfprColorScale =  rev(rainbow(length(pfprGrid), start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('green'))[1]))


pfpr_color = lapply(1:14, function(x){pfprColorScale[which( pfprGrid > (average.pfpr$x[x]))[1]]}) %>% unlist
#
plot_country = average.pfpr$Group.1
plot_country[3] = "Democratic Republic of the Congo"
plot_country[13] = "Gambia"

#("Bangladesh", "Cambodia"        "DR_of_the_Congo" "Ghana"
# [5] "Guinea"          "Laos"            "Malawi"          "Mali"
# [9] "Myanmar"         "Nigeria"         "Senegal"         "Thailand"
#[13] "The_Gambia"      "Vietnam"

#a =
pdf("mymap.pdf", width = 20, height = 9)
par(mar=c(1,1,1,1))
layout(matrix(c(rep(1,9),2,3,4,
                rep(1,9),5,6,7), 2,12, byrow=T))
a = map('world', need_countries, plot=F)
plot(a, type = "l", axes = F, xlab = "", ylab = "")
for ( country.idx in 1:14){
    map('world', regions = plot_country[country.idx], col = pfpr_color[country.idx], fill=T, add=TRUE)
}

effkGrid = seq(min(mean_eff_k), max(mean_eff_k)+0.01, by = 0.01)
effkColorScale =  rainbow(length(effkGrid), start=rgb2hsv(col2rgb('blue'))[1], end=rgb2hsv(col2rgb('red'))[1])
effk_color = lapply(1:14, function(x){effkColorScale[which( effkGrid > (mean_eff_k[x]))[1]]}) %>% unlist


IBDGrid = seq(min(mean_ibd), max(mean_ibd)+0.01, by = 0.01)
IBDColorScale =  rev(rainbow(length(IBDGrid), start=rgb2hsv(col2rgb('green'))[1], end=rgb2hsv(col2rgb('blue'))[1]))
IBD_color = lapply(1:14, function(x){IBDColorScale[which( IBDGrid > (mean_ibd[x]))[1]]}) %>% unlist


for ( country.idx in 1:14){
    lines(c(coord[country.idx,2], off_coord[country.idx,2]), c(coord[country.idx,1], off_coord[country.idx,1]), lwd = 2)
    points(off_coord[country.idx,2], off_coord[country.idx,1], pch = 21, cex = n.sample[country.idx]/20, bg = effk_color[country.idx], col = "black")
    points(off_coord[country.idx,2], off_coord[country.idx,1], pch = 24, cex = n.sample[country.idx]/32, bg = IBD_color[country.idx], col = "black")
}


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

#    dev.new(width=1.75, height=5)
    plot(c(0,5), c(min,max), xlim = c(0,10),type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab="", main=title, cex.lab = 2)
    axis(2, ticks, las=1, cex.axis = 1.5)
    for (i in 1:(length(lut)-1)) {
    	y = (i-1)/scale + min
    	rect(0,y,5,y+1/scale, col=lut[i], border=NA)
    }
}

plot(c(0,0),c(1,1), type="n", axes=F)
plot(c(0,0),c(1,1), type="n", axes=F)
plot(c(0,0),c(1,1), type="n", axes=F)
#lut = rev(rainbow(100, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('blue'))[1]))
color.bar(pfprColorScale, min(pfprGrid), max(pfprGrid),title = "Prevalence")
color.bar(effkColorScale, min(effkGrid), max(effkGrid), title = "Effective K")
color.bar(IBDColorScale, min(IBDGrid), max(IBDGrid), title = "Relatedness")


dev.off()
