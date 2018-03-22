rm(list=ls())
library(maps)
library(mapdata)
library(mapproj)
library(dplyr)

deg2rad <- function(deg) {
  return (deg * (pi/180))
}

getDistanceFromLatLonInKm <- function (lat1,lon1,lat2,lon2) {
  R = 6371
#  #// Radius of the earth in km
  dLat = deg2rad(lat2-lat1)  #// deg2rad below
  dLon = deg2rad(lon2-lon1)
  a = sin(dLat/2)^2 +
      cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2)^2
  c = 2 * atan2(sqrt(a), sqrt(1-a))
  d = R * c
  return (d)
}

country_site = read.table("site_coord.txt.csv", header=T, stringsAsFactors=F, sep = "\t")
sampling_radius = 150
segments = 40
angles <- ((-segments):segments)  * pi/segments
unit.circle <- cbind(cos(angles), sin(angles)) * sampling_radius/6371*60

for (country in unique(country_site[,1])){
    a = read.csv(paste("PfPR_Data_", country, ".csv", sep=""))
    png(paste(country, ".png", sep=""))
    if ( country == "Congo" ){
      map('world', c("Democratic Republic of the Congo"), fill=T, col = "grey95")
    } else{
      map('world', c(country), fill=T, col = "grey95")
    }
    points(a$long, a$latitude, cex = 1, col="grey80", pch = 1)

    tmpTab = country_site[country_site$country == country,]
    years = sort(unique(a$year_end))
    years = years[years > 2006]
    colors = 1:20
    year_i = 0
    for ( year in years ){
        year_i = year_i + 1
        a.groupby.year = a[a$year_end == year,]
#            idx = which(getDistanceFromLatLonInKm(tmpTab$lat[i], tmpTab$long[i], a.groupby.year$latitude, a.groupby.year$longitude)<sampling_radius)
        points(a.groupby.year$long, a.groupby.year$latitude, cex = 1, col=colors[year_i], pch = 16)
#cat(year, tmpTab$site[i], "(", length(a.groupby.year$pf_pr), ")", " mean pfpr = ", sum(a.groupby.year$pf_pr, na.rm=T) / length(a.groupby.year$pf_pr) , ", mean nonzero pfpr = ", mean(a.groupby.year$pf_pr[a.groupby.year$pf_pr>0], na.rm=T),"\n")

    }


    for ( i in 1:dim(tmpTab)[1] ){
        lines(tmpTab$long[i]+unit.circle[,1], tmpTab$lat[i]+unit.circle[,2], lty= 2)
#        points(tmpTab$long[i], tmpTab$lat[i], pch = "x", col = "black", cex = 3)
        text(tmpTab$long[i], tmpTab$lat[i], labels = tmpTab$site[i], col = "black", cex = 2)
    }
    legend("topright", legend=years, col = colors[1:length(years)], pch = 16)
    dev.off()
}


new.data = data.frame( country=character(), year=character(), n.sample=numeric(), mean_pfpr=numeric(), mean_nonzero_pfpr=numeric())

for (country in unique(country_site[,1])){
    a = read.csv(paste("PfPR_Data_", country, ".csv", sep=""))
    cat(country, "(", length(a$pf_pr), ")", " mean pfpr = ", sum(a$pf_pr, na.rm=T) / length(a$pf_pr) , ", mean nonzero pfpr = ", mean(a$pf_pr[a$pf_pr>0], na.rm=T),"\n")

#    tmpTab = country_site[country_site$country == country,]
    a$year_end = as.character(a$year_end)
    a$year_end[is.na(a$year_end)] = "NA"
    years = sort(unique(a$year_end))
#    years = years[years > 2006]
    year_i = 0
    for ( year in years ){
        year_i = year_i + 1
        a.groupby.year = a[which(a$year_end == year),]
        n.sample=length(a.groupby.year$pf_pr)
        cat(year, country, "(", n.sample, ")", " mean pfpr = ", sum(a.groupby.year$pf_pr, na.rm=T) / n.sample , ", mean nonzero pfpr = ", mean(a.groupby.year$pf_pr[a.groupby.year$pf_pr>0], na.rm=T),"\n")

        new.data = rbind(new.data, data.frame( country=country, year=year, n.sample = n.sample, mean_pfpr=sum(a.groupby.year$pf_pr, na.rm=T) / n.sample, mean_nonzero_pfpr=mean(a.groupby.year$pf_pr[a.groupby.year$pf_pr>0], na.rm=T)))
    }
}


pdf("pfpr_rate_over_time.pdf", width = 13, height = 20)
colors = c("blue", "blue4", "purple", "purple4", "red", "grey", "pink", "orange4",
           "green", "yellow", "gold", "brown", "black", "wheat")
par(mfrow=c(4,1))
plot(c(1985, 2016), c(0.006,1), type="n", log="y", xlab = "year", ylab = "pfpr", main = "mean of pfpr")

for (i in 1:14){
  country = unique(new.data$country)[i]
  tmp.data = new.data[new.data$country == country,]
  tmp.data$year = as.numeric(as.character(tmp.data$year))
  lines(tmp.data$year, tmp.data$mean_pfpr, col = colors[i], lwd=2)
}
abline(h=0.1, lty=2)
abline(h=0.01, lty=2)
#abline(h=0.001, lty=2)
legend("right", legend=unique(new.data$country), col = colors,lty = 1, lwd=2)


plot(c(1985, 2016), c(0.006,1), type="n", log="y", xlab = "year", ylab = "pfpr", main = "mean of nonzero pfpr")
for (i in 1:14){
  country = unique(new.data$country)[i]
  tmp.data = new.data[new.data$country == country,]
  tmp.data$year = as.numeric(as.character(tmp.data$year))
  lines(tmp.data$year, tmp.data$mean_nonzero_pfpr, col = colors[i], lwd=2)
}
legend("right", legend=unique(new.data$country), col = colors,lty = 1, lwd=2)
abline(h=0.1, lty=2)
abline(h=0.01, lty=2)
#abline(h=0.001, lty=2)

plot(c(1985, 2016), c(0.006,1), type="n", log="y", xlab = "year", ylab = "pfpr", main = "mean of nonzero pfpr, lowess")
for (i in 1:14){
  country = unique(new.data$country)[i]
  tmp.data = new.data[new.data$country == country,]
  tmp.data$year = as.numeric(as.character(tmp.data$year))
  good = which(!is.na(tmp.data$year) & tmp.data$n.sample>0)
  lines(stats::lowess(tmp.data$year[good], tmp.data$mean_nonzero_pfpr[good], f=.8), col = colors[i], lwd=2)
}
legend("right", legend=unique(new.data$country), col = colors,lty = 1, lwd=2)
abline(h=0.1, lty=2)
abline(h=0.01, lty=2)
#abline(h=0.001, lty=2)

plot(c(1985, 2016), c(0.006,.6), type="n", xlab = "year", ylab = "pfpr", main = "mean of nonzero pfpr, lowess")
mytable = data.frame(label = character(), year = numeric(), pfpr = numeric())
for (i in 1:14){
  country = unique(new.data$country)[i]
  tmp.data = new.data[new.data$country == country,]
  tmp.data$year = as.numeric(as.character(tmp.data$year))
  good = which(!is.na(tmp.data$year) & tmp.data$n.sample>0)
  pfpr = tmp.data$mean_nonzero_pfpr[good]
  good.year = tmp.data$year[good]
  if ( length(good) > 1){
    mymodel = loess( pfpr ~ good.year, span = 1 )
    lines(good.year, predict(mymodel, good.year), col = colors[i], lwd=2)
    year = 2007:2014
    mytable = rbind(mytable, data.frame(label = paste(country, ".", year, sep=""), year = year, pfpr = predict(mymodel, year)))
  }
#
}
legend("right", legend=unique(new.data$country), col = colors,lty = 1, lwd=2)
abline(h=0.1, lty=2)
abline(h=0.01, lty=2)

dev.off()



write.table(new.data, file="pfpr_rate_over_time.txt", quote=F, row.names=F, sep="\t")

write.table(mytable, file="pfpr_rate_inferred.txt", quote=F, row.names=F, sep="\t")
