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
    points(a$long, a$latitude, cex = 1, col="grey85", pch = 21, lwd = .3)

    tmpTab = country_site[country_site$country == country,]
    years = sort(unique(a$year_end))
    years = years[years > 2006]
    colors = c("blue", "red", "purple", "green", "yellow", "brown")
    year_i = 0
    for ( year in years ){
        year_i = year_i + 1
        a.groupby.year = a[a$year_end == year,]
        col.pfpr = ifelse(is.na(a.groupby.year$pf_pr), 0, a.groupby.year$pf_pr)
        points(a.groupby.year$long, a.groupby.year$latitude, cex = 1, col = "grey85", lwd = .3, bg = lapply(col.pfpr, function(x){adjustcolor(colors[year_i], alpha.f = x)}) %>% unlist, pch = 21)
    }


    for ( i in 1:dim(tmpTab)[1] ){
        lines(tmpTab$long[i]+unit.circle[,1], tmpTab$lat[i]+unit.circle[,2], lty= 2)
        text(tmpTab$long[i], tmpTab$lat[i], labels = tmpTab$site[i], col = "black", cex = 2)
    }
    legend("topleft", legend=years, col = colors[1:length(years)], pch = 16)
    dev.off()
}

#####png(paste("world.png", sep=""), width = 1500, height = 1500)
#####par(mfrow = c(4,4))
#####for (country in unique(country_site[,1])){
#####    a = read.csv(paste("PfPR_Data_", country, ".csv", sep=""))
#####    if ( country == "Congo" ){
#####      map('world', c("Democratic Republic of the Congo"), fill=T, col = "grey95")
#####    } else{
#####      map('world', c(country), fill=T, col = "grey95")
#####    }
#####    points(a$long, a$latitude, cex = 1, col="grey85", pch = 21, lwd = .3)

#####    tmpTab = country_site[country_site$country == country,]
#####    years = sort(unique(a$year_end))
######    years = years[years > 2006]
#####    colors = c("blue", "red", "purple", "green", "yellow", "brown")
#####    year_i = 0
#####    for ( year in years ){
#####        year_i = year_i + 1
#####        a.groupby.year = a[a$year_end == year,]
#####        col.pfpr = ifelse(is.na(a.groupby.year$pf_pr), 0, a.groupby.year$pf_pr)
#####        points(a.groupby.year$long, a.groupby.year$latitude, cex = 1, col = "grey85", lwd = .3, bg = lapply(col.pfpr, function(x){adjustcolor(colors[year_i], alpha.f = x)}) %>% unlist, pch = 21)
#####    }


######    for ( i in 1:dim(tmpTab)[1] ){
######        lines(tmpTab$long[i]+unit.circle[,1], tmpTab$lat[i]+unit.circle[,2], lty= 2)
######        text(tmpTab$long[i], tmpTab$lat[i], labels = tmpTab$site[i], col = "black", cex = 2)
######    }
######    legend("topleft", legend=years, col = colors[1:length(years)], pch = 16)
#####}
#####dev.off()


new.data = data.frame( country=character(), year=character(), n.sample=numeric(),
                        mean_pfpr=numeric(), mean_nonzero_pfpr=numeric())

#                        mn.post.pi = numeric(), mn.post.c = numeric(), mn.post.wt.pi = numeric())

for (country in unique(country_site[,1])){
    a = read.csv(paste("PfPR_Data_", country, ".csv", sep=""))
    cat(country, "(", length(a$pf_pr), ")", " mean pfpr = ", sum(a$pf_pr, na.rm=T) / length(a$pf_pr) , ", mean nonzero pfpr = ", mean(a$pf_pr[a$pf_pr>0], na.rm=T),"\n")

#    tmpTab = country_site[country_site$country == country,]
#    a$year_end = as.character(a$year_end)
#    a$year_end[is.na(a$year_end)] = "NA"
    years = sort(unique(a$year_end))
#    years = years[years > 2006]
    year_i = 0
    for ( year in years ){
        year_i = year_i + 1
        a.groupby.year = a[which(a$year_end == year),]
        n.sample=length(a.groupby.year$pf_pr)
#        tmp.fit = fit.beta(a.groupby.year$examined, a.groupby.year$pf_pos, prop.c=.1, alpha = 1)
        cat(year, country, "(", n.sample, ")", " mean pfpr = ", sum(a.groupby.year$pf_pr, na.rm=T) / n.sample , ", mean nonzero pfpr = ", mean(a.groupby.year$pf_pr[a.groupby.year$pf_pr>0], na.rm=T))
#        cat(tmp.fit$mn.post.pi, tmp.fit$mn.post.c, tmp.fit$mn.post.wt.pi)
        cat("\n")
        new.data = rbind(new.data,
          data.frame( country=country, year=year, n.sample = n.sample,
            mean_pfpr=sum(a.groupby.year$pf_pr, na.rm=T) / n.sample,
            mean_nonzero_pfpr=mean(a.groupby.year$pf_pr[a.groupby.year$pf_pr>0], na.rm=T)
#            mn.post.pi = tmp.fit$mn.post.pi,
#            mn.post.c = tmp.fit$mn.post.c,
#            mn.post.wt.pi = tmp.fit$mn.post.wt.pi
            ))
    }
}

save(new.data, file="tmp.new.data.Rdata")
