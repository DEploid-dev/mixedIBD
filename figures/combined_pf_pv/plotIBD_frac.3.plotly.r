rm(list=ls())
library(plotly)


#years = read.table("year-all.tab", header=T)

#library(beeswarm)
##a = read.table("ibd.fraction.txt", stringsAsFactors=F, header=T)
#ibd = read.table("ibd.fraction.3.tmp.txt", stringsAsFactors=F, header=T)


axis <- function(title) {
  list(
    title = title,
    titlefont = list(
      size = 10
    ),
    tickfont = list(
      size = 15
    ),
    tickcolor = 'rgba(0,0,0,0)',
    ticklen = 5
  )
}

#ibd.3strain = data.frame(samples, nonIBD = ibd.prob.3strain.nonIBD,
#  allIBD = ibd.prob.3strain.allIBD,
#  someIBD = ibd.prob.3strain.someIBD,
#  country = targeted.countries
#)

#p.3strain <- ibd.3strain %>%
#  plot_ly() %>%
#  add_trace(
#    type = 'scatterternary',
#    mode = 'markers',
#    a = ~nonIBD,
#    b = ~someIBD,
#    c = ~allIBD,
#    color = ~country,
#    text = ~samples,
#    marker = list(
#      symbol = 100,
##      color = '#DB7365',
#      size = 14,
#      line = list('width' = 2)
#    )
#  ) %>%
#  layout(
#    title = "Fraction of IBD",
#    ternary = list(
#      aaxis = axis('A: nonIBD'),
#      baxis = axis('B: someIBD'),
#      caxis = axis('C: allIBD')
#    )
#  )

#ibd.4strain = data.frame(samples, nonIBD = ibd.prob.4strain.nonIBD,
#  allIBD = ibd.prob.4strain.allIBD,
#  someIBD = ibd.prob.4strain.someIBD,
#  country = targeted.countries
#)


#p.4strain <- ibd.4strain %>%
#  plot_ly() %>%
#  add_trace(
#    type = 'scatterternary',
#    mode = 'markers',
#    a = ~nonIBD,
#    b = ~someIBD,
#    c = ~allIBD,
#    color = ~country,
#    text = ~samples,
#    marker = list(
#      symbol = 100,
##      color = '#DB7365',
#      size = 14,
#      line = list('width' = 2)
#    )
#  ) %>%
#  layout(
#    title = "Fraction of IBD",
#    ternary = list(
#      aaxis = axis('A: nonIBD'),
#      baxis = axis('B: someIBD'),
#      caxis = axis('C: allIBD')
#    )
#  )


## Create a shareable link to your chart
## Set up API credentials: https://plot.ly/r/getting-started
#chart_link = api_create(p, filename="ternary/basic")
#chart_link

country.order = c("Thailand", "Cambodia", "Vietnam", "Laos", "Myanmar", "Bangladesh", "Senegal", "The Gambia", "Malawi", "DR of the Congo", "Mali", "Nigeria", "Guinea", "Ghana", "pv.Thailand", "pv.Cambodia", "pv.Indonesia")

load("pv-snap.Rdata")

ibd.prob.2strain.nonIBD = 1-ibd.prob.2strain.allIBD
ibd.prob.mixed.nonIBD = ifelse(is.na(ibd.prob.2strain.nonIBD), ibd.prob.3strain.nonIBD, ibd.prob.2strain.nonIBD)
ibd.prob.mixed.nonIBD = ifelse(is.na(ibd.prob.mixed.nonIBD), ibd.prob.4strain.nonIBD, ibd.prob.mixed.nonIBD)

ibd.prob.mixed.allIBD = ifelse(is.na(ibd.prob.2strain.allIBD), ibd.prob.3strain.allIBD, ibd.prob.2strain.allIBD)
ibd.prob.mixed.allIBD = ifelse(is.na(ibd.prob.mixed.allIBD), ibd.prob.4strain.allIBD, ibd.prob.mixed.allIBD)

ibd.prob.mixed.someIBD = ifelse(!is.na(ibd.prob.2strain.allIBD), 0, ibd.prob.3strain.someIBD)
ibd.prob.mixed.someIBD = ifelse(is.na(ibd.prob.mixed.someIBD), ibd.prob.4strain.someIBD, ibd.prob.mixed.someIBD)

pv.trim.idx = which(targeted.countries %in% c("Thailand", "Cambodia", "Indonesia"))
pv.ibd.prob.mixed.nonIBD = ibd.prob.mixed.nonIBD[pv.trim.idx]
pv.ibd.prob.mixed.allIBD = ibd.prob.mixed.allIBD[pv.trim.idx]
pv.ibd.prob.mixed.someIBD = ibd.prob.mixed.someIBD[pv.trim.idx]
targeted.countries = paste("pv.", targeted.countries[pv.trim.idx], sep="")

pv.order.countries = factor(targeted.countries, levels = country.order)
pv.samples = samples[pv.trim.idx]

load("snap.Rdata")

ibd.prob.2strain.nonIBD = 1-ibd.prob.2strain.allIBD
ibd.prob.mixed.nonIBD = ifelse(is.na(ibd.prob.2strain.nonIBD), ibd.prob.3strain.nonIBD, ibd.prob.2strain.nonIBD)
ibd.prob.mixed.nonIBD = ifelse(is.na(ibd.prob.mixed.nonIBD), ibd.prob.4strain.nonIBD, ibd.prob.mixed.nonIBD)

ibd.prob.mixed.allIBD = ifelse(is.na(ibd.prob.2strain.allIBD), ibd.prob.3strain.allIBD, ibd.prob.2strain.allIBD)
ibd.prob.mixed.allIBD = ifelse(is.na(ibd.prob.mixed.allIBD), ibd.prob.4strain.allIBD, ibd.prob.mixed.allIBD)

ibd.prob.mixed.someIBD = ifelse(!is.na(ibd.prob.2strain.allIBD), 0, ibd.prob.3strain.someIBD)
ibd.prob.mixed.someIBD = ifelse(is.na(ibd.prob.mixed.someIBD), ibd.prob.4strain.someIBD, ibd.prob.mixed.someIBD)

pf.ibd.prob.mixed.nonIBD = ibd.prob.mixed.nonIBD
pf.ibd.prob.mixed.allIBD = ibd.prob.mixed.allIBD
pf.ibd.prob.mixed.someIBD = ibd.prob.mixed.someIBD


continent = ifelse(targeted.countries %in% c("Thailand", "Cambodia", "Vietnam", "Laos", "Bangladesh", "Myanmar"), "Asia", "Africa")
continent_color = ifelse(targeted.countries %in% c("Thailand", "Cambodia", "Vietnam", "Laos", "Bangladesh", "Myanmar"), "coral", "aquamarine")

continent_color = c(continent_color, rep("purple", length(pv.samples)))
#pf.targeted.countries.fac = factor(targeted.countries, levels = unique(targeted.countries))
pf.order.countries = factor(targeted.countries, levels = country.order)

targeted.countries.fac = factor(c(as.character(pf.order.countries),
                           as.character(pv.order.countries)), levels = country.order)

continent_color

ibd.mixed = data.frame(samples = c(samples, pv.samples),
  nonIBD = c(pf.ibd.prob.mixed.nonIBD, pv.ibd.prob.mixed.nonIBD),
  allIBD = c(pf.ibd.prob.mixed.allIBD, pv.ibd.prob.mixed.allIBD),
  someIBD = c(pf.ibd.prob.mixed.someIBD, pv.ibd.prob.mixed.someIBD),
  country = targeted.countries.fac
#  continent = continent
)


p.mixed <- ibd.mixed %>%
  plot_ly() %>%
  add_trace(
    type = 'scatterternary',
    mode = 'markers',
    a = ~someIBD,
    b = ~nonIBD,
    c = ~allIBD,
    color = ~country,
    text = ~samples,
    marker = list(
      symbol = 100,
      color = continent_color,
      size = 4,
      line = list('width' = 2)
    )
  ) %>%
  layout(
#    title = "Fraction of IBD",
    title = "",
    ternary = list(
      aaxis = axis('A: someIBD'),
      baxis = axis('B: nonIBD'),
      caxis = axis('C: allIBD')
    )
  )
