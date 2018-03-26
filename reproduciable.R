rm(list=ls())
library(dplyr)

load("new_Table.Rdata")
load("compiled_data.Rdata")

cat("################################# abstract ##########################\n")
cat("#####################################################################\n")
cat("#####################################################################\n")
################## abstract
groupByCountry = group_by(new.data, new.data$Population)
summarise(groupByCountry, sum(Est_K>1)/length(Est_K))[,2] %>% range %>% cat("range",., "\n")

good_table = filter(new_table, nsam > 15)
print("effective k vs prevalence")
cor.test(good_table$eff_k, good_table$pfpr_nonzero) %>% print(.)
print("relatedness vs prevalence")
cor.test(good_table$relatedness, good_table$pfpr_nonzero) %>% print(.)


print(sum(new.data$Est_K>2)/sum(new.data$Est_K>1))

sum(new.data$cluster=="Sib")/sum(new.data$Est_K>1)

sum(new.data$cluster%in%c("Sib", "HighSib"))/sum(new.data$Est_K>1)
# 54%???
sum(new.data$Est_K==2) / sum(new.data$Est_K>1) %>%

cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Unsure about 54%\n")
#!!!!!!!!!!!!!!!!!!!!!!!!!!

cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Unsure about 71%\n")
#!!!!!!!!!!!!!!!!!!!!!!!!!!
# 71%???

cat("################################# Validation ##########################\n")
cat("#####################################################################\n")
cat("#####################################################################\n")
################### Validation
load("all.error.IBD.Rdata")
load("all.error.nonIBD.Rdata")

group_by(all.error, case_panel) %>% summarise(., mean(genotypError)) %>% print
group_by(all.nonIBD.error, case_panel) %>% summarise(., mean(genotypError)) %>% print


cat("################################# Figure2 ##########################\n")
cat("#####################################################################\n")
cat("#####################################################################\n")
#################### Figure 2
print( 1/(.25^2+.75^2) )

print( 1/(.55^2+.45^2) )

cat("################################# Results ##########################\n")
cat("#####################################################################\n")
cat("#####################################################################\n")
#################### Results
print(summary(new.data$Depth)%>% round)

print("% of mixture")
summarise(groupByCountry, percent_mix = sum(Est_K>1)/length(Est_K)) %>% print
print("% of mixture > 2 ")
summarise(groupByCountry, percent_mix = sum(Est_K>2)/length(Est_K)) %>% print

print("Dual infection Vs mixed infection")
b = summarise(groupByCountry, percent_mix = sum(Est_K>2)/length(Est_K))
a = summarise(groupByCountry, percent_mix = sum(Est_K==2)/length(Est_K))
plot(a$percent_mix, b$percent_mix)
cor.test(a$percent_mix, b$percent_mix)

print("Relatedness")
print(summarise(groupByCountry, mean(relatedness, na.rm=T)))
dual_infection = filter(new.data, Est_K==2)

print("Asia dual infection relatedness")
filter(dual_infection, group %in% levels(dual_infection$group)[c(1:3)])$relatedness %>% mean(., na.rm=T) %>% print
print("Africa dual infection relatedness")
filter(dual_infection, group %in% levels(dual_infection$group)[c(4:7)])$relatedness %>% mean(., na.rm=T) %>% print

complex_infection = filter(new.data, Est_K>2)

print("Asia complex infection relatedness")
filter(complex_infection, group %in% levels(complex_infection$group)[c(1:3)])$relatedness %>% mean(., na.rm=T) %>% print
print("Africa complex infection relatedness")
filter(complex_infection, group %in% levels(complex_infection$group)[c(4:7)])$relatedness %>% mean(., na.rm=T) %>% print


a = dual_infection%>% group_by(., Population) %>% summarise(., rho = mean(relatedness, na.rm=T))
b = complex_infection%>% group_by(., Population) %>% summarise(., rho = mean(relatedness, na.rm=T))

cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!weighted linear model?????????\n")
cor.test(a$rho, b$rho)
#!!!!!!!!!!!!!!!!!!!!!!!!!!

cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Overall, XXXX% of all mixed infections have 30% IBD\n")
a = new.data$relatedness[new.data$Est_K>1] > .3
print(sum(a, na.rm=T)/length(a))

cat("################################# Results Characteristics  ##########################\n")
cat("#####################################################################\n")
cat("#####################################################################\n")
################################## Characteristics of mixed infection correlate with local prevalence
print("prevalence range")
range(new_table$pfpr_nonzero, na.rm=T)

filter(new_table, grepl("Asia", group))$pfpr_nonzero %>% mean(., na.rm=T) %>% cat("Asia prevalence", ., "\n")

filter(new_table, grepl("Africa", group))$pfpr_nonzero %>% mean(., na.rm=T) %>% cat("Africa prevalence", ., "\n")

cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!page 10 correlations \n")
as = filter(good_table, grepl("Asia", group))
af = filter(new_table, grepl("Africa", group))
cor.test(good_table$pfpr_nonzero, good_table$eff_k)
cor.test(as$pfpr_nonzero, as$eff_k)
cor.test(af$pfpr_nonzero, af$eff_k)


cor.test(good_table$pfpr_nonzero, good_table$relatedness)
cor.test(as$pfpr_nonzero, as$relatedness)
cor.test(af$pfpr_nonzero, af$relatedness)


cat("################################# Discussion  ##########################\n")
cat("#####################################################################\n")
cat("#####################################################################\n")
############################### Discussion
summarise(groupByCountry, sum(Est_K>1)/length(Est_K))[,2] %>% range %>% cat("range",., "\n")
print(sum(new.data$Est_K>2)/sum(new.data$Est_K>1))

vietnam = filter(new.data, Population=="Vietnam")

 sum(vietnam$Est_K>1) / 97
sum(vietnam$cluster=="Sib") /sum(vietnam$Est_K>1)
cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!include high sib??\n")
sum(vietnam$cluster %in% c("Sib", "HighSib")) /sum(vietnam$Est_K>1)


