rm(list=ls())
library(dplyr)

load("new_Table.Rdata")
load("compiled_data.Rdata")

################## abstract
groupByCountry = group_by(new.data, new.data$Population)
summarise(groupByCountry, sum(Est_K>1)/length(Est_K))[,2] %>% range %>% cat("range",., "\n")

print(cor.test(new_table$eff_k[new_table$nsam>15], new_table$pfpr_nonzero[new_table$nsam>15]))
print(cor.test(new_table$relatedness[new_table$nsam>15], new_table$pfpr_nonzero[new_table$nsam>15]))


print(sum(new.data$Est_K>2)/sum(new.data$Est_K>1))

sum(new.data$cluster=="Sib")/sum(new.data$Est_K>1)

sum(new.data$cluster%in%c("Sib", "HighSib"))/sum(new.data$Est_K>1)
# 54%???
sum(new.data$Est_K==2) / sum(new.data$Est_K>1)

print("Unsure about 54%")
!!!!!!!!!!!!!!!!!!!!!!!!!!

print("Unsure about 71%")
!!!!!!!!!!!!!!!!!!!!!!!!!!
# 71%???


################### Validation
load("all.error.IBD.Rdata")
load("all.error.nonIBD.Rdata")

group_by(all.error, case_panel) %>% summarise(., mean(genotypError)) %>% print
group_by(all.nonIBD.error, case_panel) %>% summarise(., mean(genotypError)) %>% print

#################### Figure 2
print( 1/(.25^2+.75^2) )

print( 1/(.55^2+.45^2) )

#################### Results
print(summary(new.data$Depth)%>% ceiling)

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

print("weighted linear model?????????")
!!!!!!!!!!!!!!!!!!!!!!!!!!
cor.test(a$rho, b$rho)

print("Overall, XXXX% of all mixed infections have 30% IBD")
a = new.data$relatedness[new.data$Est_K>1] > .3
print(sum(a, na.rm=T)/length(a))


################################## Characteristics of mixed infection correlate with local prevalence
print("prevalence range")
range(new_table$pfpr_nonzero, na.rm=T)

filter(new_table, grepl("Asia", group))$pfpr_nonzero %>% mean(., na.rm=T) %>% cat("Asia prevalence", ., "\n")

filter(new_table, grepl("Africa", group))$pfpr_nonzero %>% mean(., na.rm=T) %>% cat("Africa prevalence", ., "\n")

print("page 10 correlations")
!!!!!!!!!!!!!!!
good_table = filter(new_table, nsam > 15)
as = filter(good_table, grepl("Asia", group))
af = filter(new_table, grepl("Africa", group))
cor.test(good_table$pfpr_nonzero, good_table$eff_k)
cor.test(as$pfpr_nonzero, as$eff_k)
cor.test(af$pfpr_nonzero, af$eff_k)


cor.test(good_table$pfpr_nonzero, good_table$relatedness)
cor.test(as$pfpr_nonzero, as$relatedness)
cor.test(af$pfpr_nonzero, af$relatedness)


############################### Discussion
summarise(groupByCountry, sum(Est_K>1)/length(Est_K))[,2] %>% range %>% cat("range",., "\n")
print(sum(new.data$Est_K>2)/sum(new.data$Est_K>1))

vietnam = filter(new.data, Population=="Vietnam")

 sum(vietnam$Est_K>1) / 97
sum(vietnam$cluster=="Sib") /sum(vietnam$Est_K>1)
print("include high sib???")
sum(vietnam$cluster %in% c("Sib", "HighSib")) /sum(vietnam$Est_K>1)


