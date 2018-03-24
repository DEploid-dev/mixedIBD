rm(list=ls())

#To get statistics from the summary data

#setwd("C:/Users/mcvean/Dropbox/Joe/MixedIBD/fig3_data");
load("compiled_data.Rdata");

pops<-unique(new.data$Pop);
cont<-new.data$group[match(pops, new.data$Pop)];
af<-pops[grep("Afr", cont)];
as<-pops[grep("Asi", cont)];
cont.label<-rep(1, length(pops));
cont.label[pops %in% af]<-2;

#Average MI rate per population
op<-array(0, c(length(pops), 10));
colnames(op)<-c("Mean.K", "Mean.KG1", "Mean.KG2", "Mean.KE2",
	"Mean.rel.dual", "Mean.rel.triplus", "N", "N2", "NG2", "Mean.IBD");
for (i in 1:length(pops)) {
	wi<-which(new.data$Pop==pops[i]);
	tmp<-new.data[wi,];
	op[i,1:4]<-c(mean(new.data$Est_K[wi]), mean(new.data$Est_K[wi]>1), mean(new.data$Est_K[wi]>2), mean(new.data$Est_K[wi]==2));
	op[i,5]<-mean(tmp$rel[tmp$Est_K==2]);
	op[i,6]<-mean(tmp$rel[tmp$Est_K>2]);
	op[i,7]<-nrow(tmp);
	op[i,8]<-sum(tmp$Est_K==2);
	op[i,9]<-sum(tmp$Est_K>2);
	op[i,10]<-mean(tmp$rel, na.rm=T);
}

#Range of MI in Africa
wp<-which(cont.label==2);
mn<-which.min(op[wp,2]);
mx<-which.max(op[wp,2]);
cat("\nRange of MI in Africa = ", op[wp[mn],2], ", (", pops[wp[mn]], "), to ", op[wp[mx],2], ", (", pops[wp[mx]], ")", sep="");


#Range of MI in Asia
wp<-which(cont.label==1);
mn<-which.min(op[wp,2]);
mx<-which.max(op[wp,2]);
cat("\nRange of MI in Asia = ", op[wp[mn],2], ", (", pops[wp[mn]], "), to ", op[wp[mx],2], ", (", pops[wp[mx]], ")", sep="");


#Look at per year effects
for (i in 1:length(pops)) {
	wi<-which(new.data$Pop==pops[i]);
	tmp<-new.data[wi,];
	yy<-unique(tmp$year);
	cat("\n\n************** Pop:", pops[i], "**************\n\n");
	if (length(yy)>1) print(summary(lm(tmp$Est_K ~ as.factor(tmp$year))));
}


#Range of >2 infection
mn<-which.min(op[,3]);
mx<-which.max(op[,3]);
cat("\n\nRange of >2 infection = ", op[mn,3], ", (", pops[mn], "), to ", op[mx,3], ", (", pops[mx], ")", sep="");

#Correlation dual and higher rates infection
lma<-print(summary(lm(op[,3]~op[,4]+cont.label)));
cat("\n\nCorrelation between dual and higher rates of infection = ", cor(op[,4], op[,3], method="pearson"), sep="");


#IBD in dual infections
mn<-which.min(op[,5]);
mx<-which.max(op[,5]);
cat("\n\nRange of IBD in dual infection = ", op[mn,5], ", (", pops[mn], "), to ", op[mx,5], ", (", pops[mx], ")", sep="");

#Asia v Africa IBD
cat("\n\nMean IBD.dual in Asia = ", mean(op[cont.label==1,5]), " and mean IBD.dual in Africa = ", mean(op[cont.label==2,5]), sep="");


#IBD in triple and more
mn<-which.min(op[,6]);
mx<-which.max(op[,6]);
cat("\n\nRange of IBD in >2 infection = ", op[mn,6], ", (", pops[mn], "), to ", op[mx,6], ", (", pops[mx], ")", sep="");


#Asia v Africa IBD
cat("\n\nMean IBD.G2 in Asia = ", mean(op[cont.label==1,6]), " and mean IBD.G2 in Africa = ", mean(op[cont.label==2,6]), sep="");


#Sib level relatedness
tmp<-sum(new.data$rel>0.3, na.rm=T)/sum(new.data$rel>0.0001, na.rm=T);
cat("\n\nOverall, ", 100*signif(tmp, 3), "% of all mixied infections have IBD>0.3", sep="");

#Correlations and summary plot
labs<-c("THA", "CAM", "VIE", "LAO", "BAN", "MYA", "MALA", "DRC", "GHA", "MALI",
	"NIG", "SEN", "GAM", "GUI");

pdf(file="mixInfSum.pdf", height=6, width=5);
par(mar = c(7, 5, 1, 1))
plot(op[,2], op[,10], type="n", col=c("blue", "red")[cont.label], pch=19,
	xlab="Fraction mixed infection",
	ylab="Mean IBD in mixed infections",
	xlim=c(min(op[,2])*0.85, max(op[,2])*1.15),
	ylim=c(min(op[,10])*0.85, max(op[,10])*1.15),
	bty="n",
	main="", cex.lab = 1.5
);

sex<-sqrt(op[,2]*(1-op[,2])/op[,7]);
sey<-sqrt(op[,10]*(1-op[,10])/(op[,8]+op[,9]));
fac<-1;
segments(op[,2]-fac*sex, op[,10], op[,2]+fac*sex, op[,10], col=grey(0.75));
segments(op[,2], op[,10]-fac*sey, op[,2], op[,10]+fac*sey, col=grey(0.75));
points(op[,2], op[,10], col=c("royalblue", "orangered1")[cont.label], pch=19)
text(x=op[,2], y=op[,10], labels=labs, cex=1.1, pos=1);
slm<-summary(lm(op[-12,10]~op[-12,2], weights=op[-12,7]));
#slm<-summary(lm(op[,10]~op[,2]+as.factor(cont.label), weights=op[,7]));
print(slm)
#abline(slm$coef[1,1], slm$coef[2,1], lty="dotted", col="royalblue");
#abline(slm$coef[1,1]+slm$coef[3,1], slm$coef[2,1], lty="dotted", col="orangered1");
abline(slm$coef[1,1], slm$coef[2,1], lty="dotted", col="black");

dev.off();




#Plot K results
tab<-array(0, c(length(pops), 4));
for (i in 1:length(pops)) {
	tmp<-new.data[new.data$Pop == pops[i],];
	tab[i,]<-c(mean(tmp$Est_K==1), mean(tmp$Est_K==2), mean(tmp$Est_K==3), mean(tmp$Est_K==4));
}
c1<-which(cont.label==1);
c2<-which(cont.label==2);
oo<-c(c1[order(tab[c1,1], decreasing=T)], c2[order(tab[c2,1], decreasing=T)]);

par(mfrow=c(1,2));
sp<-c(rep(0.1, length(c1)), 1, rep(0.1, length(c2)-1));

pdf(file="mixK.pdf", height=6, width=6);
par(mar = c(14, 1, 1, 1))
layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
v<-barplot(t(tab[oo,]),
	names=pops[oo],
	cex.names=2, las=2, yaxt="n",
	col=c("papayawhip", "orange", "red", "grey50"),
	ylim=c(0,1.1),
	border=NA,
	space=sp
);

plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,10), ylim=c(0,10), bty="n");
legend("top", legend=c("Clonal", "Dual", "Triple", "More than 3"), fill=c("papayawhip", "orange", "red", "grey50"),
	bty="n", border=NA, cex = 2);
legend("bottom", legend=c("r < 0.1", "0.1 < r < 0.3", "0.3 < r < 0.7", "r > 0.7"), fill=c("darkorange4", "darkorange3", "orange", "moccasin"),
	bty="n", border=NA, cex = 2);

dev.off();



#Plot IBD results
tab<-array(0, c(length(pops), 4));
for (i in 1:length(pops)) {
	tmp<-new.data[new.data$Pop == pops[i],];
	tab[i,]<-c(mean(tmp$rel<=0.1, na.rm=T),
		mean(tmp$rel>0.1 & tmp$rel<=0.3, na.rm=T),
		mean(tmp$rel>0.3 & tmp$rel<=0.7, na.rm=T),
		mean(tmp$rel>0.7, na.rm=T));
}

pdf(file="mixIBD.pdf", height=6, width=4.05);
par(mar = c(9, 1, 0, 0))
v<-barplot(t(tab[oo,]),
	names=pops[oo],
	cex.names=1.25, las=2, yaxt="n",
	col=c("darkorange4", "darkorange3", "orange", "moccasin"),
	ylim=c(0,1.1),
	border=NA,
	space=sp
);
dev.off();




#Legend

#pdf(file="legend.pdf", height=5, width=4);
#dev.off();














