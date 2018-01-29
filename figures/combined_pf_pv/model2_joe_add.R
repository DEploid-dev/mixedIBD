rm(list=ls())

#Looking at multiple infections

#Parameters
#I[1:3] = numbers of uninfected, singly infected and multiply infected hosts
#V[1:3] = numbers of uninfected, singly infected and multiply infected vectors
#gamma = host clearance rate
#lambda.per.mos = rate of mosquito biting per day
#pi = probability that an infected mosquito establishes an infection
#nu = probability that an infected host leads to an infection in the vector
#alpha = rate at which a multiply infected host reverts to singly infected host
#phi = probability that a multiply infected vector creates a multple infection in a host (given at least some)
#eta = rate of mosquito clearance (replacement)
#theta = probability that an infection in an already singly infected mosquito leads to potential for transmission diversity
#beta = probability that a multiply infected host creates a multiply infected vector
#dt = step size rescaling parameter
#nsamp = number of samples to take
#samp.rate = rate at which to take samples

simulate.epi<-function(I, V, gamma=0.05, lambda.per.mos=0.1, pi=0.5, nu=1, alpha=0.2,
	phi=0.5, eta=0.1, beta=0.5, theta=0.2, dt=0.1, nsamp=200, samp.rate=10, plot=FALSE) {

	#Rescale biting rate to per mosquito per host
	lambda<-lambda.per.mos/sum(I);

	#Initialise
	In<-I;
	Vn<-V;
	samp<-0;
	stats<-array(0, c(nsamp, 6));
	gen<-0;

	while (samp<nsamp) {

		gen<-gen+1;

		In[1]<-((I[2]+I[3])*gamma-I[1]*(V[2]+V[3])*lambda*pi)*dt;
		In[2]<-(I[1]*(V[2]+V[3]*(1-phi))*lambda*pi+I[3]*alpha-I[2]*(gamma+lambda*pi*(V[2]+V[3])))*dt;
		In[3]<-(I[2]*(V[2]+V[3])*lambda*pi+I[1]*V[3]*phi*lambda*pi-I[3]*(alpha+gamma))*dt;

		Vn[1]<-((V[2]+V[3])*eta-V[1]*lambda*nu*(I[2]+I[3]))*dt;
		Vn[2]<-(V[1]*(I[2]+I[3]*(1-theta))*lambda*nu-V[2]*eta-V[2]*lambda*nu*beta*(I[2]+I[3]))*dt;
		Vn[3]<-(V[1]*I[3]*lambda*theta*nu+V[2]*lambda*nu*beta*(I[2]+I[3])-V[3]*eta)*dt;

		I<-I+In;
		V<-V+Vn;

		if (gen %% samp.rate == 0) {
			samp<-samp+1;
			stats[samp,]<-c(I, V);
		}
	}

	if (plot) {
		par(mfrow=c(1,2));
		plot(stats[,1], ylim=c(0,sum(I)), col="blue", xlab="Time", ylab="Proportion", type="l", main="host");
		lines(stats[,2], col="orange");
		lines(stats[,3], col="red");
		lines(stats[,3]*sum(I)/(stats[,2]+stats[,3]), lty="dotted");
		sf<-stats[,6]*phi*stats[,1]/(stats[,6]*phi*stats[,1]+(stats[,5]+stats[,6])*stats[,2]);
		lines(sf*sum(I), lty="dotted", col="pink", lwd=2);
		legend(x=0, y=sum(I), legend=c("Uninfected", "Singly", "Multiple", "M/S", "Sib rate"),
			fill=c("blue", "orange", "red", "black", "pink"), border=NA, bty="n");


		plot(stats[,4], ylim=c(0,sum(V)), col="blue", xlab="Time", ylab="Proportion", type="l", main="vector");
		lines(stats[,5], col="orange");
		lines(stats[,6], col="red");
		lines(stats[,6]*sum(V)/(stats[,5]+stats[,6]), lty="dotted");
		legend(x=0, y=sum(V), legend=c("Uninfected", "Singly", "Multiple", "M/S"),
			fill=c("blue", "orange", "red", "black"), border=NA, bty="n");

	}
	return(c(I,V));
}

simulate.epi(I=c(100,0,0), V=c(900,100,0));

#nn<-c(10,50,100,150, 200,300, 400, 500,750,1000,2000);
nn<-seq(1, 500, length.out = 100)
phi.sim<-0.3;
ss<-array(0,c(length(nn), 6));
for (i in 1:length(nn)) ss[i,]<-simulate.epi(I=c(100,0,0), V=c(0.9*nn[i], 0.1*nn[i], 0), phi=phi.sim);

pdf("prev_vs_fractions.pdf")
par(mar = c(5,5,2,2))
par(mfrow=c(1,1));

plot((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]), ss[,3]/(ss[,2]+ss[,3]), xlab="Prevalence", ylim=c(0,0.3), xlim=c(0,.5), ylab="Fraction", type="l", col="blue", cex.lab = 1.5);
lines((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]),(ss[,6]*phi.sim*ss[,1])/(ss[,6]*phi.sim*ss[,1]+(ss[,5]+ss[,6])*ss[,2]), col="orange", type="l");
legend("top",legend=c("Mixed infection", "Sib-infection"), bty="n", border=NA, fill=c("blue", "orange"));


load("snap.Rdata")
library(dplyr)

inferred.k.all.ibd_ld[adjusted.k.all.ibd_ld == 1] = 1


pfpr = read.csv("latest_pfpr_summary.txt", header=F, stringsAsFactors=F)
meta_with_yr = read.csv("pf3k_release_5_metadata_20170728_cleaned.csv", header=T, stringsAsFactors=F)
unrelated_frac = c()
sib_frac = c()
n.sample = c()
country.list = c()
for ( i in 1:dim(pfpr)[1] ){
    country = pfpr$V1[i] %>% gsub(".2.*$","",.)
    country.list = c(country.list, country)
#    print(country)
    year = pfpr$V2[i]
    a = which(meta_with_yr$country == country & meta_with_yr$year == year)
    tmpSamples = meta_with_yr$sample_name[a]
    sampleIdx = which(samples %in% tmpSamples)
    num_unrelated = sum(ibd.prob.2strain.allIBD[sampleIdx] < .1, na.rm=T) +
                    sum(ibd.prob.3strain.nonIBD[sampleIdx] > 0.9, na.rm=T)
    unrelated_frac = c(unrelated_frac, num_unrelated/length(sampleIdx))

    num_sib = sum(ibd.prob.2strain.allIBD[sampleIdx] < .7 & ibd.prob.2strain.allIBD[sampleIdx] > .3, na.rm=T) +
          sum((ibd.prob.3strain.someIBD[sampleIdx] + ibd.prob.3strain.allIBD[sampleIdx])> 0.9, na.rm=T)
    sib_frac = c(sib_frac, num_sib/length(sampleIdx))
    n.sample = c(n.sample, length(sampleIdx))
}


idx = which(n.sample > 15)
p.pch = rep("x", length(pfpr$V4))
p.pch[country.list %in% c("Thailand", "Cambodia", "Bangladesh", "Vietnam", "Myanmar", "Laos")] = "o"
points( pfpr$V4[idx], unrelated_frac[idx], col="blue", pch = p.pch[idx])
legend( "topright", legend = c("Asia countries", "Africa countries"), pch = c("o", "x"),bty="n", border=NA,)
#text( pfpr$V4[idx], unrelated_frac[idx], labels = paste(pfpr$V1, "(",n.sample,")"), col="blue")

points( pfpr$V4[idx], sib_frac[idx], col="orange", pch = p.pch[idx])
dev.off()
names(pfpr) = c("group", "year", "pfpr.low", "pfpr", "pfpr.up")
pfpr$sib_frac = sib_frac
pfpr$unrelated_frac = unrelated_frac
pfpr$n.sample = n.sample
write.table(pfpr, "pfpr_with_sib_unrelated_frac.txt", row.names=F, quote=F, sep="\t")
#End of script



pdf("prev_vs_fractions_blowup.pdf")
par(mfrow=c(1,1));
#par(mar = c(5,5,2,2))

plot((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]), ss[,3]/(ss[,2]+ss[,3]), xlab="", ylim=c(0,0.3), xlim=c(0,.005), ylab="", type="l", col="blue", cex.axis = 2);
lines((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]),(ss[,6]*phi.sim*ss[,1])/(ss[,6]*phi.sim*ss[,1]+(ss[,5]+ss[,6])*ss[,2]), col="orange", type="l");
#legend("top",legend=c("Mixed infection", "Sib-infection"), bty="n", border=NA, fill=c("blue", "orange"));


load("snap.Rdata")
library(dplyr)

inferred.k.all.ibd_ld[adjusted.k.all.ibd_ld == 1] = 1


pfpr = read.csv("latest_pfpr_summary.txt", header=F, stringsAsFactors=F)
meta_with_yr = read.csv("pf3k_release_5_metadata_20170728_cleaned.csv", header=T, stringsAsFactors=F)
unrelated_frac = c()
sib_frac = c()
n.sample = c()
country.list = c()
for ( i in 1:dim(pfpr)[1] ){
    country = pfpr$V1[i] %>% gsub(".2.*$","",.)
    country.list = c(country.list, country)
#    print(country)
    year = pfpr$V2[i]
    a = which(meta_with_yr$country == country & meta_with_yr$year == year)
    tmpSamples = meta_with_yr$sample_name[a]
    sampleIdx = which(samples %in% tmpSamples)
    num_unrelated = sum(ibd.prob.2strain.allIBD[sampleIdx] < .1, na.rm=T) +
                    sum(ibd.prob.3strain.nonIBD[sampleIdx] > 0.9, na.rm=T)
    unrelated_frac = c(unrelated_frac, num_unrelated/length(sampleIdx))

    num_sib = sum(ibd.prob.2strain.allIBD[sampleIdx] < .7 & ibd.prob.2strain.allIBD[sampleIdx] > .3, na.rm=T) +
          sum((ibd.prob.3strain.someIBD[sampleIdx] + ibd.prob.3strain.allIBD[sampleIdx])> 0.9, na.rm=T)
    sib_frac = c(sib_frac, num_sib/length(sampleIdx))
    n.sample = c(n.sample, length(sampleIdx))
}


idx = which(n.sample > 15)
p.pch = rep("x", length(pfpr$V4))
p.pch[country.list %in% c("Thailand", "Cambodia", "Bangladesh", "Vietnam", "Myanmar", "Laos")] = "o"
points( pfpr$V4[idx], unrelated_frac[idx], col="blue", pch = p.pch[idx], cex = 2)
#legend( "topright", legend = c("Asia countries", "Africa countries"), pch = c("o", "x"),bty="n", border=NA,)

points( pfpr$V4[idx], sib_frac[idx], col="orange", pch = p.pch[idx], cex = 2)
dev.off()



pdf("prev_vs_fractions.pdf")
par(mfrow=c(1,1));
par(mar = c(5,5,2,2))

plot((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]), ss[,3]/(ss[,2]+ss[,3]), xlab="Prevalence", ylim=c(0,0.5), xlim=c(0,.5), ylab="Fraction", type="l", col="blue", cex.lab = 1.5);
lines((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]),(ss[,6]*phi.sim*ss[,1])/(ss[,6]*phi.sim*ss[,1]+(ss[,5]+ss[,6])*ss[,2]), col="orange", type="l");
legend("top",legend=c("Mixed infection", "Sib-infection"), bty="n", border=NA, fill=c("blue", "orange"));


load("snap.Rdata")
library(dplyr)

inferred.k.all.ibd_ld[adjusted.k.all.ibd_ld == 1] = 1


pfpr = read.csv("latest_pfpr_summary.txt", header=F, stringsAsFactors=F)
meta_with_yr = read.csv("pf3k_release_5_metadata_20170728_cleaned.csv", header=T, stringsAsFactors=F)
unrelated_frac = c()
sib_frac = c()
n.sample = c()
country.list = c()
for ( i in 1:dim(pfpr)[1] ){
    country = pfpr$V1[i] %>% gsub(".2.*$","",.)
    country.list = c(country.list, country)
#    print(country)
    year = pfpr$V2[i]
    a = which(meta_with_yr$country == country & meta_with_yr$year == year)
    tmpSamples = meta_with_yr$sample_name[a]
    sampleIdx = which(samples %in% tmpSamples)
    num_unrelated = sum(ibd.prob.2strain.allIBD[sampleIdx] < .1, na.rm=T) +
                    sum(ibd.prob.3strain.nonIBD[sampleIdx] > 0.9, na.rm=T)
    unrelated_frac = c(unrelated_frac, num_unrelated/length(sampleIdx))

    num_sib = sum(ibd.prob.2strain.allIBD[sampleIdx] < .7 & ibd.prob.2strain.allIBD[sampleIdx] > .3, na.rm=T) +
          sum((ibd.prob.3strain.someIBD[sampleIdx] + ibd.prob.3strain.allIBD[sampleIdx])> 0.9, na.rm=T)
    sib_frac = c(sib_frac, num_sib/length(sampleIdx))
    n.sample = c(n.sample, length(sampleIdx))
}


idx = which(n.sample > 15)
p.pch = rep("x", length(pfpr$V4))
p.pch[country.list %in% c("Thailand", "Cambodia", "Bangladesh", "Vietnam", "Myanmar", "Laos")] = "o"
points( pfpr$V4[idx], unrelated_frac[idx], col="blue", pch = p.pch[idx])
legend( "topright", legend = c("Asia countries", "Africa countries"), pch = c("o", "x"),bty="n", border=NA,)
#text( pfpr$V4[idx], unrelated_frac[idx], labels = paste(pfpr$V1, "(",n.sample,")"), col="blue")

points( pfpr$V4[idx], sib_frac[idx], col="orange", pch = p.pch[idx])
dev.off()


#nn<-seq(1, 500, length.out = 100)
#phi.sim<-0.3;
#ss<-array(0,c(length(nn), 6));
#for (i in 1:length(nn)) ss[i,]<-simulate.epi(I=c(100,0,0), V=c(0.9*nn[i], 0.1*nn[i], 0), phi=phi.sim, alpha = .8);

#plot((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]), ss[,3]/(ss[,2]+ss[,3]), xlab="Prevalence", ylim=c(0,0.3), xlim=c(0,.005), ylab="Fraction", type="l", col="blue", cex.lab = 1.5);
#lines((ss[,2]+ss[,3])/(ss[,1]+ss[,2]+ss[,3]),(ss[,6]*phi.sim*ss[,1])/(ss[,6]*phi.sim*ss[,1]+(ss[,5]+ss[,6])*ss[,2]), col="orange", type="l");
#legend("top",legend=c("Mixed infection", "Sib-infection"), bty="n", border=NA, fill=c("blue", "orange"));

#points( pfpr$V4[idx], sib_frac[idx], col="orange", pch = p.pch[idx])
