rm(list=ls()); dEploidRootDir="/home/joezhu/DEploid"
#library(DEploid)
source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))
source(paste(dEploidRootDir, "/utilities/DEploidR.R", sep=""))

args = c("-ref", "PD0577-C_ref.txt",
         "-alt", "PD0577-C_alt.txt",
         "-exclude", "asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt",
         "-dEprefix", "PD0577-C_newFilter_seed3_IBD",
         "-o", "myring")

myInput = fun.parse ( args )




myCoverageInfo = fun.extract.coverage ( myInput )

myExcludeInfo = fun.extract.exclude (myInput$excludeFileName, myInput$excludeBool)

#mypdf = function(pdfname, mypattern = "MYTEMPPNG", ...) {
#    fname = paste0(mypattern, "%05d.png")
#    gpat = paste0(mypattern, ".*\\.png")
#    takeout = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
#    if (length(takeout) > 0)
#        file.remove(takeout)
#    pngname = file.path(tempdir(), fname)
#    png(pngname, ...)
#    return(list(pdfname = pdfname, mypattern = mypattern))
#}
## copts are options to sent to convert
#mydev.off = function(pdfname, mypattern, copts = "") {
#    dev.off()
#    gpat = paste0(mypattern, ".*\\.png")
#    pngs = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
#    mystr = paste(pngs, collapse = " ", sep = "")
#    system(sprintf("convert %s -quality 100 %s %s", mystr, pdfname, copts))
#}

#res = mypdf("altVsRef_small.pdf", res = 600, height = 7, width = 7, units = "in")
#plotAltVsRef(myCoverageInfo$refCount, myCoverageInfo$altCount, cex.lab=1, title="");
#mydev.off(pdfname = res$pdfname, mypattern = res$mypattern)

png("altVsRef_small.png", res = 600, height = 7, width = 7, units = "in")
par(mar=c(5,5,2,2))
plotAltVsRef(myCoverageInfo$refCount, myCoverageInfo$altCount, cex.lab=2, title="");
dev.off()
system("convert altVsRef_small.png -quality 100 altVsRef_small.pdf")

#coverage = extractCoverageFromTxt(refFileName="PD0577-C_ref.txt", altFileName="PD0577-C_alt.txt")

#pdf("altVsRef.pdf");
postscript("altVsRef.ps");
plotAltVsRef(myCoverageInfo$refCount, myCoverageInfo$altCount, cex.lab=1, title="");
dev.off()


obsWSAF = computeObsWSAF (myCoverageInfo$altCount , myCoverageInfo$refCount)


dEploidOutput = fun.dEploidPrefix ( "PD0577-C_newFilter_seed3_IBD" )
tmpProp = read.table(dEploidOutput$propFileName, header=F)

lastProp = as.numeric(tmpProp[dim(tmpProp)[1],])

hapInfo = read.table(dEploidOutput$hapFileName, header=T)
hapChrom = hapInfo[,1]
hap = as.matrix(hapInfo[,-c(1,2)])
expWSAF = hap %*% lastProp

probs = read.table("PD0577-C_newFilter_seed3_IBD.ibd.postProb.txt", header=T)

pdf("hap_post_prob_short.pdf", height = 7, width = 18)
par(mar=c(0,0,0,0))
layout(matrix(c(1,rep(2,10),rep(3,10), 4, rep(5,10),rep(6,10)), 2, 21, byrow=T))
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext("WSAF", side=4, line=-1.5, cex = 1.5)
name1=factor("Pf3D7_11_v3", level = levels(myCoverageInfo$CHROM))
#chrom.frac = factor(c("Pf3D7_11_v3","Pf3D7_12_v3"), level = levels(myCoverageInfo$CHROM))
tmpCoveragePos = myCoverageInfo$POS[myCoverageInfo$CHROM == name1]
plot(1:length(hapInfo$POS[hapChrom == name1]), obsWSAF[myCoverageInfo$CHROM == name1][tmpCoveragePos %in% hapInfo$POS[hapChrom == name1]], col="red", pch = 16, ylab = "", yaxt='n', xaxt="n", axes=F)
points(1:length(hapInfo$POS[hapChrom == name1]), expWSAF[hapChrom == name1], col="blue", pch = 16)


name2=factor("Pf3D7_12_v3", level = levels(myCoverageInfo$CHROM))
tmpCoveragePos = myCoverageInfo$POS[myCoverageInfo$CHROM == name2]
#tmpExcludePos = myExcludeInfo$excludeTable$POS[myExcludeInfo$excludeTable$CHROM == name1]
#excludeLogic = ( tmpCoveragePos %in% tmpExcludePos )
#plotIndex = which(!excludeLogic)


plot(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16, ylab = "", type="n", yaxt='n', xaxt="n", axes=F)
rect(min(tmpCoveragePos), 0, max(tmpCoveragePos), 1, col="grey")
points(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16)
points(hapInfo$POS[hapChrom == name2], expWSAF[hapChrom == name2], col="blue", pch = 16)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n")
axis(4, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext("IBD probability", side=4, line=-1.5, cex = 1.5)

haplotypePainter(probs[probs$CHROM==name1,c(-1,-2)], labelScaling=1)
haplotypePainter(probs[probs$CHROM==name2,c(-1,-2)], labelScaling=1)
dev.off()

prop = lastProp

pdf("haps11_and_12.pdf", height = 7, width = 18)
#svg("haps11_and_12.svg", height = 7, width = 18)
par(mar=c(0,0,0,0))
par(mfrow=c(2,1))
hap = hapInfo[hapInfo$CHROM == name1, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.lab = 2, cex.main = 1.4, cex.axis=2, xaxt="n", axes=F)
    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
#    lab.at = c(1, 400, 800, 1200, 1600, 2000, 2369)
#    axis(1, at=lab.at, labels=lab.at, las=1, lwd = 0, cex=2, cex.axis=2.4)

hap = hapInfo[hapInfo$CHROM == name2, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.lab = 2, cex.main = 1.4, cex.axis=2, xaxt="n", axes=F)
    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
dev.off()


prop = lastProp

#pdf("hap_post_prob_hap.pdf", height = 7, width = 18)
#svg("hap_post_prob_hap.svg", height = 7, width = 18)
postscript("hap_post_prob_hap.ps", height = 7, width = 18)
par(mar=c(0,0,0,0))
layout(matrix(c(1,rep(2,10),rep(3,10),1,rep(2,10),rep(3,10), 1,rep(2,10),rep(3,10),
                4, rep(5,10),rep(6,10),4, rep(5,10),rep(6,10),4, rep(5,10),rep(6,10),
                7, rep(8,10),rep(9,10),7, rep(8,10),rep(9,10),7, rep(8,10),rep(9,10),
                10, rep(11,10),rep(12,10)
                ), 10, 21, byrow=T))
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext("WSAF", side=4, line=-1.5, cex = 1.5)
name1=factor("Pf3D7_11_v3", level = levels(myCoverageInfo$CHROM))
#chrom.frac = factor(c("Pf3D7_11_v3","Pf3D7_12_v3"), level = levels(myCoverageInfo$CHROM))
tmpCoveragePos = myCoverageInfo$POS[myCoverageInfo$CHROM == name1]
plot(1:length(hapInfo$POS[hapChrom == name1]), obsWSAF[myCoverageInfo$CHROM == name1][tmpCoveragePos %in% hapInfo$POS[hapChrom == name1]], col="red", pch = 16, ylab = "", yaxt='n', xaxt="n", axes=F)
points(1:length(hapInfo$POS[hapChrom == name1]), expWSAF[hapChrom == name1], col="blue", pch = 16)


name2=factor("Pf3D7_12_v3", level = levels(myCoverageInfo$CHROM))
tmpCoveragePos = myCoverageInfo$POS[myCoverageInfo$CHROM == name2]
#tmpExcludePos = myExcludeInfo$excludeTable$POS[myExcludeInfo$excludeTable$CHROM == name1]
#excludeLogic = ( tmpCoveragePos %in% tmpExcludePos )
#plotIndex = which(!excludeLogic)


plot(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16, ylab = "", type="n", yaxt='n', xaxt="n", axes=F)
rect(min(tmpCoveragePos), 0, max(tmpCoveragePos), 1, col="grey", border=NA)
points(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16)
points(hapInfo$POS[hapChrom == name2], expWSAF[hapChrom == name2], col="blue", pch = 16)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext("IBD probability", side=4, line=-1.5, cex = 1.5)

haplotypePainter <- function (posteriorProbabilities, title = "", labelScaling,
                        numberOfInbreeding = 0){
    rainbowColorBin <- 16
    rainbowColors = rainbow(rainbowColorBin)
    if ( numberOfInbreeding > 0 ){
        panelSize <- dim(posteriorProbabilities)[2]-numberOfInbreeding
        rainbowColors <- c(rep("#46a8e1", panelSize),
                           rep("#f34747", numberOfInbreeding))
    }
    barplot(t(posteriorProbabilities), beside = F, border = NA,
        col = rainbowColors, space = 0, xlab = "",
        ylab = "", main = title, cex.axis = labelScaling / 5,
        cex.lab = labelScaling / 6, cex.main = labelScaling / 5,
        xaxt = "n", yaxt = "n", axes=F)
#    newXaxt = round(seq(1, dim(posteriorProbabilities)[1], length.out = 6))
#    axis(1, at = newXaxt, labels = as.character(newXaxt),
#        cex.axis= labelScaling / 7)
#    newYaxt = seq(0, 1, length.out = 3)
#    axis(2, at = newYaxt, labels = as.character(newYaxt),
#        cex.axis= labelScaling / 7)
}

haplotypePainter(probs[probs$CHROM==name1,c(-1,-2)], labelScaling=1)
haplotypePainter(probs[probs$CHROM==name2,c(-1,-2)], labelScaling=1)


plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext("Proportion", side=4, line=-1.5, cex = 1.5)

#pdf("haps11_and_12.pdf", height = 7, width = 18)
#par(mar=c(0,0,0,0))
#par(mfrow=c(2,1))
hap = hapInfo[hapInfo$CHROM == name1, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.lab = 1.5, cex.main = 1.4, cex.axis=2, yaxt="n", axes=F)
    axis(1)
    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
#    lab.at = c(1, 400, 800, 1200, 1600, 2000, 2369)
#    axis(1, at=lab.at, labels=lab.at, las=1, lwd = 0, cex=2, cex.axis=2.4)

hap = hapInfo[hapInfo$CHROM == name2, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.lab = 1.5, cex.main = 1.4, cex.axis=2, yaxt="n", axes=F)
    axis(1)
    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
#axis(4, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
#mtext("WSAF", side=4, line=-1.5, cex = 1.5)

hap = hapInfo[hapInfo$CHROM == name1, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
plot(c(0,0),c(1,1), xlim = c(0,nhap), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
#axis(3)#, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext(paste(name1, "SNP index"), side=3, line=-4, cex = 1.5)

hap = hapInfo[hapInfo$CHROM == name2, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
plot(c(0,0),c(1,1), xlim = c(0,nhap), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
#axis(3)#, at=c(0,.5,1), labels=c(0,.5,1), cex=1.5)
mtext(paste(name2, "SNP index"), side=3, line=-4, cex = 1.5)
dev.off()





pdf("hap_post_prob_hap_07.pdf", height = 7, width = 14)
par(mar=c(0,0,0,0))
layout(matrix(c(1,rep(2,10),
                3, rep(4,10),
                5, rep(6,10)), 3, 11, byrow=T))
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n")
axis(4, at=c(0,1), labels=c(0,1))
mtext("WSAF", side=4, line=-1.5, cex = 1.5)
name1=factor("Pf3D7_07_v3", level = levels(myCoverageInfo$CHROM))
#chrom.frac = factor(c("Pf3D7_11_v3","Pf3D7_12_v3"), level = levels(myCoverageInfo$CHROM))
tmpCoveragePos = myCoverageInfo$POS[myCoverageInfo$CHROM == name1]
plot(1:length(hapInfo$POS[hapChrom == name1]), obsWSAF[myCoverageInfo$CHROM == name1][tmpCoveragePos %in% hapInfo$POS[hapChrom == name1]], col="red", pch = 16, ylab = "", yaxt='n', xaxt="n")
points(1:length(hapInfo$POS[hapChrom == name1]), expWSAF[hapChrom == name1], col="blue", pch = 16)


#name2=factor("Pf3D7_12_v3", level = levels(myCoverageInfo$CHROM))
#tmpCoveragePos = myCoverageInfo$POS[myCoverageInfo$CHROM == name2]
##tmpExcludePos = myExcludeInfo$excludeTable$POS[myExcludeInfo$excludeTable$CHROM == name1]
##excludeLogic = ( tmpCoveragePos %in% tmpExcludePos )
##plotIndex = which(!excludeLogic)


#plot(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16, ylab = "", type="n", yaxt='n', xaxt="n")
#rect(min(tmpCoveragePos), 0, max(tmpCoveragePos), 1, col="grey")
#points(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16)
#points(hapInfo$POS[hapChrom == name2], expWSAF[hapChrom == name2], col="blue", pch = 16)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n")
axis(4, at=c(0,1), labels=c(0,1))
mtext("IBD posterior probability", side=4, line=-1.5, cex = 1.5)

haplotypePainter(probs[probs$CHROM==name1,c(-1,-2)], labelScaling=1)
#haplotypePainter(probs[probs$CHROM==name2,c(-1,-2)], labelScaling=1)


plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n")
axis(4, at=c(0,1), labels=c(0,1))
mtext("Proportion", side=4, line=-1.5, cex = 1.5)

#pdf("haps11_and_12.pdf", height = 7, width = 18)
#par(mar=c(0,0,0,0))
#par(mfrow=c(2,1))
hap = hapInfo[hapInfo$CHROM == name1, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.lab = 2, cex.main = 1.4, cex.axis=2, xaxt="n", yaxt="n")
    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
dev.off()
