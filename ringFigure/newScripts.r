rm(list=ls()); dEploidRootDir="/home/joezhu/DEploid"
#library(DEploid)
library("dplyr")

source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))
source(paste(dEploidRootDir, "/utilities/DEploidR.R", sep=""))

args = c("-ref", "PD0577-C_ref.txt",
         "-alt", "PD0577-C_alt.txt",
         "-exclude", "asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt",
         "-dEprefix", "PD0577-C_newFilter_seed3_IBD",
         "-o", "myring")

myInput = fun.parse ( args )


tick_size = 1.5
label_size = 2


myCoverageInfo = fun.extract.coverage ( myInput )

myExcludeInfo = fun.extract.exclude (myInput$excludeFileName, myInput$excludeBool)

plotAltVsRef <- function ( ref, alt, title = "",
                    cex.lab = label_size, cex.main = 1, cex.axis = tick_size ){
    cr <- colorRampPalette(colors = c("#de2d26", "#2b8cbe"))
    colors <- cr(31)
    ratios <- ref / (ref + alt + 0.0000001)
    tmpRange <- 260
    plot ( ref, alt, xlim = c(0, tmpRange), ylim = c(0, tmpRange),
        pch = 20, col = scales::alpha(colors[ceiling(ratios * 30) + 1], 0.7),
        xlab = "", ylab = "", main = title,
        cex = 0.5, cex.lab = label_size, cex.main = cex.main, cex.axis = cex.axis, axes=F)
        tick_seq = seq(0,260, by  =65)
axis(1, at=tick_seq, labels=tick_seq, cex.axis=tick_size)
mtext("No. REF reads", side=1, line=2.8, cex = label_size)

axis(2, at=tick_seq, labels=tick_seq, cex.axis=tick_size)
mtext("No. ALT reads", side=2, line=2.2, cex = label_size)

    legend("topright", legend = c("100% ALT", "100% REF", "50/50"),
        fill = colors[c(1, 31, 15)], cex = label_size, border = NA, box.lwd = 0,
        box.col = "white", bg = NA)
    abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
#    points(exclude.ref, exclude.alt, col = "red")
    abline(v = 65, untf = FALSE, lty = 2)
    abline(h = 65, untf = FALSE, lty = 2)

    abline(h = 130, untf = FALSE, lty = 2)
    abline(v = 130, untf = FALSE, lty = 2)

}


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




#png("altVsRef_small.png", res = 600, height = 7, width = 7, units = "in")
#par(mar=c(5,5,2,2))
#plotAltVsRef(myCoverageInfo$refCount, myCoverageInfo$altCount, cex.lab=2, title="");
#dev.off()
#system("convert altVsRef_small.png -quality 100 altVsRef_small.pdf")

#coverage = extractCoverageFromTxt(refFileName="PD0577-C_ref.txt", altFileName="PD0577-C_alt.txt")

alt.nonzero.idx = which(myCoverageInfo$altCount > 0)
set.seed(1)
alt.zero.idx = sample(which(myCoverageInfo$altCount == 0), 5000)
plot.idx = c(alt.nonzero.idx, alt.zero.idx)

pdf("altVsRef.pdf", width = 6, height = 6);
par(mar = c(5,5,2,2))
plotAltVsRef(myCoverageInfo$refCount[plot.idx], myCoverageInfo$altCount[plot.idx]);
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


prop = lastProp

tick_size = 1.5*1.2
label_size = 2


pdf("hap_post_prob_hap.pdf", height = 6, width = 14)
par(mar=c(0,0,0,0))
layout(matrix(c(14, rep(15,10),rep(16,10),rep(13,3),
                1,rep(2,10),rep(3,10),rep(13,3),
                1,rep(2,10),rep(3,10),rep(13,3),
                1,rep(2,10),rep(3,10),rep(13,3),
                4, rep(5,10),rep(6,10),rep(13,3),
                4, rep(5,10),rep(6,10),rep(13,3),
                4, rep(5,10),rep(6,10),rep(13,3),
                7, rep(8,10),rep(9,10),rep(13,3),
                7, rep(8,10),rep(9,10),rep(13,3),
                7, rep(8,10),rep(9,10),rep(13,3),
                10, rep(11,10),rep(12,10),rep(13,3),
                10, rep(11,10),rep(12,10),rep(13,3)
                ), 12, 24, byrow=T))

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1),line=-.5, cex.axis=tick_size)
mtext("WSAF", side=4, line=-2, cex = tick_size)
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
#rect(min(tmpCoveragePos), 0, max(tmpCoveragePos), 1, col="grey", border=NA)
points(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16)
points(hapInfo$POS[hapChrom == name2], expWSAF[hapChrom == name2], col="blue", pch = 16)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), line = -.5,cex.axis=tick_size)
mtext("IBD state", side=4, line=-2, cex = tick_size)

haplotypePainter <- function (posteriorProbabilities, title = "", labelScaling,
                        numberOfInbreeding = 0){
#    rainbowColorBin <- 16
    rainbowColors = rainbowColors = c("brown1", "yellow1", "gold", "orange", "chartreuse")
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
axis(4, at=c(0,.5,1), labels=c(0,.5,1), line = -.5, cex.axis=tick_size)
mtext("Proportion", side=4, line=-2, cex = tick_size)

#pdf("haps11_and_12.pdf", height = 7, width = 18)
#par(mar=c(0,0,0,0))
#par(mfrow=c(2,1))
hap = hapInfo[hapInfo$CHROM == name1, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.axis=tick_size, yaxt="n", axes=F)
    axis(1, cex.axis = tick_size*1.1)
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
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.axis=tick_size, yaxt="n", axes=F)
    axis(1, cex.axis = tick_size*1.1)
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
plot(c(0,0),c(1,1), xlim = c(0,haplength), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("SNP Index"), side=3, line=-5.5, cex = label_size)

hap = hapInfo[hapInfo$CHROM == name2, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
plot(c(0,0),c(1,1), xlim = c(0,nhap), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("SNP Index"), side=3, line=-5.5, cex = label_size)

plot(c(0,1), c(0,1), type = "n", axes = F, xlab = "", ylab= "")

#legend("left", legend = c("No IBD", "All IBD", "1-2 IBD", "1-3 IBD", "2-3 IBD"), col = c("brown1", "chartreuse", "yellow1", "gold", "orange"), pch = 15,  bty = "n", cex = 3, ncol = 1)
legend("left", legend = c("No IBD", "All IBD", "1-2 IBD", "1-3 IBD", "2-3 IBD"), col = c("brown1", "chartreuse", "yellow1", "gold", "orange"),
    pch = 15,  bty = "n", cex = label_size*1.2, ncol = 1)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("Chromosome 11"), side=3, line=-2.5, cex = label_size)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("Chromosome 12"), side=3, line=-2.5, cex = label_size)
dev.off()

tmp.state = read.table("PD0577-C_newFilter_seed3_IBD.ibd.viterbi.txt", header=T)

pdf("hap_viterbi_hap.pdf", height = 8, width = 14)
par(mar=c(0,0,0,0))
layout(matrix(c(18, rep(19,10),rep(20,10),rep(15,3),
                1,rep(2,10),rep(3,10),rep(15,3),
                1,rep(2,10),rep(3,10),rep(15,3),
                1,rep(2,10),rep(3,10),rep(15,3),
                4, rep(5,10),rep(6,10),rep(15,3),
                4, rep(5,10),rep(6,10),rep(15,3),
                4, rep(5,10),rep(6,10),rep(15,3),
                7, rep(8,10),rep(9,10),rep(15,3),
                7, rep(8,10),rep(9,10),rep(15,3),
                7, rep(8,10),rep(9,10),rep(15,3),
                10, rep(11,10),rep(12,10),rep(15,3),
                10, rep(11,10),rep(12,10),rep(15,3),
                13, rep(14,10),rep(14,10),rep(15,3),
                13, rep(14,10),rep(14,10),rep(15,3),
                13, rep(14,10),rep(14,10),rep(15,3),
                16, rep(17,10),rep(17,10),rep(15,3),
                16, rep(17,10),rep(17,10),rep(15,3)
                ), 17, 24, byrow=T))

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1),line=-.5, cex.axis=tick_size)
mtext("WSAF", side=4, line=-2, cex = tick_size)
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
#rect(min(tmpCoveragePos), 0, max(tmpCoveragePos), 1, col="grey", border=NA)
points(tmpCoveragePos, obsWSAF[myCoverageInfo$CHROM == name2], col="red", pch = 16)
points(hapInfo$POS[hapChrom == name2], expWSAF[hapChrom == name2], col="blue", pch = 16)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), line = -.5,cex.axis=tick_size)
mtext("IBD state", side=4, line=-2, cex = tick_size)

haplotypePainter(probs[probs$CHROM==name1,c(-1,-2)], labelScaling=1)
#chromBlock = a[tmp.state$CHROM == name1,]
#barplot(t(as.matrix(chromBlock)), col =  c("brown1", "yellow1", "gold", "orange", "chartreuse"), beside=F, border=NA, space=0,xaxt = "n", yaxt = "n", axes=F)
text(c(1000,3000), c(.5, .5), labels = c(1,2), cex = tick_size)

haplotypePainter(probs[probs$CHROM==name2,c(-1,-2)], labelScaling=1)
#chromBlock = a[tmp.state$CHROM == name2,]
#barplot(t(as.matrix(chromBlock)), col =  c("brown1", "yellow1", "gold", "orange", "chartreuse"), beside=F, border=NA, space=0, xaxt = "n", yaxt = "n", axes=F)
text(c(1000,2250, 2850, 3600), c(.5, .5, .5, .5), labels = c(3,4,5,6), cex = tick_size)


plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(4, at=c(0,.5,1), labels=c(0,.5,1), line = -.5, cex.axis=tick_size)
mtext("Proportion", side=4, line=-2, cex = tick_size)

#pdf("haps11_and_12.pdf", height = 7, width = 18)
#par(mar=c(0,0,0,0))
#par(mfrow=c(2,1))
hap = hapInfo[hapInfo$CHROM == name1, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, haplength)
    yrange = c(0, 1)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.axis=tick_size, yaxt="n", axes=F)
    axis(1, cex.axis = tick_size*1.1)
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
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", xlab = "", cex.axis=tick_size, yaxt="n", axes=F)
    axis(1, cex.axis = tick_size*1.1)
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
plot(c(0,0),c(1,1), xlim = c(0,haplength), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("SNP Index"), side=3, line=-5, cex = label_size)

hap = hapInfo[hapInfo$CHROM == name2, -c(1,2)]
    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
plot(c(0,0),c(1,1), xlim = c(0,nhap), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("SNP Index"), side=3, line=-5, cex = label_size)


#13
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)

#chromBlock = a[tmp.state$CHROM == name1,]
#barplot(t(as.matrix(chromBlock)), col =  c("brown1", "yellow1", "gold", "orange", "chartreuse"), beside=F, border=NA, space=0,xaxt = "n", yaxt = "n", axes=F)


aa = tmp.state[tmp.state$CHROM %in% factor(c(as.character(name1), as.character(name2)), levels = levels(tmp.state$CHROM)),3]
segs = list()
segs[[1]] = aa[1:2018]; segs[[2]] = aa[2019:3954]; segs[[3]] = aa[3955:5996]; segs[[4]] = aa[5997:6433]; segs[[5]] = aa[6434:7161]; segs[[6]] = aa[7162:8006]


#lapply(segs, length) %>% unlist %>% sort.int(., index.return=T)
#[1]  437  728  845 1936 2018 2042
#[1] 4 5 6 2 1 3
new.aa = list()
new.aa[[1]] = segs[[4]]
new.aa[[2]] = segs[[5]]
new.aa[[3]] = segs[[6]]
new.aa[[4]] = segs[[2]]
new.aa[[5]] = segs[[1]]
new.aa[[6]] = segs[[3]]
new.aa = unlist(new.aa)

a = matrix(0.0, ncol = 5, nrow = length(new.aa))
a[cbind(1:length(new.aa), new.aa)] = 1.0

#14
barplot(t(as.matrix(a)), col =  c("brown1", "yellow1", "gold", "orange", "chartreuse"), beside=F, border=NA, space=0,xaxt = "n", yaxt = "n", axes=F)
bb = lapply(segs, length) %>% unlist %>% sort.int(., index.return=T)
bbx=bb$x / 2
text(cumsum(c(0,bb$x[1:5]))+bbx, 0.5, labels = c("4", "5", "6", "2", "1 (N50)", "3"), cex = tick_size)

# 15
plot(c(0,1), c(0,1), type = "n", axes = F, xlab = "", ylab= "")
#legend("left", legend = c("No IBD", "All IBD", "1-2 IBD", "1-3 IBD", "2-3 IBD"), col = c("brown1", "chartreuse", "yellow1", "gold", "orange"), pch = 15,  bty = "n", cex = 3, ncol = 1)
legend("left", legend = c("No IBD", "All IBD", "1-2 IBD", "1-3 IBD", "2-3 IBD"), col = c("brown1", "chartreuse", "yellow1", "gold", "orange"),
    pch = 15,  bty = "n", cex = label_size*1.2, ncol = 1)



# 16 17
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
plot(c(0,0),c(1,1), xlim = c(0,100), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
axis(side = 3, at = c(0,50,100), labels= c("0", "50", "100%"), cex.axis = tick_size)
mtext(paste("Ordered IBD blocks"), side=3, line=-4.5, cex = label_size)

# 18 19 20
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("Chromosome 11"), side=3, line=-2.5, cex = label_size)
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(paste("Chromosome 12"), side=3, line=-2.5, cex = label_size)
dev.off()



