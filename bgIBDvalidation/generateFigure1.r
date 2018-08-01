rm(list=ls())
library(dplyr)
#sample pair of PG0055-C PG0054-C
source("src.r")
source("~/DEploid/utilities/dEploidTools.r")
computeNx <- function (x, threshold){
    if (length(x) > 0 ){
        sort.x = sort(x)
        sum.x = sum(x)
        cumsum.x = cumsum(sort.x/sum.x)
        idx = which(cumsum.x > threshold)[1]
        return (sort.x[idx])
    } else {
        return (NA)
    }
}

sample1 = "PG0058-C"
sample2 = "PG0071-C"

sample1.single0 = read.table(paste("/home/joezhu/recombMapPaper/figures/figure_validation_by_crosses/paintings/", sample1, ".single0", sep = ""), header=T, stringsAsFactors=F)
sample2.single0 = read.table(paste("/home/joezhu/recombMapPaper/figures/figure_validation_by_crosses/paintings/", sample2, ".single0", sep = ""), header=T, stringsAsFactors=F)
post.chrom.pos = paste(sample1.single0$CHROM, sample1.single0$POS, sep = "-")

sample1.sample2.viterbi = read.table(paste("/home/joezhu/pf3k-IBD-deconvolution-clean-version/crosses/PG0058-C-PG0071-C.wg.ibd.viterbi.txt", sep =""), header=T, stringsAsFactors=F)
viterbi.chrom.pos = paste(sample1.sample2.viterbi$CHROM, sample1.sample2.viterbi$POS, sep = "-")

post.keeping.idx = which(post.chrom.pos %in% viterbi.chrom.pos)

threshold = .5
all.together = data.frame(CHROM = sample1.sample2.viterbi$CHROM,
                          POS = sample1.sample2.viterbi$POS,
                          sample1 = sample1.single0[post.keeping.idx,3]> threshold,
                          sample2 = sample2.single0[post.keeping.idx,3]> threshold,
                          painting = ((sample1.single0[post.keeping.idx,3] > threshold) == (sample2.single0[post.keeping.idx,3] > threshold))*1,
                          viterbi = sample1.sample2.viterbi[,3]
                          )

tick_size = 1.5*1.2
label_size = 2


plot.at = round(seq(1, 1056812, length.out = 1000))
plotting.mat = all.together[plot.at, -c(2)]

pdf("bgIBD.pdf", width = 14, height = 8)

layout(matrix(c(1, rep(2,20), rep(10,5),
                1, rep(2,20), rep(10,5),
                3, rep(4,20), rep(10,5),
                3, rep(4,20), rep(10,5),
                5, rep(6,20), rep(10,5),
                5, rep(6,20), rep(10,5),
                7, rep(8,20), rep(10,5),
                7, rep(8,20), rep(10,5),
                11, rep(9,20), rep(10, 5)), 9, 26, byrow=T))

par(mar=c(0,0,0,0))
plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(sample1, side=4, line=-2, cex = tick_size)
barplot(t(cbind(plotting.mat[,2], 1-plotting.mat[,2])), col = c("red", "blue"), space = 0, border=F, axes=F)
abline(v = which(diff(as.numeric(plotting.mat[,1])) > 0), lwd = 2 )


plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext(sample2, side=4, line=-2, cex = tick_size)
barplot(t(cbind(plotting.mat[,3], 1-plotting.mat[,3])), col = c("red", "blue"), space = 0, border=F, axes=F)
abline(v = which(diff(as.numeric(plotting.mat[,1])) > 0), lwd = 2 )

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext("IBD states", side=4, line=-2, cex = tick_size)
barplot(t(plotting.mat[,4]), col = c("grey80", "white"),space = 0, border=F, axes=F)
abline(v = which(diff(as.numeric(plotting.mat[,1])) > 0), lwd = 2 )

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
mtext("IBD states", side=4, line=-2, cex = tick_size)
barplot(t(plotting.mat[,5]-1),col = c("grey50", "white"), space = 0, border=F, axes=F)
abline(v = which(diff(as.numeric(plotting.mat[,1])) > 0), lwd = 2 )

#sample1.sample2.viterbi = read.table(paste("viterbis/", sample1, "-", sample2, ".wg.ibd.viterbi.txt", sep =""), header=T, stringsAsFactors=F)
#sample1.sample2.viterbi = read.table(paste("viterbis/snp20k.filter/", sample1, "-", sample2, "-20k.filter.ibd.viterbi.txt", sep =""), header=T, stringsAsFactors=F)

plot(c(0,1000),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
tmp.tab = table(plotting.mat[,1])
text(cumsum(tmp.tab) - round( tmp.tab / 2), .5, labels = 1:14, cex = 2)

plot(c(0,0),c(1,1), ylim = c(0,1), type="n", ylab="", xlab="", yaxt="n", xaxt="n", axes=F)
legend("right", legend = c("3D7 ancestry", "HB3 ancestry", "Compare ancestry", "DEploidIBD", "non-IBD"), cex = 2, bty = "n",
                   fill = c("red", "blue", "grey80", "grey50", "white"))
dev.off()

computeNx <- function (x, threshold){
    if (length(x) > 0 ){
        sort.x = sort(x)
        sum.x = sum(x)
        cumsum.x = cumsum(sort.x/sum.x)
        idx = which(cumsum.x > threshold)[1]
        return (sort.x[idx])
    } else {
        return (NA)
    }
}


block.length.by.posterior.all = c()
block.length.by.posterior.chrom = c()

block.length.by.viterbi.all = c()
block.length.by.viterbi.chrom = c()

for ( chrom in unique(sample1.single0$CHROM)){

#    png(paste(sample1, "-", sample2, "-", chrom, ".png", sep =""), width = 2000, height = 600)
    chromBlock_i = 1*(sample1.single0[sample1.single0$CHROM == chrom,3] > threshold)
    chromBlock_j = 1*(sample2.single0[sample2.single0$CHROM == chrom,3] > threshold)

    blocks = exact_chunk_length_vec(sample1.single0[sample1.single0$CHROM == chrom,2], 1*(chromBlock_i==chromBlock_j),T)
    block.length.by.posterior.chrom = c(block.length.by.posterior.chrom, computeNx(blocks, .5))
    block.length.by.posterior.all = c(block.length.by.posterior.all, blocks)

    subseq = round(seq(1, length(chromBlock_i), length.out = 1000))
#    par(mfrow = c(4,1))
#    barplot(t(sample1.single0[sample1.single0$CHROM == chrom,-c(1,2)][subseq,]), main = paste(sample1, "posterior prob"), col = c("red", "blue"), border=F)
#    barplot(t(sample2.single0[sample2.single0$CHROM == chrom,-c(1,2)][subseq,]), main = paste(sample2, "posterior prob"), col = c("red", "blue"), border=F)
#    barplot(t(chromBlock_j == chromBlock_i), main = paste(sample1, sample2, "IBD"))

    viterbi.block = sample1.sample2.viterbi[sample1.sample2.viterbi$CHROM == chrom, ]
#    barplot(t(viterbi.block[,3]-1), main = paste(sample1, sample2, "viterbi"))

    block.viterbi = exact_chunk_length_vec(viterbi.block$POS, 1*(viterbi.block[,3] == 2))
    block.length.by.viterbi.all = c(block.length.by.viterbi.all, block.viterbi)
    block.length.by.viterbi.chrom = c(block.length.by.viterbi.chrom, computeNx( block.viterbi, .5))

#    dev.off()
}

#hmmibd = read.table("tmp.hmm.txt", header=T)
#hmmibd_len = hmmibd$end[hmmibd$different==0] - hmmibd$start[hmmibd$different==0]

pdf("bgN50.pdf", width = 8, height = 8)
par(mar = c(5,5,5,3))
x = lapply(seq(0, .9, by = .1), function(x){computeNx(block.length.by.posterior.all, x)}) %>% unlist
y = lapply(seq(0, .9, by = .1), function(x){computeNx(block.length.by.viterbi.all, x)}) %>% unlist
#hmmibd_len_Nx = lapply(seq(0, .9, by = .1), function(x){computeNx(hmmibd_len, x)}) %>% unlist
plot(x/1000,y/1000, xlim=c(0, 1.2e3), ylim=c(0,1.2e3), xlab = "Compare ancestry IBD block NX (kb)", ylab = "DEploidIBD IBD block NX (kb)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, main = paste("Correlation: ", round(cor(x,y), digits = 2)))
#points(x/1000, hmmibd_len_Nx/1000, col = "red")
abline(0,1)
dev.off()









