rm(list=ls()); dEploidRootDir="/home/joezhu/DEploid"
# DESCRIPTION:
#
# USAGE:
#    ./interpretDEploid.r -vcf FILE -plaf FILE -dEprefix STRING -o STRING
#    R --slave "--args -vcf FILE -plaf FILE -dEprefix STRING -o STRING " < utilities/interpretDEploid.r
#
# EXAMPLE:
#    ./interpretDEploid.r -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel " < utilities/interpretDEploid.r
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanelExclude -o PG0390-CNopanelExclude -exclude data/testData/labStrains.test.exclude.txt " < utilities/interpretDEploid.r

if (!exists("dEploidRootDir")){
print("dEploidRootDir undefined, try make dEploid again!")
}

library(methods)
source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))
source(paste(dEploidRootDir, "/utilities/DEploidR.R", sep=""))

#args = (commandArgs(TRUE))

args = c("-ref", "PD0577-C_ref.txt",
         "-alt", "PD0577-C_alt.txt",
         "-exclude", "asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt",
         "-dEprefix", "PD0577-C_newFilter_seed3_IBD",
         "-o", "myring")

myInput = fun.parse ( args )




myCoverageInfo = fun.extract.coverage ( myInput )

myExcludeInfo = fun.extract.exclude (myInput$excludeFileName, myInput$excludeBool)

plot.wsaf.vs.index.ring <- function ( coverage, expWSAF = c(), expWSAFChrom = c(), exclude, titlePrefix = "" ){
    chromCol = (as.numeric(1:length(levels(coverage$CHROM)) %% 2 ))
    chromCol[chromCol==1] = NA
    chromCol[chromCol==0] = 8

    circlize::circos.trackPlotRegion(factor = expWSAFChrom, ylim=c(0,1), track.height = 0.18, bg.col = chromCol, panel.fun=function(x,y){
        name = circlize::get.cell.meta.data("sector.index")
        xlim = circlize::get.cell.meta.data("xlim")
        ylim = circlize::get.cell.meta.data("ylim")
        chromRegion = coverage[coverage$CHROM==name,]

        ref = chromRegion$refCount
        alt = chromRegion$altCount
        obsWSAF = computeObsWSAF ( alt, ref )

        nSnp = dim(chromRegion)[1]

        if (exclude$excludeBool){
            tmpCoveragePos = coverage$POS[coverage$CHROM==name]
            tmpExcludePos = exclude$excludeTable$POS[exclude$excludeTable$CHROM==name]
            excludeLogic = ( tmpCoveragePos %in% tmpExcludePos )
            plotIndex = which(!excludeLogic)
        } else {
            plotIndex = c(1:nSnp)
        }
        circlize::circos.points(1:length(plotIndex), obsWSAF[plotIndex], col="red", pch = 16)
        circlize::circos.points(1:length(plotIndex), expWSAF[expWSAFChrom == name], col="blue", pch = 16)
    })
}


png(paste(myInput$outPrefix, ".ring.png", sep = ""), width = 3500, height = 3500, bg = "transparent")
    cexSize = 2.5

dEploidOutput = fun.dEploidPrefix ( myInput$dEploidPrefix )
tmpProp = read.table(dEploidOutput$propFileName, header=F)

lastProp = as.numeric(tmpProp[dim(tmpProp)[1],])

probs = read.table("PD0577-C_newFilter_seed3_IBD.ibd.postProb.txt", header=T)

fun.ring.plot.initialize(probs$CHROM)
    hapInfo = read.table(dEploidOutput$hapFileName, header=T)
    hapChrom = hapInfo[,1]
    hap = as.matrix(hapInfo[,-c(1,2)])
    expWSAF = hap %*% lastProp

plot.wsaf.vs.index.ring ( myCoverageInfo, expWSAF, hapChrom, myExcludeInfo )

circlize::circos.trackPlotRegion(factor = probs$CHROM, ylim=c(0,1), track.height = 0.3,
    panel.fun=function(x,y){
        name = circlize::get.cell.meta.data("sector.index")
        xlim = circlize::get.cell.meta.data("xlim")
        ylim = circlize::get.cell.meta.data("ylim")
        cat(".")
        chromRegion = probs[probs$CHROM==name,]
            rainbowColorBin <- 16
            rainbowColors = rainbow(rainbowColorBin)
        nSnp = dim(chromRegion)[1]
        print(nSnp)
        nhap = dim(chromRegion)[2] - 2
        cumProb = rep(0, nSnp)
        for ( i in 1:nhap ){
            circlize::circos.rect(0:(nSnp-1), cumProb, 1:nSnp, cumProb + chromRegion[,i+2], col=rainbowColors[i], border=NA)
            cumProb = cumProb + chromRegion[,i+2]
        }
    }
)

circlize::circos.clear();
dev.off()

system("convert myring.ring.png -quality 100 myring.ring.pdf")
