use PD0577-C IBD seed 3

scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/refAndAlt/PD0577-C_ref.txt .
scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/refAndAlt/PD0577-C_alt.txt .
scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/asiaGroup1_PLAF.txt .

scp rescomp:/well/mcvean/joezhu/pf3k-IBD-deconvolution/asiaGroup1/PD0577-C/*seed3_IBD* .

use the new filter ones

~/pfPvRecombMap/R/plotIBD.r PD0577-C_newFilter_seed3_IBD  asiaGroup1_PLAF.txt PD0577-C_ref.txt PD0577-C_alt.txt asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt

making linear plots for tracks


dEploid -plaf asiaGroup1_PLAF.txt -ref PD0577-C_ref.txt -alt PD0577-C_alt.txt -exclude asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt -noPanel -ibdPainting -o post -initialP .22 .52 .26 -o ibdpost-Gdefault
dEploid -plaf asiaGroup1_PLAF.txt -ref PD0577-C_ref.txt -alt PD0577-C_alt.txt -exclude asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt -noPanel -ibdViterbi -o viterbi -initialP .22 .52 .26 -o viterbi-Gdefault

dEploid -plaf asiaGroup1_PLAF.txt -ref PD0577-C_ref.txt -alt PD0577-C_alt.txt -exclude asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt -noPanel -ibdPainting -o post -initialP .22 .52 .26 -o ibdpost-G0dot1 -G 0.1
dEploid -plaf asiaGroup1_PLAF.txt -ref PD0577-C_ref.txt -alt PD0577-C_alt.txt -exclude asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt -noPanel -ibdViterbi -o viterbi -initialP .22 .52 .26 -o viterbi-G0dot1 -G 0.1

dEploid -plaf asiaGroup1_PLAF.txt -ref PD0577-C_ref.txt -alt PD0577-C_alt.txt -exclude asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt -noPanel -ibdPainting -o post -initialP .22 .52 .26 -o ibdpost-G0dot05 -G 0.05
dEploid -plaf asiaGroup1_PLAF.txt -ref PD0577-C_ref.txt -alt PD0577-C_alt.txt -exclude asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt -noPanel -ibdViterbi -o viterbi -initialP .22 .52 .26 -o viterbi-G0dot05 -G 0.05

rm(list= ls())
suffix = "-Gdefault"


#suffix = "-G0dot1"
suffix = "-G0dot05"

postFile = paste("ibdpost", suffix, ".ibd.probs", sep = "")
post = read.table(postFile, header=T, stringsAsFactors = F)

tmpFile = paste("viterbi", suffix, ".viterbi", sep = "")
tmp.state = read.table(tmpFile, header=T, stringsAsFactors = F)

    nstate = 2
    stateSpace =  unique(tmp.state[,3])
    if (length(stateSpace) != 2){
        nstate = 5
    }

    a = matrix(0.0, ncol = nstate, nrow = length(tmp.state[, 3]))
    a[cbind(1:length(tmp.state[,3]), tmp.state[,3]+1)] = 1.0


    viterbi.ibd.only.blocks = c()

png(paste("PD0577-C", suffix, ".post.png", sep = ""), width = 2000, height = 1200)
    par(mfrow = c(7,2))
        for ( chrom in unique(tmp.state$CHROM) ){
            chromBlock = post[post$CHROM == chrom,-c(1,2)]
            barplot(t(as.matrix(chromBlock)), col =  c("brown1", "orange", "gold", "chartreuse", "yellow1"), beside=F, border=NA, space=0)
        }

    dev.off()

png(paste("PD0577-C", suffix, ".viterbi.png", sep = ""), width = 2000, height = 1200)
    par(mfrow = c(7,2))
        for ( chrom in unique(tmp.state$CHROM) ){
            chromBlock = a[tmp.state$CHROM == chrom, ]
            barplot(t(as.matrix(chromBlock)), col =  c("brown1", "orange", "gold", "chartreuse", "yellow1"), beside=F, border=NA, space=0)
        }

    dev.off()
