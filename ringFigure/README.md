use PD0577-C IBD seed 3

scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/refAndAlt/PD0577-C_ref.txt .
scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/refAndAlt/PD0577-C_alt.txt .
scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/asiaGroup1_PLAF.txt .

scp rescomp:/well/mcvean/joezhu/pf3k-IBD-deconvolution/asiaGroup1/PD0577-C/*seed3_IBD* .

use the new filter ones

~/pfPvRecombMap/R/plotIBD.r PD0577-C_newFilter_seed3_IBD  asiaGroup1_PLAF.txt PD0577-C_ref.txt PD0577-C_alt.txt asiaGroup1_and_pf3k_bad_snp_in_at_least_50_samples.txt

making linear plots for tracks
