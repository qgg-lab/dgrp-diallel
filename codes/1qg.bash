# =================================
# = quantitative genetic analysis =
# =================================

# 1. QG analysis to test for wolbachia effet
# ============================================================

~/software/R-3.2.2/bin/Rscript prepWolba.R diallel.csv wolba.csv sas/wolbaData.csv
sas sas/wolba.sas -print sas/wolba.lst -log sas/wolba.log &

# 2. QG analysis after adjustment of wolbachia
#    test for cross and the cockerham's bio-model
# ============================================================

~/software/R-3.2.2/bin/Rscript adjustWolba.R diallel.csv wolba.csv 5.9080 6.0108 3.6544 diallelWolbaAdjusted.csv
sas sas/varComp.sas -print sas/varComp.lst -log sas/varComp.log &
