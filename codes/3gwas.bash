# ==============================
# = run gwas and annotate SNPs =
# ==============================

# 1. prepare phenotype
# ============================================================

~/software/R-3.2.2/bin/Rscript -e 'ind.data <- read.csv("diallelWolbaAdjusted.csv", header = FALSE, as.is = TRUE); line.means <- sapply(split(ind.data[, 4], ind.data[, 1]), mean); write.table(cbind(names(line.means), names(line.means), line.means), file = "gwas/female.pheno", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ");'

# 2. test for inversion effects
#    nothing is significant
# ============================================================

~/software/R-3.2.2/bin/Rscript inversionTest.R > gwas/inversionTest.Rout

# 3. simple run without any adjustment
# ============================================================

~/software/plink-v1.90b3w/plink --silent --bfile jgil/diallel --pheno gwas/female.pheno --assoc --out gwas/diallel.gwas

# 4. get top snps and their frequency information
# ============================================================

awk '$9 < 1e-5 { print $2 }' gwas/diallel.gwas.qassoc > gwas/top.snp
~/software/plink-v1.90b3w/plink --silent --bfile jgil/diallel --extract gwas/top.snp --freq --out gwas/top.snp

# 5. get snp genotypes in all DGRP lines
# ============================================================

~/software/plink-v1.90b3w/plink --silent --bfile jgil/all --extract gwas/top.snp --recode --transpose --out gwas/top.snp.geno

# 6. pleiotropy analysis
# ============================================================

# create phenotype data on local computer
Rscript reportCode/adjustTransform.R reportData/line.id.txt reportData/adjustData.RData reportData/female.line.means.pheno reportData/male.line.means.pheno reportData/both.line.means.pheno

# get SNPs to test, test only those with similar MAF X_19897101
# ============================================================

~/qgg/software/plink-v1.90b5.3/plink --silent --tped ../jgil/all.tped --tfam ../jgil/all.tfam --keep dgrp.fam --make-bed --out dgrp &
~/qgg/software/plink-v1.90b5.3/plink --silent --bfile dgrp --freq --out dgrp.freq &
awk '$5 >= 0.275 && $5 <= 0.325 {print $2}' dgrp.freq.frq > dgrp.common.snp
~/qgg/software/plink-v1.90b5.3/plink --silent --bfile dgrp --extract dgrp.common.snp --make-bed --out dgrp.common &

~/qgg/software/plink-v1.90b5.3/plink --silent --bfile dgrp.common --pheno both.line.means.pheno --all-pheno --assoc --out both
~/qgg/software/plink-v1.90b5.3/plink --silent --bfile dgrp.common --pheno female.line.means.pheno --all-pheno --assoc --out female
~/qgg/software/plink-v1.90b5.3/plink --silent --bfile dgrp.common --pheno male.line.means.pheno --all-pheno --assoc --out male

