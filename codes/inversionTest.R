# ====================================================
# = test for inversion effects for female line means =
# ====================================================

# load phenotype data
# ============================================================

pheno <- read.table("gwas/female.pheno", header = FALSE, as.is = TRUE)[, -1]

# inversions
# ============================================================

load("~/dgrp/adjustData.RData")

# test
# ============================================================

common.line <- intersect(rownames(wolba), pheno[, 1])
rownames(pheno) <- pheno[, 1]
pheno.data <- cbind(pheno[common.line, 2], inv[common.line, c("In_2L_t", "In_2R_NS", "In_3R_P", "In_3R_K", "In_3R_Mo")])
fit <- lm(pheno.data[, 1] ~ as.factor(pheno.data$In_2L_t) + as.factor(pheno.data$In_2R_NS) + as.factor(pheno.data$In_3R_P) + as.factor(pheno.data$In_3R_Mo))
drop1(fit, test = "F")

# none is signifcant
# ============================================================
