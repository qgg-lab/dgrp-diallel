# =======================================================
# = analyze and make figures for the body size analysis =
# =======================================================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "reportData/female.pheno", "bodySize.csv", "reportData/top.snp.geno.tped", "reportData/top.snp.geno.tfam")
library("RColorBrewer")

# read data
# ============================================================

female.pheno <- read.table(args[2], header = FALSE, as.is = TRUE)
rownames(female.pheno) <- female.pheno[, 1]

body.size <- read.csv(args[3], header = TRUE, as.is = TRUE)
thorax.length <- with(subset(body.size, Sex == 1), sapply(split(Length.mm., paste("line_", Line, sep = "")), mean, na.rm = TRUE))
thorax.width <- with(subset(body.size, Sex == 1), sapply(split(Width.mm., paste("line_", Line, sep = "")), mean, na.rm = TRUE))
thorax.length <- thorax.length[female.pheno[, 1]]
thorax.width <- thorax.width[female.pheno[, 1]]

tped <- read.table(args[4], header = FALSE, as.is = TRUE, na.strings = "0")
tfam <- read.table(args[5], header = FALSE, as.is = TRUE)
geno <- tped[, seq(5, ncol(tped), 2)]
colnames(geno) <- tfam[, 1]
rownames(geno) <- tped[, 2]

d2r.geno <- unlist(geno["X_19897101", female.pheno[, 1]])

diallel.female <- 

# make plot of 
# ============================================================

# test
# ============================================================

pval <- matrix(ncol = 3, nrow = nrow(geno))

for (i in 1:nrow(geno)) {
  
  pval[i, ] <- c(anova(lm(female.pheno[, 3] ~ as.factor(unlist(geno[i, ]))))$Pr[1],
                 anova(lm(female.pheno[, 3] ~ thorax.length + as.factor(unlist(geno[i, ]))))$Pr[2], 
                 anova(lm(thorax.length ~ as.factor(unlist(geno[i, ]))))$Pr[1])
  
}