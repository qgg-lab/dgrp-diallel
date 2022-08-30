# ==================================================
# = calculate relationship matrix for the 50 lines =
# ==================================================

library("rrBLUP")

read.tped <- function(file.prefix) {
  # check file existence
  if (!file.exists(paste(file.prefix, ".tped", sep = "")) || !file.exists(paste(file.prefix, ".tfam", sep = ""))) {
    stop("cannot find tped/tfam files!")
  }
  
  # read genotypes
  tfam <- read.table(paste(file.prefix, ".tfam", sep = ""), as.is = TRUE, header = FALSE)
  n.ind <- nrow(tfam)
  tped <- matrix(scan(paste(file.prefix, ".tped", sep = ""), what = ""), ncol = 4 + 2*n.ind, byrow = TRUE)

  snp.name <- tped[, 2]
  snp.loc <- tped[, c(1, 4)]
  geno.code <- matrix(4 - (as.numeric(tped[, seq(from = 5, length = n.ind, by = 2)]) + as.numeric(tped[, seq(from = 6, length = n.ind, by = 2)])), ncol = n.ind)
  geno.code[geno.code == 4] <- NA
  
  rownames(geno.code) <- snp.name
  colnames(geno.code) <- tfam[, 1]
  
  return(geno.code)
  
}

geno <- read.tped("gwas/diallel")

# impute genotypes
# ============================================================

for (i in 1:nrow(geno)) {
  if (sum(is.na(geno[i, ])) > 0) {
    this.geno <- geno[i, ]
    this.geno[is.na(this.geno)] <- mean(this.geno, na.rm = TRUE)
    geno[i, ] <- this.geno
  }
  if (i %% 10000 == 0) {
    cat(i, "\n")
  }
}

X <- scale(t(geno))
A <- X %*% t(X)
A <- A/mean(diag(A))

# perform GWAS
# ============================================================

pheno <- read.table("gwas/female.pheno", header = FALSE, as.is = TRUE)[, -1]
pheno <- pheno[match(colnames(geno), pheno[, 1]), ]
mixed.fit <- mixed.solve(pheno[, 2], Z = diag(50), K = A)
mixed.fit

# get PCs
# ============================================================

pca <- prcomp(X)
pcs <- predict(pca)
fit <- lm(pheno[, 2] ~ pcs[, 1:10])
summary(fit)
# none is significant

# write out A's
write.table(rbind(c("var", paste(colnames(A), colnames(A), sep = " ")), cbind(paste(colnames(A), colnames(A), sep = " "), A)), file = "gwas/A.rel.mat", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
