# ============================================================================================
# = adjust line means for covariates (inversions, wolbachia, PCs), normal quantile transform =
# ============================================================================================

args <- commandArgs(TRUE)
# args <- c("reportData/line.id.txt", "reportData/adjustData.RData", "reportData/female.line.means.pheno", "reportData/male.line.means.pheno")
load(args[2])

# function to do normal quantile transformation
# ============================================================

nqt <- function(x) {
  
  return(qnorm(rank(x)/(length(x) + 1), sd = 1))
  
}

# function to do adjustment and 
adjustPheno <- function(raw.pheno) {

# raw.pheno is a two column data.frame where the first is
# line id and the second is raw phenotype

  common.line <- intersect(rownames(wolba), raw.pheno[, 1])
  rownames(raw.pheno) <- raw.pheno[, 1]
  # get data
  pheno.data <- cbind(raw.pheno[common.line, 2], wolba[common.line, 1], inv[common.line, c("In_2L_t", "In_2R_NS", "In_3R_P", "In_3R_K", "In_3R_Mo")], pcs[common.line, ])
  colnames(pheno.data)[2] <- "wolba"
  fit.form <- "pheno.data[, 1] ~ 1"
  if (length(unique(pheno.data[, "wolba"])) > 1) {
    fit.form <- paste(fit.form, " + factor(wolba)")
  }
  if (length(unique(pheno.data[, "In_2L_t"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_2L_t)")
  }
  if (length(unique(pheno.data[, "In_2R_NS"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_2R_NS)")
  }
  if (length(unique(pheno.data[, "In_3R_P"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_P)")
  }
  if (length(unique(pheno.data[, "In_3R_K"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_K)")
  }
  if (length(unique(pheno.data[, "In_3R_Mo"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_Mo)")
  }
  fit.form <- paste(fit.form, " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11", sep = "")

  lm.fit <- lm(formula(fit.form), data = pheno.data)
  
  covar.coef <- summary(lm.fit)$coefficients
  
  pheno.trans <- nqt(residuals(lm.fit))
  
  return(list(covar.coef, pheno.trans))
  
}

# main codes to process phenotypes
# ============================================================

# profile the directory
# ============================================================

line.id <- scan(file = args[1], what = "", quiet = TRUE)

female.trait.files <- list.files(path = "reportData", pattern = "line.means.*\\.female\\.csv")
male.trait.files <- list.files(path = "reportData", pattern = "line.means.*\\.male\\.csv")
both.trait.files <- list.files(path = "reportData", pattern = "line.means.*\\.both\\.csv")

female.traits <- gsub("line\\.means\\.", "", gsub("\\.female\\.csv", "", female.trait.files))
male.traits <- gsub("line\\.means\\.", "", gsub("\\.male\\.csv", "", male.trait.files))
both.traits <- gsub("line\\.means\\.", "", gsub("\\.both\\.csv", "", both.trait.files))


# loop through traits to load data
# ============================================================

female.line.means <- matrix(nrow = length(line.id), ncol = length(female.traits))
for (i in 1:length(female.traits)) {
  trait.name <- female.traits[i]
  # read data
  female.line.means[, i] <- read.csv(paste("reportData/line.means.", trait.name, ".female.csv", sep = ""), header = FALSE, as.is = TRUE, row.names = 1)[line.id, 1]
}
colnames(female.line.means) <- female.traits

male.line.means <- matrix(nrow = length(line.id), ncol = length(male.traits))
for (i in 1:length(male.traits)) {
  trait.name <- male.traits[i]
  # read data
  male.line.means[, i] <- read.csv(paste("reportData/line.means.", trait.name, ".male.csv", sep = ""), header = FALSE, as.is = TRUE, row.names = 1)[line.id, 1]
}
colnames(male.line.means) <- male.traits

both.line.means <- matrix(nrow = length(line.id), ncol = length(both.traits))
for (i in 1:length(both.traits)) {
  trait.name <- both.traits[i]
  # read data
  both.line.means[, i] <- read.csv(paste("reportData/line.means.", trait.name, ".both.csv", sep = ""), header = FALSE, as.is = TRUE, row.names = 1)[line.id, 1]
}
colnames(both.line.means) <- both.traits


# make adjustment
# ============================================================

female.adjust.pheno <- female.line.means
female.adjust.coef <- list()
for (i in 1:ncol(female.line.means)) {
  
  this.trait.result <- adjustPheno(data.frame(line <- as.character(line.id), pheno <- female.line.means[, i], stringsAsFactors = F))
  female.adjust.coef[[i]] <- this.trait.result[[1]]
  female.adjust.pheno[, i] <- this.trait.result[[2]][line.id]
  
}
colnames(female.adjust.pheno) <- paste("female.", female.traits, sep = "")
names(female.adjust.coef) <- paste("female.", female.traits, sep = "")

male.adjust.pheno <- male.line.means
male.adjust.coef <- list()
for (i in 1:ncol(male.line.means)) {
  
  this.trait.result <- adjustPheno(data.frame(line <- as.character(line.id), pheno <- male.line.means[, i], stringsAsFactors = F))
  male.adjust.coef[[i]] <- this.trait.result[[1]]
  male.adjust.pheno[, i] <- this.trait.result[[2]][line.id]
  
}
colnames(male.adjust.pheno) <- paste("male.", male.traits, sep = "")
names(male.adjust.coef) <- paste("male.", male.traits, sep = "")

both.adjust.pheno <- both.line.means
both.adjust.coef <- list()
for (i in 1:ncol(both.line.means)) {
  
  this.trait.result <- adjustPheno(data.frame(line <- as.character(line.id), pheno <- both.line.means[, i], stringsAsFactors = F))
  both.adjust.coef[[i]] <- this.trait.result[[1]]
  both.adjust.pheno[, i] <- this.trait.result[[2]][line.id]
  
}
colnames(both.adjust.pheno) <- paste("both.", both.traits, sep = "")
names(both.adjust.coef) <- paste("both.", both.traits, sep = "")



# write data
# ============================================================

write.table(rbind(c("FID", "IID", colnames(female.adjust.pheno)), cbind(line.id, line.id, female.adjust.pheno)), file = args[3], sep = " ", col.names = F, row.names = F, quote = F)
write.table(rbind(c("FID", "IID", colnames(male.adjust.pheno)), cbind(line.id, line.id, male.adjust.pheno)), file = args[4], sep = " ", col.names = F, row.names = F, quote = F)
write.table(rbind(c("FID", "IID", colnames(both.adjust.pheno)), cbind(line.id, line.id, both.adjust.pheno)), file = args[5], sep = " ", col.names = F, row.names = F, quote = F)
