# ================================
# = analyze the D2R diallel data =
# ================================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "reportData/prod.csv", "reportData/female.pheno", "reportData/top.snp.geno.tped", "reportData/top.snp.geno.tfam", "report/figureD2Rdiallel.pdf")
library("RColorBrewer")

# prepare file
# ============================================================

file.width = 45 # in mm
cairo_pdf(file = args[6], width = file.width/25.4, height = file.width/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(1.5, 2, 0.8, 0.2), ps = 7, lwd = 0.5, xpd = TRUE)

# read data
# ============================================================

diallel.raw.data <- read.csv(args[2], header = TRUE, as.is = TRUE, na.strings = "")
diallel.raw.data[is.na(diallel.raw.data)] <- as.integer(0)

# remove data when there are female or male dead
diallel.prod.data <- subset(diallel.raw.data, (d1_female_dead_at_transfer == 0 & d1_male_dead_at_transfer == 0 &
                                               d2_female_dead_at_transfer == 0 & d2_male_dead_at_transfer == 0))

diallel.prod.data <- within(diallel.prod.data, { female <- paste("line_", female, sep = ""); male <- paste("line_", male, sep = ""); })

# calculate total number of eggs and adults
diallel.prod.data <- within(diallel.prod.data, { total_eggs <- d1_eggs + d2_eggs;
                                                 total_adults <- d1_d10female + d1_d10male + d1_d13female + d1_d13male +
                                                 d1_d16female + d1_d16male + d1_dead_female + d1_dead_male +
                                                 d2_d10female + d2_d10male + d2_d13female + d2_d13male +
                                                 d2_d16female + d2_d16male + d2_dead_female + d2_dead_male; })

diallel.prod.data <- subset(diallel.prod.data, (total_adults <= total_eggs))


female.pheno <- read.table(args[3], header = FALSE, as.is = TRUE)
rownames(female.pheno) <- female.pheno[, 1]

tped <- read.table(args[4], header = FALSE, as.is = TRUE, na.strings = "0")
tfam <- read.table(args[5], header = FALSE, as.is = TRUE)
geno <- tped[, seq(5, ncol(tped), 2)]
colnames(geno) <- tfam[, 1]
rownames(geno) <- tped[, 2]

d2r.geno <- unlist(geno["X_19897101", female.pheno[, 1]])

diallel.female.egg.mean <- with(diallel.prod.data, sapply(split(total_eggs, female), mean))
diallel.female.adult.mean <- with(diallel.prod.data, sapply(split(total_adults, female), mean))

diallel.female.egg.se <- with(diallel.prod.data, sapply(split(total_eggs, female), function(x) {return(sqrt(var(x)/length(x)))}))
diallel.female.adult.se <- with(diallel.prod.data, sapply(split(total_adults, female), function(x) {return(sqrt(var(x)/length(x)))}))

diallel.d2r.geno <- unlist(geno["X_19897101", names(diallel.female.egg.mean)])

# plot data
# ============================================================

plot.col <- rep(brewer.pal(9, "Reds")[9], 10)
plot.col[diallel.d2r.geno == "T"] <- brewer.pal(9, "Blues")[9]
plot.col[names(diallel.female.egg.mean) == "line_427"] = brewer.pal(9, "Blues")[4]
plot(c(rep(1, 10), rep(2, 10)), c(diallel.female.egg.mean, diallel.female.adult.mean), xlim = c(0.5, 2.5), ylim = c(0, 200), pch = 21, bg = c(plot.col, plot.col), col = NA, axes = FALSE, xlab = "", ylab = "", cex = 0.5)
segments(rep(1, 10), diallel.female.egg.mean, rep(2, 10), diallel.female.adult.mean, col = plot.col)

axis(side = 1, at = c(1, 2), labels = c("Eggs", "Adult offspring"), lwd = 0.5, mgp = c(0.8, 0, 0), cex.axis = 6/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.4, 0), cex.axis = 6/par("ps")/par("cex"))
box(bty = "l")
title(ylab = "Count", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

legend("topright", pch = 21, lty = 1, pt.bg = c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9]), col = c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9]), legend = c("T", "G"), cex = 0.5, bty = "n", x.intersp = 0.5, y.intersp = 0.8, title = expression(paste(italic("D2R"), " SNP")))

dev.off()
