# =========================================
# = check productivity for all DGRP lines =
# =========================================

tfam <- read.table("~/work/ncsu/projects/diallel2/reportData/top.snp.geno.tfam", as.is = TRUE, header = FALSE)
n.ind <- nrow(tfam)
tped <- matrix(scan("~/work/ncsu/projects/diallel2/reportData/top.snp.geno.tped", what = "", na.strings = "0"), ncol = 4 + 2*n.ind, byrow = TRUE)

geno <- tped[, seq(5, ncol(tped), 2)]
dgrp.lines <- tfam[, 1]
colnames(geno) <- dgrp.lines
rownames(geno) <- tped[, 2]

# prod data
prod <- read.csv("~/work/ncsu/projects/diallel2/reportData/dgrp.prod.data.csv", header = TRUE, as.is = TRUE, na.strings = ".")
prod[is.na(prod[, 3]), 3] <- 0

prod25 <- prod[prod[, 4] == 25 & prod[, 5] == 25, ]
dgrp.prod25 <- sapply(split(prod25[, 3], paste("line_", prod25[, 1], sep = "")), mean)

diallel.raw.data <- read.csv("~/work/ncsu/projects/diallel2/diallel.csv", header = TRUE, as.is = TRUE, na.strings = "")
diallel.self.cross <- with(subset(diallel.raw.data, female.parent == male.parent), sapply(split(n.total, female.parent), mean))
names(diallel.self.cross) <- paste("line_", names(diallel.self.cross), sep = "")

# common line
common.line <- intersect(names(dgrp.prod25), names(diallel.self.cross))
cor.test(dgrp.prod25[common.line], diallel.self.cross[common.line])

# dgrp.geno
dgrp.geno <- geno[, names(dgrp.prod25)]

for (i in 1:nrow(dgrp.geno)) {
  
  cat(rownames(dgrp.geno)[i], anova(lm(dgrp.prod25 ~ dgrp.geno[i, ]))$Pr[1], "\n")
  
}

# dop2r expression from Everett et al 2020
dop2r <- read.table(file = "~/work/ncsu/projects/diallel2/reportData/dop2r.exp.everett2020.csv", header = T, as.is = T)

rownames(dop2r) <- paste("line_", gsub("_F", "", dop2r[, 1]), sep = "")

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "report/figureDGRP.pdf")


library("RColorBrewer")
library("mediation")
library("lsmeans")

# prepare file
# ============================================================

file.width = 89 # in mm
cairo_pdf(file = args[2], width = file.width/25.4, height = file.width/25.4*0.5, family = args[1])
par(mfrow = c(1, 2), las = 1, tcl = -0.2, mar = c(1.8, 1.8, 0.8, 0.2), ps = 7, lwd = 0.5, xpd = TRUE)

# plot DGRP productivity versus genotype
prod.fit <- aov(y ~ x, data = data.frame(y = dgrp.prod25, x = factor(dgrp.geno[15, ])))
prod.lsmeans <- summary(lsmeans(prod.fit, ~ x))

plot(c(0.5, 2.5), c(0, 210), type = "n", xlab = "", axes = FALSE, ylab = "")
points(runif(sum(dgrp.geno[15, ] == "G", na.rm = T), -0.1, 0.1) + 1, dgrp.prod25[which(dgrp.geno[15, ] == "G")], cex = 0.3)
points(runif(sum(dgrp.geno[15, ] == "T", na.rm = T), -0.1, 0.1) + 2, dgrp.prod25[which(dgrp.geno[15, ] == "T")], cex = 0.3)


segments(x0 = c(1, 2), y0 = prod.lsmeans$lower.CL,
         x1 = c(1, 2), y1 = prod.lsmeans$upper.CL, col = c("red", "blue"), lwd = 0.5)

segments(x0 = 1:2 - 0.1, y0 = prod.lsmeans$lsmean,
         x1 = 1:2 + 0.1, y1 = prod.lsmeans$lsmean, col = c("red", "blue"), lwd = 2)

axis(side = 1, at = 1:2, labels = c("G", "T"), lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 6/par("ps")/par("cex"))
axis(side = 2, at = seq(0, 200, 50), lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 6/par("ps")/par("cex"))
title(xlab = expression(paste(italic("Dop2R"), " SNP allele")), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(ylab = "Number of adult offspring", cex.lab = 7/par("ps")/par("cex"), mgp = c(1, 0, 0))
box(bty = "l")
text(1.5, 200, expression(paste(italic("P"), " = 0.0069")), cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05 , from = "inches", to = "user") , grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1)

# expression association with productivity
# ============================================================

common.line <- intersect(rownames(dop2r), colnames(dgrp.geno))
common.exp <- dop2r[common.line, 2]
common.geno <- dgrp.geno[15, common.line]
common.prod <- dgrp.prod25[common.line]

plot(c(1.1, 2.7), c(20, 210), type = "n", xlab = "", axes = FALSE, ylab = "")
points(common.exp, common.prod, cex = 0.5)
cor.test(common.exp, common.prod)

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 6/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 6/par("ps")/par("cex"))
title(xlab = expression(paste(italic("Dop2R"), " expression (log", {}[2], "FPKM)")), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(ylab = "Number of adult offspring", cex.lab = 7/par("ps")/par("cex"), mgp = c(1, 0, 0))
box(bty = "l")
text(1.8, 200, expression(paste(italic("r"), " = 0.34 (", italic("P"), " = 9.81e-7)")), cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user") , grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1)


# mediation analysis
# ============================================================

this.data <- na.omit(data.frame(prod = common.prod, geno = common.geno, exp = common.exp))

m <- lm(exp ~ geno, data = this.data)
y <- lm(prod ~ exp + geno, data = this.data)
res <- mediate(m, y, sims = 10000, treat = "geno", mediator = "exp")
summary(res)

dev.off()
