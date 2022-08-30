# ===========================
# = test for sex ratio bias =
# ===========================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "diallel.csv", "report/figureSexRatio.pdf", "reportData/SexRatioChisq.RData")
library("RColorBrewer")

# prepare file
# ============================================================

file.width = 89 # in mm
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width/25.4*0.5, family = args[1])
par(las = 1, tcl = -0.2, mfrow = c(1, 2), mar = c(1.5, 1.5, 1, 1), ps = 7, lwd = 0.5, xpd = TRUE)

# read and process data
# ============================================================

diallel <- read.csv(args[2], header = TRUE, as.is = TRUE)
diallel <- within(diallel, { female.parent <- factor(paste("DGRP_", sprintf("%03d", female.parent), sep = ""));
                             male.parent   <- factor(paste("DGRP_", sprintf("%03d", male.parent), sep = "")); })

# make plot
# ============================================================

with(diallel, plot(n.female, n.male, pch = 3, axes = FALSE, xlim = c(0, 130), xlab = "", ylab = "", cex = 0.2, lwd = 0.2, asp = 1))
abline(a = 0, b = 1, xpd = FALSE, lwd = 0.5, lty = 2)
axis(side = 1, at = seq(0, 120, 40), mgp = c(2, -0.2, 0), cex.axis = 6/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, at = seq(0, 120, 40), mgp = c(2, 0.3, 0), cex.axis = 6/par("ps")/par("cex"), lwd = 0.5)
title(xlab = "Number of female progeny", mgp = c(0.3, 0, 0), cex.lab = 6/par("ps")/par("cex"))
title(ylab = "Number of male progeny", mgp = c(0.8, 0, 0), cex.lab = 6/par("ps")/par("cex"))
box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user") , grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1)


# test for significance

par(mar = c(1.5, 1.8, 1, 1))
chisq.p <- with(diallel, apply(cbind(n.female, n.male), 1, function(x) { chisq.test(x, p = c(0.5, 0.5), simulate.p.value = TRUE, B = 10)$p.value }))
hist(chisq.p, breaks = seq(0, 1, 0.05), axes = FALSE, xlab = "", ylab = "", main = "", col = "grey50")
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(2, -0.2, 0), cex.axis = 6/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), cex.axis = 6/par("ps")/par("cex"), lwd = 0.5)
title(xlab = expression(paste(italic(P), " value")), mgp = c(0.3, 0, 0), cex.lab = 6/par("ps")/par("cex"))
title(ylab = "Frequency", mgp = c(1.2, 0, 0), cex.lab = 6/par("ps")/par("cex"))
box(bty = "l")
save(chisq.p, file = args[4])

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user") , grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1)

dev.off()
