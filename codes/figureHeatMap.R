# ===================================
# = make heat map for diallel cross =
# ===================================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "diallel.csv", "report/figureHeatMap.pdf")
library("RColorBrewer")

# prepare file
# ============================================================

file.width = 89 # in mm
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width/25.4*0.92, family = args[1])
par(las = 1, tcl = -0.2, mar = c(0, 2.2, 1.5, 1), ps = 7, lwd = 0.5, xpd = TRUE)

# read and process data
# ============================================================

diallel <- read.csv(args[2], header = TRUE, as.is = TRUE)
diallel <- within(diallel, { female.parent <- factor(paste("DGRP_", sprintf("%03d", female.parent), sep = ""));
                             male.parent   <- factor(paste("DGRP_", sprintf("%03d", male.parent), sep = "")); })

cross.mean <- tapply(diallel$n.total, diallel[, c("male.parent", "female.parent")], mean)
line.order <- levels(diallel$female.parent)
cross.mean <- cross.mean[line.order, rev(line.order)]

# make plot
# ============================================================

image(1:50, 1:50, cross.mean, axes = F, xlab = "", ylab = "", xlim = c(0, 55), breaks = seq(0, 200, 4), col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(50), asp = 1)
segments(0.5, seq(0.5, 50.5, 1), 50.5, seq(0.5, 50.5, 1), lwd = 0.2)
segments(seq(0.5, 50.5, 1), 0.5, seq(0.5, 50.5, 1), 50.5, lwd = 0.2)

axis(side = 3, at = 1:50, labels = rownames(cross.mean), las = 2, cex.axis = 4/par("ps")/par("cex"), mgp = c(0, 0.2, 0), tck = -0.01, pos = 50.5, lwd = 0.5)
axis(side = 2, at = 1:50, labels = colnames(cross.mean), pos = 0.5, cex.axis = 4/par("ps")/par("cex"), mgp = c(0, 0.2, 0), tck = -0.01, lwd = 0.5)
mtext("Male parent", side = 3, at = 25, line = 0.8)
mtext("Female parent", side = 2, at = 25, line = 1.5, las = 3)
rect(0.5, 0.5, 50.5, 50.5, lwd = 0.5)

image(c(52, 53), 1:50, matrix(c(seq(0, 196, 4) + 1, seq(0, 196, 4) + 1), nrow = 2, byrow = TRUE), breaks = seq(0, 200, 4), col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(50), add = T)
rect(51.5, 0.5, 53.5, 50.5, lwd = 0.5)
axis(side = 4, at = seq(0.5, 50.5, 5), labels = seq(0, 200, 20), cex.axis = 5/par("ps")/par("cex"), las = 2, pos = 53.5, mgp = c(0, 0.3, 0), lwd = 0.5)
dev.off()
