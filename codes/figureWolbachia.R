# ===============================
# = figure for wolbachia effect =
# ===============================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "diallel.csv", "wolba.csv", "report/figureWolbachia.pdf")
library("RColorBrewer")

# prepare file
# ============================================================

file.width = 45 # in mm
cairo_pdf(file = args[4], width = file.width/25.4, height = file.width/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(1.5, 1.5, 0.8, 0.2), ps = 7, lwd = 0.5, xpd = TRUE)

# read and process data
# ============================================================

diallel <- read.csv(args[2], header = TRUE, as.is = TRUE)
wolba <- read.csv(args[3], header = TRUE, as.is = TRUE)

diallel <- within(diallel, {female.parent <- paste("line_", female.parent, sep = "");
                            male.parent <- paste("line_", male.parent, sep = "")})

rownames(wolba) <- wolba[, 1]
diallel$female.parent.wolba <- wolba[diallel$female.parent, 2]
diallel$male.parent.wolba <- wolba[diallel$male.parent, 2]
diallel$wolbaint <- paste(diallel$female.parent.wolba, diallel$male.parent.wolba, sep = "")
diallel <- within(diallel, wolbaint <- factor(wolbaint, levels = c("nn", "ny", "yn", "yy")))


# make plot
# ============================================================

with(diallel, boxplot(n.total ~ wolbaint, at = c(1, 2, 4, 5), axes = FALSE, xlab = "", ylab = "", ylim = c(0, 275), range = 0, col = c("white", "grey50", "white", "grey50")))
axis(side = 1, at = c(1.5, 4.5), labels = c("-", "+"), mgp = c(2, -0.1, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), cex.axis = 6/par("ps")/par("cex"), lwd = 0.5)
title(xlab = expression(paste(italic("Wolbachia"), " infection of female parent")), mgp = c(0.5, 0, 0), cex.lab = 6/par("ps")/par("cex"))
title(ylab = "Number of progeny", mgp = c(0.9, 0, 0), cex.lab = 6/par("ps")/par("cex"))
legend(0, 310, pch = 22, pt.bg = c("white", "grey50"), legend = c("-", "+"), bty = "n", x.intersp = 0.5, horiz = TRUE, y.intersp = 0.7)
text(0, 290, expression(paste(italic("Wolbachia"), " infection of male parent")), pos = 4, cex = 5/par("ps")/par("cex"))
box(bty = "l")

dev.off()
