# plot histogram after wolbachia adjustment
# ============================================================

args <- commandArgs(TRUE) # args <- c("diallel.csv", "wolba.csv", -5.91, -6.01, -3.65, "report/figureHistogram.pdf")
# the adjustments applied to yy yn ny
library("RColorBrewer")

diallel.file <- args[2]
wolba.file <- args[3]
yy <- as.numeric(args[4])
yn <- as.numeric(args[5])
ny <- as.numeric(args[6])
output.file <- args[7]

# read raw data
# ============================================================

diallel <- read.csv(diallel.file, header = TRUE, as.is = TRUE)
wolba <- read.csv(wolba.file, header = TRUE, as.is = TRUE)

diallel <- within(diallel, {female.parent <- paste("line_", female.parent, sep = "");
                            male.parent <- paste("line_", male.parent, sep = "")})


rownames(wolba) <- wolba[, 1]
diallel$female.parent.wolba <- wolba[diallel$female.parent, 2]
diallel$male.parent.wolba <- wolba[diallel$male.parent, 2]
diallel$wolbaint <- paste(diallel$female.parent.wolba, diallel$male.parent.wolba, sep = "")
diallel <- within(diallel, wolbaint <- factor(wolbaint, levels = c("nn", "ny", "yn", "yy")))

diallel$adjusted <- NA
diallel <- within(diallel, { adjusted[wolbaint == "nn"] <- n.total[wolbaint == "nn"]
                             adjusted[wolbaint == "yy"] <- n.total[wolbaint == "yy"] + yy;
                             adjusted[wolbaint == "yn"] <- n.total[wolbaint == "yn"] + yn; 
                             adjusted[wolbaint == "ny"] <- n.total[wolbaint == "ny"] + ny;
                            })

# prepare file
# ============================================================

file.width = 45 # in mm
cairo_pdf(file = output.file, width = file.width/25.4, height = file.width/25.4, family = args[1])
par(las = 1, tcl = -0.2, mar = c(1.5, 1.5, 0.8, 0.2), ps = 7, lwd = 0.5, xpd = TRUE)

hist(diallel$adjusted, breaks = seq(0, 230, 10), col = "grey90", main = "", xlab = "Productivity", ylab = "Frequency", freq = TRUE, axes = FALSE)
axis(side = 2, mgp = c(2, 0.3, 0), cex.axis = 6/par("ps")/par("cex"), lwd = 0.5)
axis(side = 1, at = seq(0, 200, 50), mgp = c(2, -0.1, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
box(bty = "l")
title(xlab = "Productivity", mgp = c(0.5, 0, 0), cex.lab = 6/par("ps")/par("cex"))
title(ylab = "Count", mgp = c(0.9, 0, 0), cex.lab = 6/par("ps")/par("cex"))

# find the area under curve for the normal curve
# ============================================================

x.seq <- seq(0, 240, 1)
y.seq <- dnorm(x.seq, mean = mean(diallel$adjusted), sd = sd(diallel$adjusted))
lines(x.seq, y.seq*length(diallel$adjusted)*10)
text(110, 800, "Normal density", pos = 4)
text(130, 600, paste("CV = ", formatC(sd(diallel$adjusted)/mean(diallel$adjusted)*100, format = "f", digits = 2), "%", sep = ""), pos = 4)
