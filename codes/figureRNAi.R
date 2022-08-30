# =========================
# = analyze the RNAi data =
# =========================

# there is one row that includes extra flies on day 20, (rep3 of CG7510 for male 375)
# those were added to the day 16 data to keep the columns consistent 

# also delete data on 5/29 for V60000, those were off
# ============================================================

args <- commandArgs(TRUE) # args <- c("Myriad Pro", "reportData/prodRNAi.csv", "report/figureRNAi.pdf")
library("RColorBrewer")
library("lsmeans")
library("multcomp")
library("agricolae")

# prepare file
# ============================================================

file.width = 89 # in mm
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width/25.4, family = args[1])
par(mfrow = c(2, 1), las = 1, tcl = -0.2, mar = c(2, 2, 0.8, 0.2), ps = 7, lwd = 0.5, xpd = TRUE)

# read data
# ============================================================

rnai.raw.data <- read.csv(args[2], header = TRUE, as.is = TRUE, na.strings = "")
rnai.raw.data[is.na(rnai.raw.data)] <- as.integer(0)

# remove data when there are female or male dead
rnai.prod.data <- subset(rnai.raw.data, (d1_female_dead_at_transfer == 0 & d1_male_dead_at_transfer == 0 &
                                         d2_female_dead_at_transfer == 0 & d2_male_dead_at_transfer == 0))

# remove gene U26, this was tested but not used for the publication
# reason being that it's not within any 2kb of SNPs
rnai.prod.data <- subset(rnai.prod.data, female != "U26")
rnai.prod.data <- within(rnai.prod.data, { female[female == "CG7510"] <- "anchor"; female[female == "D2R"] <- "Dop2R" })

# calculate total number of eggs and adults
rnai.prod.data <- within(rnai.prod.data, { total_eggs <- d1_eggs + d2_eggs;
                                           total_adults <- d1_d10female + d1_d10male + d1_d13female + d1_d13male +
                                           d1_d16female + d1_d16male + d1_d18female + d1_d18male +
                                           d2_d10female + d2_d10male + d2_d13female + d2_d13male +
                                           d2_d16female + d2_d16male;
                                           female <- factor(female, levels = c("V60000", "Dop2R", "CG17003", "Or49a", "CG30048", "CG31897", "anchor"));
                                           male <- factor(male); })
# remove those data points with more adults than eggs
rnai.prod.data <- subset(rnai.prod.data, (total_adults <= total_eggs))
                                           

# calculate hatching rate
rnai.prod.data <- within(rnai.prod.data, hatch <- total_adults/total_eggs)
                                         
# fit data
# ============================================================

adult.fit <- aov(total_adults ~ male + female, data = rnai.prod.data)
egg.fit <- aov(total_eggs ~ female, data = rnai.prod.data)
hatch.fit <- glm(cbind(total_adults, total_eggs - total_adults) ~ male + female, data = rnai.prod.data, family = "binomial")

egg.lsmeans <- summary(lsmeans(egg.fit, ~ female))
egg.hsd <- HSD.test(egg.fit, "female", group = T, console = T)
adult.lsmeans <- summary(lsmeans(adult.fit, ~ female))
adult.hsd <- HSD.test(adult.fit, "female", group = T, console = T)


plot(c(1, 7), c(0, 210), type = "n", xlab = "", axes = FALSE, ylab = "")
with(rnai.prod.data, points(as.numeric(female) + runif(nrow(rnai.prod.data), -0.1, 0.1), total_adults, cex = 0.5))

segments(x0 = 1:7, y0 = adult.lsmeans$lower.CL,
         x1 = 1:7, y1 = adult.lsmeans$upper.CL, col = "grey20", lwd = 0.5)

segments(x0 = 1:7 - 0.1, y0 = adult.lsmeans$lsmean,
         x1 = 1:7 + 0.1, y1 = adult.lsmeans$lsmean, col = "grey20", lwd = 2)

axis(side = 1, at = 1:7, labels = c("Control", parse(text = paste("italic(\"", levels(rnai.prod.data$female)[-1], "\")", sep = ""))), lwd = 0.5, mgp = c(0.8, 0.1, 0), cex.axis = 6/par("ps")/par("cex"))
axis(side = 2, at = seq(0, 200, 50), lwd = 0.5, mgp = c(0.8, 0.4, 0), cex.axis = 6/par("ps")/par("cex"))
title(xlab = "Genes knocked down in female parents", cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(ylab = "Number of adults", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))
box(bty = "l")
text(1:7, 25, labels = adult.hsd$groups$groups[match(levels(rnai.prod.data$female), gsub("[[:space:]]", "", rownames(adult.hsd$groups)))], 6/par("ps")/par("cex"), pos = 1)
text(grconvertX(0.05 , from = "inches", to = "user") , grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 8/par("ps")/par("cex"), pos = 1)



plot(c(1, 7), c(0, 300), type = "n", xlab = "", axes = FALSE, ylab = "")
with(rnai.prod.data, points(as.numeric(female) + runif(nrow(rnai.prod.data), -0.1, 0.1), total_eggs, cex = 0.5))

segments(x0 = 1:7, y0 = egg.lsmeans$lower.CL,
         x1 = 1:7, y1 = egg.lsmeans$upper.CL, col = "grey20", lwd = 0.5)

segments(x0 = 1:7 - 0.1, y0 = egg.lsmeans$lsmean,
         x1 = 1:7 + 0.1, y1 = egg.lsmeans$lsmean, col = "grey20", lwd = 2)

axis(side = 1, at = 1:7, labels = c("Control", parse(text = paste("italic(\"", levels(rnai.prod.data$female)[-1], "\")", sep = ""))), lwd = 0.5, mgp = c(0.8, 0, 0), cex.axis = 6/par("ps")/par("cex"))
axis(side = 2, at = seq(0, 300, 50), lwd = 0.5, mgp = c(0.8, 0.4, 0), cex.axis = 6/par("ps")/par("cex"))
title(xlab = "Genes knocked down in female parents", cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(ylab = "Number of eggs", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))
box(bty = "l")
text(1:7, 35, labels = egg.hsd$groups$groups[match(levels(rnai.prod.data$female), gsub("[[:space:]]", "", rownames(egg.hsd$groups)))], 6/par("ps")/par("cex"), pos = 1)
text(grconvertX(0.05 , from = "inches", to = "user") , grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 8/par("ps")/par("cex"), pos = 1)


dev.off()
