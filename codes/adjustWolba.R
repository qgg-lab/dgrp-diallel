# adjust Wolbachia effect
# ============================================================

args <- commandArgs(TRUE) # args <- c("diallel.csv", "wolba.csv", 5.91, 6.01, 3.65, "diallelWolbaAdjusted.csv")
# the adjustments applied to yy yn ny

diallel.file <- args[1]
wolba.file <- args[2]
yy <- as.numeric(args[3])
yn <- as.numeric(args[4])
ny <- as.numeric(args[5])
output.file <- args[6]

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

# output data
# ============================================================

write.table(diallel[, c("female.parent", "male.parent", "n.total", "adjusted")], file = output.file, col.names = F, row.names = F, quote = F, sep = ",")
