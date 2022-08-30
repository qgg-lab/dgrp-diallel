# ===========================================
# = variance component for the diallel data =
# ===========================================

library("regress")

# read raw data
# ============================================================

diallel <- read.csv("diallel.csv", header = TRUE, as.is = TRUE)
# diallel <- subset(diallel, female.parent != 398 & male.parent != 398)
diallel <- within(diallel, { female.parent <- factor(female.parent); male.parent <- factor(male.parent); })

# simple mixed model
# ============================================================

fm.fit <- regress(n.total ~ 1, ~ female.parent + male.parent + I(female.parent:male.parent), identity = TRUE, data = diallel)

# complex bio model
# ============================================================

z.f <- with(diallel, model.matrix(~ female.parent - 1)); g.f <- z.f %*% t(z.f)
z.m <- with(diallel, model.matrix(~ male.parent - 1)); g.m <- z.m %*% t(z.m)
z.fm <- with(diallel, model.matrix(~ female.parent:male.parent - 1)); g.fm <- z.fm %*% t(z.fm)
z.n <- with(diallel, model.matrix(~ female.parent - 1) + model.matrix(~ male.parent - 1)); g.n <- z.n %*% t(z.n); 
z.nn <- with(diallel, model.matrix(~ factor(apply(cbind(female.parent, male.parent), 1, function(x){ paste(sort(x), collapse = "_") })) - 1 )); g.nn <- z.nn %*% t(z.nn)

g.f <- g.f/mean(diag(g.f))
g.m <- g.m/mean(diag(g.m))
g.fm <- g.fm/mean(diag(g.fm))
g.n <- g.n/mean(diag(g.n))
g.nn <- g.nn/mean(diag(g.nn))

bio.fit <- regress(n.total ~ 1, ~ g.f + g.m + g.fm + g.n + g.nn, identity = TRUE, data = diallel)

save(fm.fit, bio.fit, file = "qg/varComp.RData")
