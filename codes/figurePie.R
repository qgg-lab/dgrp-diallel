# plot pie chart for variance component
# ============================================================

args <- commandArgs(TRUE) # args <- c("Myrid Pro", 721.41, 407.75, 27.63, 35.30, 8.25, 6.96, "report/figurePie.pdf")
# the order of the number is residual, f, m, fm, n, nn
library("RColorBrewer")

r <- as.numeric(args[2])
f <- as.numeric(args[3])
m <- as.numeric(args[4])
fm <- as.numeric(args[5])
n <- as.numeric(args[6])
nn <- as.numeric(args[7])

# prepare file
# ============================================================

file.width = 45 # in mm
cairo_pdf(file = args[8], width = file.width/25.4, height = file.width/25.4*0.8, family = args[1])
par(las = 1, tcl = -0.2, mar = c(0, 0, 1, 2), ps = 7, lwd = 0.5, xpd = TRUE)

# make plot
# modify pie() to return label position and add label later
# ============================================================

mypie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("'x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
                "lavender", "cornsilk")
        else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
    p.list <- list()
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
        P <- t2xy(mean(x[i + 0:1]))
        p.list[[i]] <- P
        # lab <- as.character(labels[i])
#         if (!is.na(lab) && nzchar(lab)) {
#             lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
#             text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE,
#                 adj = ifelse(P$x < 0, 1, 0), ...)
#         }
    }
    return(p.list)
    title(main = main, ...)
    invisible(NULL)
}

p.pos <- mypie(c(r, f, m, fm, n, nn), edges = 500, clockwise = TRUE, init.angle = 0 - (m + fm + n + nn)/(r + f + m + fm + n + nn)*360/2,
                col = brewer.pal(9, "Set1")[c(9, 1:5)], border = NA)

# manually add labels
lines(c(1, 1.05)*p.pos[[1]]$x, c(1, 1.05)*p.pos[[1]]$y, lwd = 0.5, col = brewer.pal(9, "Set1")[9])
text(1.2*p.pos[[1]]$x, 1.2*p.pos[[1]]$y, parse(text = paste("paste(italic(e), \": \",", formatC(r/(r + f + m + fm + n +nn)*100, format = "f", digits = 2), ", \"%\")")), cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Set1")[9])

lines(c(1, 1.05)*p.pos[[2]]$x, c(1, 1.05)*p.pos[[2]]$y, lwd = 0.5, col = brewer.pal(9, "Set1")[1])
text(1.2*p.pos[[2]]$x, 1.2*p.pos[[2]]$y, parse(text = paste("paste(italic(f), \": \",", formatC(f/(r + f + m + fm + n +nn)*100, format = "f", digits = 2), ", \"%\")")), cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Set1")[1])

lines(c(1, 1.05)*p.pos[[3]]$x, c(1, 1.05)*p.pos[[3]]$y, lwd = 0.5, col = brewer.pal(9, "Set1")[2])
text(1.1*p.pos[[3]]$x, 1.1*p.pos[[3]]$y, parse(text = paste("paste(italic(m), \": \",", formatC(m/(r + f + m + fm + n +nn)*100, format = "f", digits = 2), ", \"%\")")), cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Set1")[2], adj = 0)

lines(c(1, 1.05)*p.pos[[4]]$x, c(1, 1.05)*p.pos[[4]]$y, lwd = 0.5, col = brewer.pal(9, "Set1")[3])
text(1.1*p.pos[[4]]$x, 1.1*p.pos[[4]]$y, parse(text = paste("paste(italic(fm), \": \",", formatC(fm/(r + f + m + fm + n +nn)*100, format = "f", digits = 2), ", \"%\")")), cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Set1")[3], adj = 0)

lines(c(c(1, 1.05)*p.pos[[5]]$x, 1.05*p.pos[[5]]$x + 0.05), c(c(1, 1.05)*p.pos[[5]]$y, 1.05*p.pos[[5]]$y - 0.02), lwd = 0.5, col = brewer.pal(9, "Set1")[4])
text(1.1*p.pos[[5]]$x + 0.03, 1.1*p.pos[[5]]$y - 0.03, parse(text = paste("paste(italic(n), \": \",", formatC(n/(r + f + m + fm + n +nn)*100, format = "f", digits = 2), ", \"%\")")), cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Set1")[4], adj = 0)

lines(c(c(1, 1.05)*p.pos[[6]]$x, 1.05*p.pos[[6]]$x + 0.05), c(c(1, 1.05)*p.pos[[6]]$y, 1.05*p.pos[[6]]$y - 0.1), lwd = 0.5, col = brewer.pal(9, "Set1")[5])
text(1.1*p.pos[[6]]$x + 0.03, 1.1*p.pos[[6]]$y - 0.12, parse(text = paste("paste(italic(nn), \": \",", formatC(nn/(r + f + m + fm + n +nn)*100, format = "f", digits = 2), ", \"%\")")), cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Set1")[5], adj = 0)


text(0.4, 1.2, expression(italic(y[ijk]) * " = " * italic(mu) * " + " * phantom(italic(f[i])) * " + " * phantom(italic(m[j])) * " + " * phantom(italic(fm[ij])) * " + " *
                          phantom(italic(n[i])) * " + " * phantom(italic(n[j])) * " + " * phantom(italic(nn[ij])) * " + " * phantom(italic(e[ijk]))))
text(0.4, 1.2, expression(phantom(italic(y[ijk])) * phantom(" = ") * phantom(italic(mu)) * phantom(" + ") * italic(f[i]) * phantom(" + ") * phantom(italic(m[j])) * phantom(" + ") * phantom(italic(fm[ij])) *
                           phantom(" + ") * phantom(italic(n[i])) * phantom(" + ") * phantom(italic(n[j])) * phantom(" + ") * phantom(italic(nn[ij])) * phantom(" + ") *
                           phantom(italic(e[ijk]))), col = brewer.pal(9, "Set1")[1])
text(0.4, 1.2, expression(phantom(italic(y[ijk])) * phantom(" = ") * phantom(italic(mu)) * phantom(" + ") * phantom(italic(f[i])) * phantom(" + ") * italic(m[j]) * phantom(" + ") * phantom(italic(fm[ij])) *
                          phantom(" + ") * phantom(italic(n[i])) * phantom(" + ") * phantom(italic(n[j])) * phantom(" + ") * phantom(italic(nn[ij])) * phantom(" + ") *
                          phantom(italic(e[ijk]))), col = brewer.pal(9, "Set1")[2])
text(0.4, 1.2, expression(phantom(italic(y[ijk])) * phantom(" = ") * phantom(italic(mu)) * phantom(" + ") * phantom(italic(f[i])) * phantom(" + ") * phantom(italic(m[j])) * phantom(" + ") * italic(fm[ij]) *
                          phantom(" + ") * phantom(italic(n[i])) * phantom(" + ") * phantom(italic(n[j])) * phantom(" + ") * phantom(italic(nn[ij])) * phantom(" + ") *
                          phantom(italic(e[ijk]))), col = brewer.pal(9, "Set1")[3])
text(0.4, 1.2, expression(phantom(italic(y[ijk])) * phantom(" = ") * phantom(italic(mu)) * phantom(" + ") * phantom(italic(f[i])) * phantom(" + ") * phantom(italic(m[j])) * phantom(" + ") *
                          phantom(italic(fm[ij])) * phantom(" + ") * italic(n[i]) * phantom(" + ") * italic(n[j]) * phantom(" + ") *
                          phantom(italic(nn[ij])) * phantom(" + ") * phantom(italic(e[ijk]))), col = brewer.pal(9, "Set1")[4])
text(0.4, 1.2, expression(phantom(italic(y[ijk])) * phantom(" = ") * phantom(italic(mu)) * phantom(" + ") * phantom(italic(f[i])) * phantom(" + ") * phantom(italic(m[j])) * phantom(" + ") *
                          phantom(italic(fm[ij])) * phantom(" + ") * phantom(italic(n[i])) * phantom(" + ") * phantom(italic(n[j])) * phantom(" + ") *
                          italic(nn[ij]) * phantom(" + ") * phantom(italic(e[ijk]))), col = brewer.pal(9, "Set1")[5])
text(0.4, 1.2, expression(phantom(italic(y[ijk])) * phantom(" = ") * phantom(italic(mu)) * phantom(" + ") * phantom(italic(f[i])) * phantom(" + ") * phantom(italic(m[j])) * phantom(" + ") *
                          phantom(italic(fm[ij])) * phantom(" + ") * phantom(italic(n[i])) * phantom(" + ") * phantom(italic(n[j])) * phantom(" + ") *
                          phantom(italic(nn[ij])) * phantom(" + ") * italic(e[ijk])), col = brewer.pal(9, "Set1")[9])

dev.off()
