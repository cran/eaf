###############################################################################
#
#                          Copyright (c) 2011
#         Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
#             Marco Chiarandini <marco@imada.sdu.dk>
#
# This program is free software (software libre); you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, you can obtain a copy of the GNU
# General Public License at http://www.gnu.org/copyleft/gpl.html
#
# IMPORTANT NOTE: Please be aware that the fact that this program is
# released as Free Software does not excuse you from scientific
# propriety, which obligates you to give appropriate credit! If you
# write a scientific paper describing research that made substantive use
# of this program, it is your obligation as a scientist to (a) mention
# the fashion in which this software was used in the Methods section;
# (b) mention the algorithm in the References section. The appropriate
# citation is:
# 
#    Manuel Lopez-Ibanez, Luis Paquete, and Thomas Stuetzle.
#    Exploratory Analysis of Stochastic Local Search Algorithms in
#    Biobjective Optimization. In T. Bartz-Beielstein, M. Chiarandini,
#    L. Paquete, and M. Preuss, editors, Experimental Methods for the
#    Analysis of Optimization Algorithms, pages 209-222. Springer,
#    Berlin, Germany, 2010. doi: 10.1007/978-3-642-02538-9_9
# 
# Moreover, as a personal note, I would appreciate it if you would email
# manuel.lopez-ibanez@ulb.ac.be with citations of papers referencing this
# work so I can mention them to my funding agent and tenure committee.
#
################################################################################
#
# TODO:
#
#  * Follow this style for coding:
#    http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html
#
################################################################################
#dyn.load("../src/eaf.so")

compute.eaf <- function(data, percentiles = NULL)
{
  if (!is.numeric(data[,1]) || !is.numeric(data[,2]) )
    stop("The first and the second column must contain numeric data for point coordinates.\n")
  # The C code expect points within a set to be contiguous.
  data <- data[order(data[,3]),]
  npoints <- aggregate(data[,1],list(data[,3]),length)$x
  nsets <- length(unique(data[,3]))
  if (is.null(percentiles)) {
    percentiles <- 1:nsets * 100 / nsets
  }
  out <- .Call("compute_eaf_C",
               as.double(t(as.matrix(data[,1:2]))),
               2,
               as.integer(cumsum(npoints)),
               as.integer(nsets),
               as.integer(percentiles)
               )

  return(out)
}

compute.eaf.as.list <- function(data, percentiles = NULL)
{
  eaf <- compute.eaf (data, percentiles = percentiles)
  return (lapply(split.data.frame(eaf, as.factor(eaf[,3])),
                 function(x) { x[,-3, drop = FALSE] }))
}

compute.eafdiff <- function(DATA, intervals = 1)
{
  # The C code expects points within a set to be contiguous.
  DATA <- DATA[order(DATA[,3]),]
  # aggregate is quite slow on large data sets because of paste()
  npoints <- aggregate(DATA[,1], list(DATA[,3]),length)$x
  nsets <- length(unique(DATA[,3]))
  division <- nsets %/% 2
  nsets1 <- division
  nsets2 <- nsets - division
  DIFF <- .Call("compute_eafdiff_C",
                as.double(t(as.matrix(DATA[,1:2]))),
                2,
                as.integer(cumsum(npoints)),
                as.integer(nsets),
                as.integer(division),
                as.integer(intervals)
                )
  # FIXME: Do this computation in C code.
  eafdiff <- list(left=NULL, right=NULL)
  eafdiff$left <- unique(DIFF[ DIFF[,3] >= 1, , drop = FALSE])
  eafdiff$right <- unique(DIFF[ DIFF[,3] <= -1, , drop = FALSE])
  eafdiff$right[,3] <- -eafdiff$right[,3]
  return(eafdiff)
}

read.data.sets <- function(file, col.names)
{
  stopifnot(file.exists(file))
  out <- .Call("read_data_sets",
               as.character(file))
  if (missing(col.names))
    col.names <- paste("V", 1L:(ncol(out)-1), sep = "")
  colnames(out) <- c(col.names, "set")
  return(as.data.frame(out))
}

## Calculate the intermediate points in order to plot a staircase-like
## polygon.
## Example: given ((1,2), (2,1)), it returns ((1,2), (2,2), (2,1)).
points.steps <- function(x, decreasing = FALSE)
{
  x <- rbind(x, cbind(x[-1,1],x[-nrow(x),2]))
  if (decreasing) {
    x <- x[order(-x[,1],x[,2]),]
  } else {
    x <- x[order(x[,1],-x[,2]),]
  }
}

sciNotation <- function(x, digits = 1) {
  if (length(x) > 1) {
    return(append(sciNotation(x[1]), sciNotation(x[-1])))
  }
  if (!x) return(0)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
  as.expression(substitute(base %*% 10^exponent,
      list(base = base, exponent = exponent)))
}

eafs <- function (points = NULL, sets = NULL, groups = NULL, percentiles = NULL)
{
  # FIXME: check conditions, ncol(points) == 2, etc.
  # MARCO: done in eafplot.default, it should be enough,
  # unless we want to export also eafs, that is, make it callable directly from users
  # exportable functions are to be declared in NAMESPACE
  points <- cbind(points, sets)

  if (is.null(groups)) {
    attsurfs <- compute.eaf (points, percentiles)
  } else {
    attsurfs <- data.frame()
    groups <- factor(groups)
    for (g in levels(groups)) {
      tmp <- compute.eaf(points[groups == g,], percentiles)
      attsurfs <- rbind(attsurfs, data.frame(tmp, groups = g))
    }
  }
  return (attsurfs)
}

eafplot.default <-
  function (x, sets = NULL, groups = NULL,
            percentiles = c(0,50,100),
            attsurfs = NULL,
            xlab = "objective 1", ylab = "objective 2",
            xlim = NULL, ylim = NULL,
            log = "",
            type = "point",
            col = c("black", "darkgrey", "black", "grey40", "darkgrey"),
            lty = c("dashed", "solid", "solid", "solid", "dashed"),
            lwd = c(1.75),
            pch = NA,
            cex.pch = par("cex"),
            las = par("las"),
            legend.pos = "topright",
            legend.txt = NULL,
            # FIXME: Can we get rid of the extra. stuff? Replace it with calling points after eafplot.default in examples and eafplot.pl.
            extra.points = NULL, extra.legend = NULL,
            extra.pch = c(4:25),
            extra.lwd = 0.5,
            extra.lty = "dashed",
            extra.col = "black",
            ymaximise = FALSE,
            xaxis.side = 1, yaxis.side = 2,
            ... )
{
  type <- match.arg (type, c("point", "area"))
  points <- x

  if (ymaximise) points[,2] <- -points[,2]

  if (is.null (xlim)) xlim <- range(points[,1])
  if (is.null (ylim)) ylim <- range(points[,2])

  if (is.null (attsurfs)) {
    points <- as.matrix(points)
    if (!is.matrix(points))
      stop("'points' must be a matrix with exactly two columns")

    if (nrow(points) < 1L || ncol(points) != 2)
      stop("not enough (finite) 'points' observations or not two dimensions")
    
    if (!is.numeric(sets) || length(sets) != nrow(points))
      stop("'sets' must be a vector of length equal to the number of observations in 'points'")

    EAF <- eafs (points, sets, groups, percentiles)

    # Transform EAF matrix into attsurfs list.
    if (is.null(groups)) {
      attsurfs <- lapply(split.data.frame(EAF, as.factor(EAF[,3])),
                         function(x) { x[,-3] })
    } else {
      attsurfs <- list()
      groups <- EAF$groups
      for (g in levels(groups)) {
        tmp <- lapply(split.data.frame(EAF[groups == g,],
                                       as.factor(EAF[groups == g, 3])),
                      function(x) { x[,-3] })
        attsurfs <- c(attsurfs, tmp)
      }
    }
    # FIXME: rm(EAF) to save memory ?
  }

  # FIXME: Make this depend on xlim, ylim.
  extreme <- c(1.0e+20, 1.0e+20)

  # FIXME: Find a better way to handle different x-y scale and
  # maximization/minimization.
  yscale <- 1
  if (ymaximise) {
    #yscale <- 60
    yreverse <- -1
    attsurfs <- lapply (attsurfs, function (x)
                        { x[,2] <- yreverse * x[,2] / yscale; x })
    ylim <- yreverse * ylim / yscale
    extreme[2] <- yreverse * extreme[2]
    if (log == "y") extreme[2] <- 1
  }

  op <- par(cex = 1.0, cex.lab = 1.1, cex.axis = 1.0
            , lab = c(10,5,7)
            )
  ##'lab' A numerical vector of the form 'c(x, y, len)' which modifies
  ##     the default way that axes are annotated.  The values of 'x'
  ##     and 'y' give the (approximate) number of tickmarks on the x
  ##     and y axes and 'len' specifies the label length.  The default
  ##     is 'c(5, 5, 7)'.  Note that this only affects the way the
  ##     parameters 'xaxp' and 'yaxp' are set when the user coordinate
  ##     system is set up, and is not consulted when axes are drawn.
  ##     'len' _is unimplemented_ in R.

  plot(xlim, ylim, type="n", xlab="", ylab="",
       xlim = xlim, ylim = ylim, log = log, axes=FALSE,
       panel.first = ({

         at <- axTicks(1)
         labels <- formatC(at,format="g")
         ## tck=1 draws the vertical grid lines (grid() is seriously broken).
         axis(xaxis.side, at=at, labels=FALSE, tck=1, col='lightgray',
              ## Work-around for R bug:
              lwd=0.5, lty="26") ##  This should be instead: lty='dotted', lwd=par("lwd"))
         axis(xaxis.side, at=at, labels=labels, las = las)
         mtext(xlab, xaxis.side, line=2.1, cex=par("cex.axis"), las=0)

         at <- axTicks(2)
         if (ymaximise) { at <- c(yscale * c(1:24)) }
         labels <- formatC(at,format="g")
         if (log == "y") {
           ## Custom log axis (like gnuplot but in R is hard)
           max.pow <- 6
           at <- c(1, 5, 10, 50, 100, 500, 1000, 1500, 10^c(4:max.pow))
           labels <- c(1, 5, 10, 50, 100, 500, 1000, 1500,
                       parse(text = paste("10^", 4:max.pow, sep = "")))

           #at <- c(60, 120, 180, 240, 300, 480, 600, 900, 1200, 1440)
           #labels <- formatC(at,format="g")

           ## Now do the minor ticks, at 1/10 of each power of 10 interval
           ##at.minor <- 2:9 * rep(c(10^c(1:max.pow)) / 10, each = length(2:9))
           at.minor <- 1:10 * rep(c(10^c(1:max.pow)) / 10, each = length(1:10))
           axis (yaxis.side, at = at.minor, tcl = -0.25, labels = FALSE, las=las)
           axis (yaxis.side, at = at.minor, labels = FALSE, tck=1,
                 col='lightgray', lty='dotted', lwd=par("lwd"))
         }

         ## tck=1 draws the horizontal grid lines (grid() is seriously broken).
         axis (yaxis.side, at=at, labels=FALSE, tck=1,
               col='lightgray', lty='dotted', lwd=par("lwd"))
         axis (yaxis.side, at=at, labels=labels, las = las)
         mtext(ylab, yaxis.side, line=2.75, cex=par("cex.axis"), las=0)

         # FIXME: Perhaps have a function plot.eaf.lines that computes
         # several percentiles for a single algorithm and then calls
         # points() or polygon() as appropriate to add attainment
         # surfaces to an existing plot. This way we can factor out
         # the code below and use it in plot.eaf and plot.eafdiff

         ## Recycle values
         lwd <- rep(lwd, length=length(attsurfs))
         lty <- rep(lty, length=length(attsurfs))
         col <- rep(col, length=length(attsurfs))
         pch <- rep(pch, length=length(attsurfs))

         if (type == "area") {
           for (k in 1:length(attsurfs)) {
             poli <- points.steps(attsurfs[[k]])
             poli <- rbind(c(min(poli[,1]), extreme[2]), poli,
                           c(extreme[1], min(poli[,2])), extreme)
             polygon(poli[,1], poli[,2], border = NA, col = col[k])
           }
         } else {
           for (k in 1:length(attsurfs)) {
             tmp <- attsurfs[[k]][,1:2]
             tmp <- rbind(c(min(tmp[,1]), extreme[2]), tmp[,1:2],
                          c(extreme[1], min(tmp[,2])))

             points(tmp[,1], tmp[,2], type="p", col=col[k], pch = pch[k],
                    cex = cex.pch)
             points(tmp[,1], tmp[,2], type="s", col=col[k], lty = lty[k],
                    lwd = lwd[k])
           }
         }
       }),
       las = las, ...)


  if (!is.null (extra.points)) {
    if (!is.list (extra.points[[1]])) {
      extra.points <- list(extra.points)
    }
    ## Recycle values
    extra.length <- length(extra.points)
    extra.lwd <- rep(extra.lwd, length=extra.length)
    extra.lty <- rep(extra.lty, length=extra.length)
    extra.col <- rep(extra.col, length=extra.length)
    extra.pch <- rep(extra.pch, length=extra.length)

    for (i in 1:length(extra.points)) {
      if (any(is.na(extra.points[[i]][,1]))) {
        ## Extra points are given in the correct order so no reverse
        extra.points[[i]][,2] <- extra.points[[i]][,2] / yscale
        abline(h=extra.points[[i]][,2], lwd = extra.lwd[i], col = extra.col[i],
               lty = extra.lty[i])
        extra.pch[i] <- 0

      } else if (any(is.na(extra.points[[i]][,2]))) {

        abline(v=extra.points[[i]][,1], lwd = extra.lwd[i], col = extra.col[i],
               lty = extra.lty[i])
        extra.pch[i] <- NA

      } else {
        ## Extra points are given in the correct order so no reverse
        extra.points[[i]][,2] <- extra.points[[i]][,2] / yscale
        points (extra.points[[i]], type="p", pch=extra.pch[i],
                col=extra.col[i], cex = cex.pch)

        extra.lty[i] <- "blank"
        extra.lwd[i] <- NA
      }
      lwd <- c(lwd, extra.lwd[i])
      lty <- c(lty, extra.lty[i])
      col <- c(col, extra.col[i])
      pch <- c(pch, extra.pch[i])
    }
  }

  # Setup legend.
  if (is.null(legend.txt) && !is.null(percentiles)) {
    legend.txt <- paste (percentiles, "%", sep="")
    legend.txt <- sub("^0%$", "best", legend.txt)
    legend.txt <- sub("^50%$", "median", legend.txt)
    legend.txt <- sub("^100%$", "worst", legend.txt)

    if (!is.null(groups)) {
      legend.txt <- unlist(tapply(levels(groups), levels(groups), paste, legend.txt))
    }
  }
  legend.txt <- c(legend.txt, extra.legend)

  if (!is.null (legend.txt) && is.na(pmatch(legend.pos,"none"))) {
    if (type == "area") {
      legend(x = legend.pos, y = NULL,
             legend = rev(legend.txt), fill = c(rev(col), "#FFFFFF"),
             bg="white",bty="n", xjust=0, yjust=0, cex=0.9)
    } else {
      legend(legend.pos,
             legend = legend.txt, xjust=1, yjust=1, bty="n",
             lty = lty,  lwd = lwd, pch = pch, col = col, merge=T)
    }
  }

  box()
  par(op)
}


.plot.eafdiff.side <- function (eafdiff, attsurfs = list(),
                               col = c("#FFFFFF", "#BFBFBF","#808080","#404040","#000000"),
                               side = stop("Argument 'side' is required"),
                               type = "point",
                               xlim = NULL, ylim = NULL, log = "",
                               las = par("las"),
                               full.eaf = FALSE,
                               title = "",
                               maximise = c(FALSE, FALSE),
                               xlab = "objective 1", ylab = "objective 2",
                               ...)
{
  side <- match.arg (side, c("left", "right"))
  type <- match.arg (type, c("point", "area"))
  xaxis.side <- if (side == "left") 1 else 3
  yaxis.side <- if (side == "left") 2 else 4
  maximise <- as.logical(maximise)
  
  # Why flooring and not ceiling? If a point has value 2.05, it should
  # be painted with color 2 rather than 3.
  # +1 because col[1] is white.
  eafdiff[,3] <- floor(eafdiff[,3]) + 1
  if (length(unique(eafdiff[,3])) > length(col)) {
    stop ("Too few colors: length(unique(eafdiff[,3])) > length(col)")
  }

  extreme <- c(1.0e+25, 1.0e+25)
  if (!is.null(xlim)) extreme[1] <- ifelse(maximise[1],
                                           xlim[1] - 0.05 * diff(xlim),
                                           xlim[2] + 0.05 * diff(xlim))
  if (!is.null(ylim)) extreme[2] <- ifelse(maximise[2],
                                           ylim[1] - 0.05 * diff(ylim),
                                           ylim[2] + 0.05 * diff(ylim))
  yscale <- 1
#    yscale <- 60
  if (yscale != 1) {
    eafdiff[,2] <- eafdiff[,2] / yscale
    attsurfs <- lapply (attsurfs, function (x)
                        { x[,2] <- x[,2] / yscale; x })
    ylim <- ylim / yscale
    if (log == "y") extreme[2] <- 1
  }

  attsurfs <- lapply (attsurfs,
                      function (x) { x <- rbind(c(min(x[,1]), extreme[2]),
                                                x[,1:2], c(extreme[1], min(x[,2]))) })

  plot(xlim, ylim, type="n", xlab = "", ylab = "",
       ylim = ylim, xlim = xlim, log = log, axes=FALSE,
       panel.first = ({

         at <- axTicks(1)
         labels <- formatC(at,format="g")
         ## tck=1 draws the vertical grid lines (grid() is seriously broken).
         axis(xaxis.side, at=at, labels=FALSE, tck=1, col='lightgray',
              ## Work-around for R bug:
              ##  This should be instead: lty='dotted', lwd=par("lwd"))
              lwd=0.5, lty="26")
         axis(xaxis.side, at=at, labels=labels, las = las)
         mtext(xlab, xaxis.side, line=2.1, las = 0,
               cex = par("cex") * par("cex.axis"))

         at <- axTicks(2)
         #if (ymaximise) { at <- c(yscale * c(1:24)) }
         labels <- formatC(at,format="g")
         if (log == "y") {
           ## Custom log axis (like gnuplot but in R is hard)
           max.pow <- 6
           ##at <- c(1, 5, 10, 50, 100, 500, 1000, 1500, 10^c(4:max.pow))
           ##labels <- c(1, 5, 10, 50, 100, 500, 1000, 1500,
           ##            parse(text = paste("10^", 4:max.pow, sep = "")))

           at <- c(60, 120, 180, 240, 300, 480, 600, 900, 1200, 1440)
           labels <- formatC(at,format="g")

           ## Now do the minor ticks, at 1/10 of each power of 10 interval
           at.minor <- 2:9 * rep(c(10^c(1:max.pow)) / 10, each = length(2:9))
           at.minor <- 1:10 * rep(c(10^c(1:max.pow)) / 10, each = length(1:10))
           ##axis (yaxis.side, at = at.minor, tcl = -0.25, labels = FALSE, las=las)
           ##axis (yaxis.side, at = at.minor, labels = FALSE, tck=1,
           ##      col='lightgray', lty='dotted', lwd=par("lwd"))
         }

         ## tck=1 draws the horizontal grid lines (grid() is seriously broken).
         axis(yaxis.side, at=at, labels=FALSE, tck=1,
              col='lightgray', lty='dotted', lwd=par("lwd"))
         axis(yaxis.side, at=at, labels=labels, las = las)
         mtext(ylab, yaxis.side, line = 2.2, las = 0,
               cex = par("cex") * par("cex.axis"))

         if (nrow(eafdiff)) {
           if (type == "area") {
             for (i in 1:length(col)) {
               if (full.eaf) {
                 poli <- points.steps(eafdiff[eafdiff[,3] == i, c(1,2), drop = FALSE])
                 poli <- rbind(c(min(poli[,1]), extreme[2]), poli,
                               c(extreme[1], min(poli[,2])), c(extreme[1], extreme[2]))
                 polygon(poli[,1], poli[,2], border = NA, col = col[i])
               } else {
                 stop("The combination of parameters: type=\"area\", full.eaf=FALSE is not implemented yet.")
               }
             }
           } else {
             ## The maximum value should also be painted.
             eafdiff[eafdiff[,3] > length(col), 3] <- length(col)
             eafdiff <- eafdiff[order(eafdiff[,3], decreasing = FALSE), , drop=FALSE]
             points(eafdiff[,1], eafdiff[,2], col = col[eafdiff[,3]], type = "p", pch=20)
           }
         }

       }), las = las, ...)

  lty <- c("solid", "dashed")
  lwd <- c(1)
  if (type == "area") {
    col <- c("black", "black", "white")
  } else {
    col <- c("black")
  }

  ## Recycle values
  lwd <- rep(lwd, len = length(attsurfs))
  lty <- rep(lty, len = length(attsurfs))
  col <- rep(col, len = length(attsurfs))

  for (k in 1:length(attsurfs)) {
    tmp <- attsurfs[[k]]
    lines(tmp[,1], tmp[,2], type = "s", lty = lty[k], lwd = lwd[k], col = col[k])
  }

  mtext(title, 1, line=3.5, cex=par("cex.lab"), las = 0, font = 2)
  box()
}

eafdiffplot <-
  function(data.left, data.right,
           # FIXME: This could be constructed automatically from a nintervals argument.
           intervals = c("[0.0, 0.2)","[0.2, 0.4)","[0.4, 0.6)", "[0.6, 0.8)","[0.8, 1.0]"),
           # FIXME: This could also be constructed automatically by dividing the black-white interval.
           col = c("#FFFFFF", "#BFBFBF","#808080","#404040","#000000"),
           percentiles = c(50),
           full.eaf = FALSE,
           type = ifelse(full.eaf, "area", "point"),
           legend.pos = ifelse(full.eaf, "bottomleft", "topright"),
           title.left = deparse(substitute(data.left)),
           title.right = deparse(substitute(data.right)),
           xlim = NULL, ylim = NULL,
           cex = par("cex"), cex.lab = par("cex.lab"), cex.axis = par("cex.axis"),
           maximise = c(FALSE, FALSE),
           ...)
{
  type <- match.arg (type, c("point", "area"))
  if (length(intervals) < length(col)) {
    stop ("Less intervals than colors: length(intervals) < length(col)")
  } else if (length(intervals) > length(col)) {
    stop ("More intervals than colors: length(intervals) > length(col)")
  }
  title.left <- title.left
  title.right <- title.right
  maximise <- as.logical(maximise)

  if (!is.null(xlim) && maximise[1]) xlim <- -xlim 
  if (!is.null(ylim) && maximise[2]) ylim <- -ylim
  
  data.left[, which(maximise)] <- -data.left[, which(maximise)]
  data.right[, which(maximise)] <- -data.right[, which(maximise)]
  
  attsurfs.left <- compute.eaf.as.list (data.left, percentiles)
  attsurfs.right <- compute.eaf.as.list (data.right, percentiles)

  # Merge the data
  nruns.left <- max(data.left[,3])
  data.combined <- data.right
  data.combined[,3] <- data.combined[,3] + nruns.left
  data.combined <- rbind(data.left, data.combined)

  # FIXME: This can be avoided and just taken from the full EAF below.
  grand.attsurf <- compute.eaf.as.list (data.combined, c(0, 100))
  grand.best <- grand.attsurf[["0"]]
  grand.worst <- grand.attsurf[["100"]]
  if (full.eaf) {
    DIFF <- list()
    if (type == "area") {
      lower.boundaries <- 0:(length(intervals)-1) * 100 / length(intervals)
      DIFF$left <- compute.eaf (data.left, percentiles = lower.boundaries)
      DIFF$right <- compute.eaf (data.right, percentiles = lower.boundaries)
    } else if (type == "point") {
      DIFF$left <- compute.eaf (data.left)
      DIFF$right <- compute.eaf (data.right)
      # Since plot.eafdiff.side uses floor to calculate the color, and
      # we want color[100] == color[99].
      DIFF$left[DIFF$left[,3] == 100,3] <- 99
      DIFF$right[DIFF$right[,3] == 100,3] <- 99
    }
    DIFF$left[,3] <- DIFF$left[,3] * length(intervals) / 100
    DIFF$right[,3] <- DIFF$right[,3] * length(intervals) / 100
    #remove(data.left,data.right,data.combined) # Free memory?
  } else {
    #remove(data.left,data.right) # Free memory?
    DIFF <- compute.eafdiff (data.combined, intervals = length(intervals))
    #remove(data.combined) # Free memory?
  }

  if (is.null(xlim)) {
    xlim <- range(c(grand.best[,1], grand.worst[,1],
                    DIFF$left[,1], DIFF$right[,1]))
  }

  if (is.null(ylim)) {
    ylim <- range(c(grand.best[,2], grand.worst[,2],
                    DIFF$left[,2], DIFF$right[,2]))
  }

  layout(matrix(c(1,1,2,2), 2, 2, byrow = FALSE))
  bottommar <- 5
  topmar   <- 4
  leftmar  <- 4
  rightmar <- 4

  # cex.axis is multiplied by cex, but cex.lab is not.
  op <- par(cex = cex, cex.lab = cex.lab, cex.axis = cex.axis
            , mar = c(bottommar, leftmar, topmar, 0)
            , lab = c(10,5,7)
            , las = 0
            )

  if (maximise[1]) xlim <- range(-xlim)
  if (maximise[2]) ylim <- range(-ylim)

  DIFF$left[, which(maximise)] <- - DIFF$left[, which(maximise)]
  DIFF$right[, which(maximise)] <- - DIFF$right[, which(maximise)]

  attsurfs <- c(list(grand.best), attsurfs.left, list(grand.worst))
  attsurfs <- lapply (attsurfs, function (x)
                      { x[, which(maximise)] <- -x[, which(maximise)]; x })

  .plot.eafdiff.side (DIFF$left,
                     attsurfs = attsurfs,
                     col = col,
                     type = type, full.eaf = full.eaf,
                     title = title.left,
                     xlim = xlim, ylim = ylim,
                     side = "left", maximise = maximise, ...)

  if (nchar(legend.pos) > 0 && !(legend.pos %in% c("no", "none"))) {
    legend(x = legend.pos, y = NULL,
           rev(intervals), rev(col),
           bg = "white", bty = "n", xjust=0, yjust=0, cex=0.9)
  }

  par(mar = c(bottommar, 0, topmar, rightmar))

  attsurfs <- c(list(grand.best), attsurfs.right, list(grand.worst))
  attsurfs <- lapply (attsurfs, function (x)
                      { x[, which(maximise)] <- -x[, which(maximise)]; x })

  .plot.eafdiff.side (DIFF$right,
                      attsurfs = attsurfs,
                      col = col,
                      type = type, full.eaf = full.eaf,
                      title = title.right,
                      xlim = xlim, ylim = ylim,
                      side = "right", maximise = maximise, ...)
  par(op)
}

### Local Variables:
### ess-indent-level: 2
### ess-continued-statement-offset: 2
### ess-brace-offset: 0
### ess-expression-offset: 4
### ess-else-offset: 0
### ess-brace-imaginary-offset: 0
### ess-continued-brace-offset: 0
### ess-arg-function-offset: 2
### ess-close-brace-offset: 0
### indent-tabs-mode: nil
### ess-fancy-comments: nil
### End:
