\name{eafplot}
\alias{eafplot}
\alias{eafplot.default}
\alias{eafplot.data.frame}
\alias{eafplot.formula}
\alias{eafplot.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Empirical Attainment Function for two objectives 
}
\description{
  %% ~~ A concise (1-5 lines) description of what the function does. ~~
  Computes and plots the Empirical Attainment Function, either as
  attainment surfaces for certain percentiles or as points.  This
  function can be used to plot random sets of points like those obtained
  by different runs of biobjective stochastic optimization algorithms.
  An EAF curve represents the boundary separating points that are known
  to be attainable (that is, dominated in Pareto sense) in at least a
  fraction (quantile) of the runs from those that are not. The median
  EAF represents the curve where the fraction of attainable points is
  50\%.  In single objective optimization the function can be used to
  plot the profile of solution quality over time of a collection of runs
  of a stochastic optmizer.  } \usage{ eafplot(x, ...)
%
eafplot.default(x, sets = NULL, groups = NULL,
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
            extra.points = NULL, extra.legend = NULL,
            extra.pch = c(4:25),
            extra.lwd = 0.5,
            extra.lty = "dashed",
            extra.col = "black",
            ymaximise = FALSE,
            xaxis.side = 1, yaxis.side = 2,
            ... )
%
eafplot.data.frame(x, y = NULL, ...)
%
eafplot.formula(x, data = data.frame(), groups = NULL, subset = NULL,
            percentiles = NULL, include.extremes = FALSE, ...)
%
eafplot.list(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{either a matrix of data values, or a data frame, or a list of
    data frames of exactly three columns, or a
    a formula. Formulas of the type: \code{time + cost ~ run | instance}
    will draw \code{time} on the x-axis and \code{cost} on the y-axis. If \code{instance} is
    present the plot is conditional to the instances. 
}
  \item{y}{either a matrix of data values, or a data frame.
}
  \item{data}{a data frame containing the fields mentioned in the formula and in groups.
}
  \item{groups}{this may be used to plot profiles of different algorithms on the same plot.
}
  \item{subset}{
a vector indicating which rows of the data should be used. If left to default \code{NULL} all data in the data frame are used. 
}
\item{include.extremes}{
whether the plot has to include x and y extremes. These points are weakly dominated and hence they would not be plotted with the default \code{include.extremes=FALSE}. 
}
\item{sets}{
  a vector indicating which set each point belongs to.
}
  \item{percentiles}{
a vector indicating which percentile should be plot. The default is to plot only the median attainment curve.
}
\item{attsurfs}{
  if different than NULL, a list of attainment surfaces as ????
}

  \item{type}{string giving the type of plot desired.  The
    following values are possible, \samp{points} and \samp{area}.
}

\item{xlab, ylab, xlim, ylim, log,
    col, lty,  lwd, pch, cex.pch, las}{ graphical
    parameters, see \code{\link{plot.default}}.}

  \item{legend.pos}{the position of the legend, see \code{\link{legend}}.
}
  \item{legend.txt}{a character or expression vector to appear
    in the legend. If \code{NULL}, appropriate labels will be generated.
}

\item{extra.points}{a list of matrices or data.frames with
  two-columns. Each element of the list defines a set of points, or
  lines if one of the columns is \code{NA}.}

\item{extra.pch, extra.lwd, extra.lty, extra.col}{control the graphical
  aspect of the points. See \code{\link{points}} and
  \code{\link{lines}}.  }

\item{extra.legend}{is a character vector providing labels for the
  groups of points.
}
\item{ymaximise}{whether the y-axis denotes a maximisation problem. The
  EAF is calculated accordingly, and the y-axis is reversed for plotting. 
}
  \item{xaxis.side}{ controls on which side that xaxis is drawn. See \code{\link{axis}}.
}
  \item{yaxis.side}{controls on which side that yaxis is drawn. See \code{\link{axis}}.
}

\item{\dots}{Other graphical parameters to \code{\link{plot.default}}.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
No value is returned.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
%   %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{read.data.sets}}, \code{\link{eafdiffplot}}
}
\examples{
data(gcp2x2)
tabucol<-subset(gcp2x2, alg!="TSinN1")
tabucol$alg<-tabucol$alg[drop=TRUE]
eafplot(time+best~run,data=tabucol,subset=tabucol$inst=="DSJC500.5")

eafplot(time+best~run|inst,groups=alg,data=gcp2x2)
eafplot(time+best~run|inst,groups=alg,data=gcp2x2,
	percentiles=c(0,50,100),include.extremes=TRUE,
	cex=1.4, lty=c(2,1,2),lwd=c(2,2,2),
        col=c("black","blue","grey50"))

A1 <- read.data.sets(file.path(system.file(package = "eaf"), "extdata", "ALG_1_dat"))
A2 <- read.data.sets(file.path(system.file(package = "eaf"), "extdata", "ALG_2_dat"))
eafplot(A1, A2, percentiles = c(50))
eafplot(list(A1 = A1, A2 = A2), percentiles = c(50))
\dontrun{dev.copy2pdf(file = "eaf.pdf", onefile = TRUE, width = 5, height = 4)}

## Using extra.points
data(HybridGA)
data(SPEA2relativeVanzyl)
eafplot(SPEA2relativeVanzyl, percentiles = c(25, 50, 75), 
        xlab = expression(C[E]), ylab = "Total switches", xlim = c(320, 400),
        extra.points = HybridGA$vanzyl, extra.legend = "Hybrid GA")

data(SPEA2relativeRichmond)
eafplot (SPEA2relativeRichmond, percentiles = c(25, 50, 75),
         xlab = expression(C[E]), ylab = "Total switches",
         xlim = c(90, 140), ylim = c(0, 25),
         extra.points = HybridGA$richmond, extra.lty = "dashed",
         extra.legend = "Hybrid GA")

data(SPEA2minstoptimeRichmond)
SPEA2minstoptimeRichmond[,2] <- SPEA2minstoptimeRichmond[,2] / 60
eafplot (SPEA2minstoptimeRichmond, xlab = expression(C[E]),
         ylab = "Minimum idle time (minutes)",
         las = 1, log = "y", ymaximise = TRUE, main = "SPEA2 (Richmond)")
}
\keyword{graphs}