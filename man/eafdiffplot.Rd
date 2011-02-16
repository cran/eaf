\name{eafdiffplot}
\alias{eafdiffplot}
\title{
Empirical attainment function differences 
}
\description{
  Plot the differences between the empirical attainment functions of two
  data sets as a two-panel plot, where the left side shows the values of
  the left EAF minus the right EAF and the right side shows the
  differences in the other direction.
}
\usage{
eafdiffplot(data.left, data.right,
             intervals = c("[0.0, 0.2)","[0.2, 0.4)","[0.4, 0.6)",
                           "[0.6, 0.8)","[0.8, 1.0]"),
             col = c("#FFFFFF", "#BFBFBF","#808080","#404040","#000000"),
             percentiles = c(50),
             full.eaf = FALSE,
             type = ifelse(full.eaf, "area", "point"),
             legend.pos = ifelse(full.eaf,"bottomleft", "topright"),
             title.left = deparse(substitute(data.left)),
             title.right = deparse(substitute(data.right)),
             xlim = NULL, ylim = NULL,
             cex = par("cex"),
             cex.lab = par("cex.lab"),
             cex.axis = par("cex.axis"),
             ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.left, data.right}{data frames corresponding to the input
    data of left and right sides, respectively. Each data frame has at least three columns,
    the third one being the set of each point. See also \code{\link{read.data.sets}}.
}
\item{intervals}{a character vector that determines the levels of EAF
  differences that are plotted. \code{length(intervals)} defines in how
  many intervals the EAF differences are partitioned, whereas
  \code{intervals} gives the labels of each interval for the legend.
}
\item{col}{a character vector giving the color for each
  interval. It must be true that \code{length(col) == length(intervals)}.
}
  \item{percentiles}{the percentiles of the EAF of each side that will
    be plotted as  attainment surfaces. See \code{\link{eafplot.default}}.
}
\item{full.eaf}{whether to plot the EAF of each side instead of the
  differences between the EAFs.

}
\item{type}{whether the EAF differences are plotted as points
  (\samp{points}) or whether to color the areas that have at least a
  certain value (\samp{area}).
}
  \item{legend.pos}{the position of the legend. See
    \code{\link{legend}}.
}
\item{title.left, title.right}{title for left and right panels, respectively.

}
\item{xlim, ylim, cex, cex.lab, cex.axis}{
  graphical
    parameters, see \code{\link{plot.default}}.
}
\item{\dots}{other graphical parameters are passed down to
  \code{\link{plot.default}}.
}
}
\details{
  This function calculates the differences between the EAFs of two
  data sets, and plots on the left the differences in favour
  of the left data set, and on the right the differences in favour of
  the right data set. By default, it also plots the grand best and worst
  attainment surfaces, that is, the 0\% and 100\%-attainment surfaces
  over all data. This two surfaces delimit the area where differences
  may exist. In addition, it also plots the 50\%-attainment surface of
  each data set.
}
\value{
  No return value.
}
\references{
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
   \code{\link{read.data.sets}}, \code{\link{eafplot}}
}
\examples{
A1<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_1_dat"))
A2<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_2_dat"))
eafdiffplot(A1,A2)
eafdiffplot(A1,A2, full.eaf=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}


%% Ignore this.
% A1<-gcp2x2[gcp2x2$alg=="TSinN1" & gcp2x2$inst=="DSJC500.5",]
% A2<-gcp2x2[gcp2x2$alg=="SAKempeFI" & gcp2x2$inst=="DSJC500.5",]
% eafdiffplot(A1[,c(5,4,3)],A2[,c(5,4,3)])}
