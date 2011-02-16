\name{SPEA2relativeRichmond}
\alias{SPEA2relativeRichmond}
\docType{data}
\title{
Results of SPEA2 with relative time-controlled triggers on Richmond
water network.
}
\description{
Results of SPEA2 with relative time-controlled triggers on Richmond
water network.
The data has the only goal of providing an example of use of eafplot.
}
\usage{data(SPEA2relativeRichmond)}
\format{
  A data frame as produced by \code{\link{read.data.sets}}.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Manuel \enc{López-Ibáñez}{Lopez-Ibanez}. Operational Optimisation of Water Distribution
  Networks. PhD thesis, School of Engineering and the Built Environment,
  Edinburgh Napier University, UK, 2009.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(HybridGA)
data(SPEA2relativeRichmond)
eafplot (SPEA2relativeRichmond, percentiles = c(25, 50, 75),
         xlab = expression(C[E]), ylab = "Total switches",
         xlim = c(90, 140), ylim = c(0, 25),
         extra.points = HybridGA$richmond, extra.lty = "dashed",
         extra.legend = "Hybrid GA")
}
\keyword{datasets}
