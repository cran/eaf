\name{SPEA2relativeVanzyl}
\alias{SPEA2relativeVanzyl}
\docType{data}
\title{
Results of SPEA2 with relative time-controlled triggers on Vanzyl's
water network.
}
\description{
Results of SPEA2 with relative time-controlled triggers on Vanzyl's
water network.
The data has the only goal of providing an example of use of eafplot.
}
\usage{data(SPEA2relativeVanzyl)}
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
data(SPEA2relativeVanzyl)
eafplot(SPEA2relativeVanzyl, percentiles = c(25, 50, 75), 
        xlab = expression(C[E]), ylab = "Total switches", xlim = c(320, 400),
        extra.points = HybridGA$vanzyl, extra.legend = "Hybrid GA")
}
\keyword{datasets}
