\name{SPEA2minstoptimeRichmond}
\alias{SPEA2minstoptimeRichmond}
\docType{data}
\title{
Results of SPEA2 when minimising electrical cost and maximising the
minimum idle time of pumps on Richmond water network.
}
\description{
Results of SPEA2 when minimising electrical cost and maximising the
minimum idle time of pumps on Richmond water network.
The data has the only goal of providing an example of use of eafplot.
}
\usage{data(SPEA2minstoptimeRichmond)}
\format{
  A data frame as produced by \code{\link{read.data.sets}}. The second
  column measures time in seconds and corresponds to a maximisation problem.
}
% \details{
% %%  ~~ If necessary, more details than the __description__ above ~~
% }
\source{
  Manuel \enc{López-Ibáñez}{Lopez-Ibanez}. Operational Optimisation of Water Distribution
  Networks. PhD thesis, School of Engineering and the Built Environment,
  Edinburgh Napier University, UK, 2009.
}
% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(HybridGA)
data(SPEA2minstoptimeRichmond)
SPEA2minstoptimeRichmond[,2] <- SPEA2minstoptimeRichmond[,2] / 60
eafplot (SPEA2minstoptimeRichmond, xlab = expression(C[E]),
         ylab = "Minimum idle time (minutes)",
         las = 1, log = "y", maximise = c(FALSE, TRUE))
}
\keyword{datasets}
