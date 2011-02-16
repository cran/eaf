\name{HybridGA}
\alias{HybridGA}
\docType{data}
\title{
Results of Hybrid GA on vanzyl and Richmond water networks
}
\description{
Results of Hybrid GA on vanzyl and Richmond water networks.
The data has the only goal of providing an example of use of eafplot.
}
\usage{data(HybridGA)}
\format{
  A list with two data frames, each of them with three columns, as
  produced by \code{\link{read.data.sets}}.
  \describe{
    \item{\code{$vanzyl}}{data frame of results on vanzyl network}
    \item{\code{$richmond}}{data frame of results on Richmond
      network. The second column is filled with \code{NA}}
  }
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
print(HybridGA$vanzyl)
print(HybridGA$richmond)
}
\keyword{datasets}
