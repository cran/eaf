\name{read.data.sets}
\alias{read.data.sets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read several data.frame sets
}
\description{
     Reads a text file in table format and creates a data frame from it. The
     file may contain several sets, separated by empty lines. The function adds an
     additional column \code{set} to indicate to which set each row
     belongs.
}
\usage{
read.data.sets(file, col.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
    the name of the file which the data are to be read from.
    Each row of the table appears as one line of the file.  If it
    does not contain an \emph{absolute} path, the file name is
    \emph{relative} to the current working directory, \code{getwd()}.
    Tilde-expansion is performed where supported.
}
\item{col.names}{
  a vector of optional names for the variables.  The default is to use \samp{"V"} followed by the column number.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A data frame (\code{data.frame}) containing a representation of the
  data in the file. An extra column \code{set} is added to indicate to
  which set each row belongs. 
}
%\references{ 
%% ~put references to the literature/web site here ~
%}
\author{
  Manuel \enc{López-Ibáñez}{Lopez-Ibanez}
}
\note{%
  There are several examples of data sets in \code{file.path(system.file(package="eaf"),"extdata")}. 
}

\section{Warning}{
  A know limitation is that the input file must use newline characters
  native to the host system, otherwise they will be, possibly silently,
  misinterpreted. In GNU/Linux the program \code{dos2unix} may be used
  to fix newline characters.
}

\seealso{
\code{\link{read.table}}, \code{\link{eafplot}}, \code{\link{eafdiffplot}}
}
\examples{
A1<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_1_dat"))
str(A1)
A2<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_2_dat"))
str(A2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
