\name{flatQQNorm}
\alias{flatQQNorm}
\alias{flatQQNormY}
\title{ QQ difference plot }
\description{
standard normal transforms to horizontal line at 0
}
\usage{
flatQQNorm(y)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ numeric vector }
}
\value{
list with elements x and y describing the trace of the
qq differrence plot
}

\examples{
set.seed(1234)
plot(flatQQNorm(rnorm(40)),ylim=c(-3,3),ylab="QQNorm - line of identity ")
}
\keyword{ models }

