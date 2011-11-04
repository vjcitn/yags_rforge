
\name{centerScale}
\alias{centerScale}
\alias{ctr}
\title{center and scale a vector to zero median, unit mad  }
\description{center and scale a vector to zero median, unit mad  }
\usage{
centerScale(x)
}
\arguments{
  \item{x}{ a numeric vector}
}
%\details{
%}
\value{ a shifted and scaled version of x with zero median, unit mad
}
%\references{  }
\author{VJ Carey <stvjc@channing.harvard.edu>  }
%\note{  }
%
%
%\seealso{  }

\examples{
summary(centerScale(runif(200)))
}
\keyword{ models }
