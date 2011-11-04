\name{makeCandmat.raw}
\alias{makeCandmat.raw}
\alias{rmixnorm}
\alias{pmix1n}
\alias{dmix1n}
\alias{qmix1n}
\alias{pmix2n}
\alias{dmix2n}
\alias{qmix2n}
\title{create and store reference distributions for edd   }
\description{
create and store reference distributions for edd
}
\usage{
makeCandmat.raw (nPerRow=20, nRowPerCand=20, dists=
  eddDistList, centerScale=TRUE)
}
\arguments{
  \item{nPerRow}{ size of each reference sample}
  \item{nRowPerCand}{ number of samples per candidate}
  \item{dists}{ list of eddDist objects specifying
reference distributions }
  \item{centerScale}{ logical indicating that
simulated data should be centered and scaled}
}
%\details{
%includes various utilities for quantiles of two gaussian mixtures
%}
\value{
matrix with rows realized from reference distributions
}
%\references{  }
\author{VJ Carey <stvjc@channing.harvard.edu>  }
%\note{  }


%\seealso{  }

\examples{
makeCandmat.raw(5,2,eddDistList[1:3])
}
\keyword{ models }

