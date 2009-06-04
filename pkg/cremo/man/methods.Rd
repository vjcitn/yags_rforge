\name{cremoMethods}
\alias{getVarn}
\alias{getP}
\alias{getSE}
\alias{getCoef}
\alias{getExtra}
\alias{getVarn,alr-method}
\alias{getP,alr-method}
\alias{getSE,alr-method}
\alias{getCoef,alr-method}
\alias{getExtra,alr-method}
\alias{getVarn,gee-method}
\alias{getP,gee-method}
\alias{getSE,gee-method}
\alias{getCoef,gee-method}
\alias{getExtra,gee-method}
\alias{getVarn,geese-method}
\alias{getP,geese-method}
\alias{getSE,geese-method}
\alias{getCoef,geese-method}
\alias{getExtra,geese-method}
\alias{getVarn,lm-method}
\alias{getP,lm-method}
\alias{getSE,lm-method}
\alias{getCoef,lm-method}
\alias{getExtra,lm-method}
\alias{getVarn,yagsResult-method}
\alias{getP,yagsResult-method}
\alias{getSE,yagsResult-method}
\alias{getCoef,yagsResult-method}
\alias{getExtra,yagsResult-method}
\alias{getVarn,clogit-method}
\alias{getP,clogit-method}
\alias{getSE,clogit-method}
\alias{getCoef,clogit-method}
\alias{getExtra,clogit-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ render model outputs for comparison }
\description{
support methods render a collection of related models
}
\usage{
getVarn(x)
getP(x)
getSE(x)
getCoef(x)
getExtra(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ instance of R model fit class }
}
\details{
These generic methods return information appropriate
for lm, gee, coxph models for use in \code{concatMods}.
See the vignette for details on how to override these methods. 
}
\value{
a numeric vector
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Vince Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }


\seealso{ \code{\link[xtable]{xtable}} }
\examples{
data(stackloss)
lm1 <- lm( stack.loss ~ Air.Flow, data=stackloss)
getCoef(lm1)
getExtra(lm1)
}
\keyword{ models }% at least one, from doc/KEYWORDS
