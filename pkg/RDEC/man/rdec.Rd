\name{rdec}
\alias{rdec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ regression with damped exponential correlation }
\description{
ML fitting of a multivariate gaussian model for longitudinal data
that generalizes AR(1), CS, and MA(1) covariance structures
}
\usage{
rdec(formula, id, S, data = sys.parent(), subset, 
   na.action = na.fail, omega.init = c(0, 0), 
   omega.low = c(0.001, 0), omega.high = c(0.95, 0.95), 
   ltol = 0.01, contrasts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ standard linear modeling formula }
  \item{id}{ vector of cluster discriminators }
  \item{S}{ vector of observation times }
  \item{data}{ resolve model formula elements in this data frame }
  \item{subset}{ logical subsetting vector }
  \item{na.action}{ action to take for NA in formula element }
  \item{omega.init}{ starting value for covariance parameter vector (dim=2) }
  \item{omega.low}{ lower bound for covariance parameter vector }
  \item{omega.high}{ upper bound for covariance parameter vector }
  \item{ltol}{ tolerance for convergence of likelihood function }
  \item{contrasts}{ probably not used }
  \item{verbose}{ print -2 likelihood value per iteration }
}
\details{
}
\value{
a structure with model fitting results
}
\references{ Munoz Carey Rosner et al, Biometrics 1992 }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
\note{ }

\seealso{ gee, geepack, nlme, reshape}
\examples{
data(rdtest)
summary( rdec(rdec.y~rdec.x, id=rdec.id, S=rdec.s, omega.init = c(.5,.5), omega.low=c(.01,.01),
   omega.high=c(.95, .95), data=rdtest ))
}
\keyword{ models }% at least one, from doc/KEYWORDS
