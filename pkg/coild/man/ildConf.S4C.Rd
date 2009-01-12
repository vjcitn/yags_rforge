\name{ildConf.S4C}
\alias{ildConf.S4C}
\alias{coild}
\title{ configuration of an incomplete longitudinal dataset }
\description{
establishes structure and models underlying
an incomplete longitudinal dataset
}
\usage{
#new("ildConf", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{ number of subjects, all of whom contribute at least
one observation }
  \item{K}{ maximum number of observations per subject }
  \item{subjectModelName}{ name of function that creates
incomplete data for a given subject }
  \item{subjectModelParms}{ arguments to the function
in slot \code{subjectModelName} }
  \item{subjectNonrespFun}{ a function that will be applied
to the complete data generated through slot \code{subjectRealizer}
to create incomplete data }
  \item{monotone}{ logical indicating whether or not missingness
is monotone }
  \item{subjectRealizer}{ function that applies the data generation
and response models to obtain the observed data }
}
\details{
This class exists to formalize the process of generating
incomplete longitudinal data with specified
data and response models.
}
\value{
an object of class \code{ildConf}, to which method \code{sampleFrom}
may be applied to obtain audited simulation series
}
\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey }

\examples{
makeC1Reg <- function(a,b,k,sig)
 {
 x <- 1:k
 a+b*x+rnorm(k,0,sig)
 }
C1parms <- list(a=0,b=1,k=5,sig=.5)

demoNRfun <- function(x,ub=1.5)
 {
 n <- length(x)
 xx <- (x - median(x))/mad(x)
 if (!any(z <- xx[2:n]>ub)) return(x)
 drop <- min((1:n)[c(FALSE,z)])
 x[(drop):n] <- NA
 x
 }

demoSubjectRealizer <- function(this, ...)
 {
 tmp <- do.call( this@subjectModelName, this@subjectModelParms )
 this@subjectNonrespFun( tmp )
 }

C1 <- new("ildConf", n=100, K=5, subjectModelName="makeC1Reg",
 subjectModelParms = C1parms, subjectNonrespFun = demoNRfun,
 subjectRealizer=demoSubjectRealizer)
SC1 <- sampleFrom(C1)

}
\keyword{ missing data }
\keyword{ longitudinal data }
\keyword{models}

