\name{is.monomiss}
\alias{is.monomiss}
\title{ helper function testing for monotone missingness in a longitudinal dataset }
\description{
checks for monotone missingness (e.g., once a subject is lost to
followup, they do not return) in a matrix representation of a
longitudinal dataset in which each row is a subject's response
sequence and columns represents repeated observations with earliest
observation to the left.  A generic function predicate also 
exists, which applies to elements of class \code{longimat}
}
\usage{
is.monomiss.mat(mat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mat}{ ~~Describe \code{mat} here~~ }
}
\details{
}
\value{
logical value
}
\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey }
\examples{
data(demod)
is.monomiss.mat(demod)
# TRUE
demod[1,1:3] <- c(1,NA,1)
is.monomiss.mat(demod)
# FALSE
}
\keyword{models}

