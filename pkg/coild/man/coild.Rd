\name{coild}
\alias{coild}
\alias{demoNRfun}
\alias{demoSubjectRealizer}
\alias{evalSlot}
\alias{.initClasses}
\alias{is.monomiss.mat}
\alias{is.monomiss.vec}
\alias{last.only}
\alias{lmisspatt}
\alias{makeC1Reg}
\alias{make.monotone.vec}
\alias{WuFollman4}
\alias{demod}
\alias{ctimes}
\alias{demod.id}
\alias{demod.times}
\alias{crephalo}
\alias{crephalo.df}
\non_function{}
\title{ Classification Of Incomplete Longitudinal Data }
\description{An S4 class, associated methods and function useful
for exploring data-based evaluations of missing data processes
operating with longitudinal datasets.}
\usage{
#library(coild)
#data(coild)
# sets up the matrices demod and demod.times for K=5 nominal visits
#demoLD <- makeLongimat( mat=demod, times=demod.times )
# instantiates the longimat class
#diggle.test(demoLD, function(x)x, stat=mean, NSIM=19, inseed=123)
# obtains montecarlo p-values for hypothesis that missing values are
# MCAR at K-1 visits
}
\source{
 Diggle Liang and Zeger 1994, sec 11.2
}
\references{
 Diggle Liang and Zeger 1994, sec 11.2
}
\examples{
data(demod)
dim(demod)
data(demod.times)
demod.times
}
\keyword{packages}
\keyword{models}

