\name{diggle.test}
\alias{diggle.test}
\title{test for missingness `completely at random'}
\description{
For an N x K data structure with K nominal observations,
obtain the K-1 p-values corresponding to the hypotheses that responses among
subjects about to drop out at time k<K are a random sample from
those among subjects who remain on study at time k.
}
\usage{
diggle.test(obj, scorefn, stat=mean, NSIM=19, inseed=123, 
  method=c("MC", "enum")[1], extreme=c("low", "high")[1],
  require.monotone=TRUE )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{obj}{ an object of class longimat, assumed to represent
nominally balanced incomplete longitudinal data with monotone missingness }
  \item{scorefn}{ vector-to-scalar function (denoted h in
DLZ 1994) measuring contribution
of observed outcomes to prediction of missingness.  The R
function \code{sum} can be used if all observed data up to
time k are to be weighted equally in predicting missingness
at time k+1.  A function that puts weight 0 on all responses
prior to visit k-1 and weight 1 
(\code{function(x) {wt <- rep(0,length(x)); wt[length(x)]<-1; wt*x}}
can be used to test for MCAR under an assumption that response
values just prior to dropout are informative.}
that
  \item{stat}{ function computing the statistic to be evaluated
on the values of h for imminent dropouts and
non dropouts.  \code{mean} is default. }
  \item{NSIM}{ MC tests of randomness involve NSIM simulations }
  \item{inseed}{ PRNG seed will be initialized with the
value of this  parameter }
  \item{method}{ one of "enum" or "MC" }
  \item{extreme}{ "high" or "low" }
  \item{require.monotone}{ logical -- currently needs to be TRUE!}
}
\details{
%Currently obtains a two-tailed MC-based p-value for each visit
%k=1, ..., K-1.  Probably should be one-tailed.
}
\value{
returns a vector of up to K-1 method-based %two-tailed 
p-values.
}
\references{ DLZ 1994, 11.2 ,
Diggle PJ, Testing for random dropouts in repeated
measures data, Biometrics v45, 1989, 1255-58.}
\author{ VJ Carey }

\examples{
data(demod)
data(demod.times)
demoLD <- makeLongimat( mat=demod, times=demod.times )
p1 <- diggle.test( demoLD, sum )
p1p <-  diggle.test( demoLD, sum, inseed=432 )
#
# reproduce the Biometrics 1989 results
data(crephalo)
data(ctimes)
CREP1 <- makeLongimat( crephalo[1:11,], ctimes[1:11,] )
CREP2 <- makeLongimat( crephalo[12:21,], ctimes[12:21,] )
CREP3 <- makeLongimat( crephalo[22:32,], ctimes[1:11,] )
CREP4 <- makeLongimat( crephalo[33:43,], ctimes[1:11,] )
CREP5 <- makeLongimat( crephalo[44:54,], ctimes[1:11,] )
p1 <- diggle.test( CREP1, last.only, extr="high", meth="enum" )
p2 <- diggle.test( CREP2, last.only, extr="high", meth="enum" )
p3 <- diggle.test( CREP3, last.only, extr="high", meth="enum" )
p4 <- diggle.test( CREP4, last.only, extr="high", meth="enum" )
p5 <- diggle.test( CREP5, last.only, extr="high", meth="enum" )
allp <- na.omit(c(p1,p2,p3,p4,p5,.95)) # stick on
           # final p-value in group 5 (which
           # had no completers) to get 16 p-values
ks.test(allp,"punif",alternative="gr")
}
\keyword{ longitudinal data }
\keyword{ diagnostics }
\keyword{models}

