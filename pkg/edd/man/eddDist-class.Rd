\name{eddDist-class}
\docType{class}
\alias{eddDist-class}
% [1] "B28"              "B82"              "CDFname"          "CS1"
% [5] "eddDistList"      "edd.unsupervised" "genName"          "LN01"
% [9] "Mad"              "med"              "MIXN1"            "MIXN2"
%[13] "N01"              "parms"            "qfName"           "qfun"
%[17] "qmixnorm"         "rmixnorm.alt"     "stub"             "T3"
%[21] "tag"              "testVec"          "U01"
\alias{plotlim}
\alias{qfName}
\alias{qfun}
\alias{stub}
\alias{tag}
\alias{testVec}
\alias{B28}
\alias{T3}
\alias{B82}
\alias{CS1}
\alias{CDFname}
\alias{CDFname,eddDist-method}
\alias{eddDistList}
\alias{genName}
\alias{genName,eddDist-method}
\alias{LN01}
\alias{U01}
\alias{Mad}
\alias{Mad,eddDist-method}
\alias{med}
\alias{med,eddDist-method}
\alias{parms}
\alias{parms,eddDist-method}
\alias{qfName,eddDist-method}
\alias{qfun,eddDist-method}
\alias{stub,eddDist-method}
\alias{tag,eddDist-method}
\alias{latexTag}
\alias{latexTag,eddDist-method}
\alias{plotlim,eddDist-method}
\alias{testVec,eddDist-method}
\alias{testVec,numeric,eddDist,logical-method}
\alias{qmixnorm}
\alias{dmixnorm}
\alias{pmixnorm}
\alias{rmixnorm.alt}
\alias{MIXN1}
\alias{MIXN2}
\alias{N01}
\title{Class "eddDist" }
\description{  objects from this class can be used to
simulate or tabulate reference distributions for edd}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("eddDist", ...)}.
These objects include a simple stub (like "norm", which can
be modified to obtain the name of a generator (prepend "r"),
cdf (prepend "p"), etc.) in R.
}
\section{Slots}{
  \describe{
    \item{\code{stub}:}{Object of class \code{"character"} stub
of a distribution system in R, to which "r" is prepended
to get a generator, "p" to get a cdf, "q" to get a quantile
function... }
    \item{\code{parms}:}{Object of class \code{"numeric"} named
      vector of parameters for a member of the family }
    \item{\code{median}:}{Object of class \code{"numeric"} median
      of the distribution (sometimes has to be computed by simulation) }
    \item{\code{mad}:}{Object of class \code{"numeric"}  MAD
      of the distribution (sometimes has to be computed by simulation)}
    \item{\code{tag}:}{Object of class \code{"character"}  an
      informative character tag }
    \item{\code{latexTag}:}{Object of class \code{"character"}  an
      informative character tag in latex format, which can use
      subscripts, greek letters, etc. }
    \item{\code{plotlim}:}{Object of class \code{"numeric"} an
      informative numeric indicating the upper and lower bounds of the
      distribution. }
  }
}

\section{Methods}{
  \describe{
    \item{CDFname}{\code{signature(x = "eddDist")}: prepend "p" to stub(x) }
    \item{genName}{\code{signature(x = "eddDist")}: prepend "r" to stub(x) }
    \item{Mad}{\code{signature(x = "eddDist")}: accessor }
    \item{med}{\code{signature(x = "eddDist")}: accessor }
    \item{parms}{\code{signature(x = "eddDist")}: accessor }
    \item{qfName}{\code{signature(x = "eddDist")}: prepend "q" to stub(x) }
    \item{qfun}{\code{signature(e = "eddDist")}: construct a quantile
function from the object }
    \item{stub}{\code{signature(x = "eddDist")}: accessor }
    \item{tag}{\code{signature(x = "eddDist")}: accessor }
    \item{latexTag}{\code{signature(x = "eddDist")}: accessor }
    \item{plotlim}{\code{signature(x = "eddDist")}: accessor }
    \item{testVec}{\code{signature(x = "numeric", eddd = "eddDist", is.centered = "logical")}: apply ks.test to a given vector x against the dist specified
by the eddDist components }
  }
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
% or\code{\link{CLASSNAME-class}} for links to other classes }

\examples{
library(edd)
eddDistList[[1]]
qfun(eddDistList[[1]])
}
\keyword{classes}
