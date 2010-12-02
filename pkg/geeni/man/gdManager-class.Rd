\name{gdManager-class}
\Rdversion{1.1}
\docType{class}
\alias{gdManager-class}
\alias{getGrp,gdManager,numeric-method}

\title{Class "gdManager"}
\description{
illustrative class for a data provider based on nlme groupedData, but with an ff back end
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gdManager", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{discrim}:}{Object of class \code{"ff_vector"} a representation of the
run length encoding of the grouped structure label, which changes values as we traverse
records from different groups, like the `id' variable often used in GEE or NLME applications}
    \item{\code{numdat}:}{Object of class \code{"ff_matrix"} numerical data for modeling,
stored in ff matrix format}
    \item{\code{formula}:}{Object of class \code{"formula"} used to describe the nlme
groupedData structure}
  }
}
\section{Methods}{
  \describe{
    \item{getGrp}{\code{signature(gd = "gdManager", ind = "numeric")}: ... }
returns the numerical data for a specified group; ind is a sequential index.
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
Vince Carey <stvjc@channing.harvard.edu>
}
\note{
This is very experimental.  We want to extend the \code{data=} pattern often used for
R modeling functions, where it is expected that an R data frame is bound to \code{data}.
Instead we want to think of a data provider, that will provide numerical data for
analysis upon request.  For GEE, the natural type of request is to acquire the data
for a cluster of observations.  The \texttt{getGrp} method accomplishes this.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("gdManager")
}
\keyword{classes}
