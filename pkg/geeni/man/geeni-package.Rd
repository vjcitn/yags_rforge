\name{geeni-package}
\alias{geeni-package}
\alias{geeni}
\docType{package}
\title{
fit a gee in bounded memory, data from disk etc. -- demonstration package
}
\description{
illustrate concepts of deployable statistical methodology
using GEE solution -- NOTE: will not run "out of box" -- you need to
run the buildff vignette, which will drop something in the installed package image,
that allows the rungee vignette to work.
}
\details{
\tabular{ll}{
Package: \tab geeni\cr
Version: \tab 0.0.0\cr
Suggests: \tab multicore\cr
Depends: \tab ff\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
Built: \tab R 2.13.0; ; 2010-11-30 11:54:39 UTC; unix\cr
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{buildff} \tab build an ff representation of grouped data instance (source, pdf)\cr
\code{rungee} \tab using a gdManager to deploy GEE (source)\cr
}

After the buildff vignette is run by the individual who installed the package an ff
archive will be generated in the mgrs folder of the installed package, which is
employed by the code in the rungee vignette

}
\author{
VJ Carey <stvjc@channing.harvard.edu>

Maintainer: Vince Carey <stvjc@channing.harvard.edu>
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
related in spirit to Lumley's biglm, bigglm, but statistics are incremented by data groups
rather than records
}
\examples{
# this just shows the basic solution process
pargee
}
