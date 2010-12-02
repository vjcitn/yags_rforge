\name{pargee}
\alias{pargee}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
demonstrate a parallelized, out-of-memory based GEE solver
}
\description{
demonstrate a parallelized, out-of-memory based GEE solver
}
\usage{
pargee(gd, family, binit, maxit = 20, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gd}{
instance of \code{\link[geeni]{gdManager-class}}
}
  \item{family}{
GLM stats \code{family} instance, e.g., \code{\link{gaussian}}
}
  \item{binit}{
initial value for regression parameter vector
}
  \item{maxit}{
maximum number of iterations allowed
}
  \item{tol}{
convergence criterion -- will stop when largest maximum absolute size
of any component of the iteration step is smaller than \code{tol}
}
}
\details{
very preliminary demonstration, mildly validated on orthodontic data for gaussian case
}
\value{
regression parameter vector
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
V Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## to obtain multicore performance, please do
# library(doMC)
# registerDoMC([num cores])  # give number of cores you want to use
## before running this example
 if (file.exists(fd <<- system.file("mgrs/flatOrth.rda", package="geeni"))) {
  load(fd)
  pargee(flatOrth, gaussian, c(0,0,0))
  } else {
  cat("you need to run the buildff vignette to update support information in the installed image of geeni")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
