\name{concatMods}
\alias{concatMods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ render model outputs for comparison }
\description{
render a collection of related models
}
\usage{
concatMods(x, mstub = "mod", type = c("se", "p")[1], dig = 3, clean=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ list of model objects, currently lm, geese, alr, yags objects
supported }
  \item{mstub}{  column header token to prepend 1, 2, ...
for column numbering }
  \item{type}{ "se" means coefficients
reported over (SE=...), "p" means coefficients
reported over (p=...) }
  \item{dig}{ number of digits for report rounding }
  \item{clean}{logical -- if TRUE, removes apparently needless NA from output}
}
\details{
The bottom part of the table is generated by a model-class-specific
\code{getExtra}
method.  See the vignette for details on how to modify this. 
}
\value{
a character matrix
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Vince Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }


\seealso{ \code{\link[xtable]{xtable}} }
\examples{
data(stackloss)
lm1 <- lm( stack.loss ~ Air.Flow, data=stackloss)
lm2 <- lm( stack.loss ~ Water.Temp, data=stackloss)
lm3 <- lm( stack.loss ~ Acid.Conc., data=stackloss)
lm4 <- lm( stack.loss ~ ., data=stackloss)
library(cremo)
print(concatMods(list(lm1, lm2, lm3, lm4)))
}
\keyword{ models }% at least one, from doc/KEYWORDS
