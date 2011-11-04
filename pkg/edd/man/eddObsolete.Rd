\name{eddObsolete}
\alias{eddObsolete}
\alias{fq.matrows}
\alias{makeCandmat.theor}
\alias{s.rmix1norm}
\alias{s.rmix2norm}
\alias{mkt}
\alias{testcl}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Expression Density Diagnostics}
\description{
Classify cohort distributions of gene expression values.
}
\usage{
eddObsolete(eset, 
   ref=c("multiCand", "uniCand", "test", "nnet")[1], 
   k=10, l=6, nnsize=6, nniter=200)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{eset}{instance of Biobase class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.} 
  \item{ref}{one of 'multiCand', 'uniCand', 'test' or 'nnet'. see details.}
  \item{k}{k setting for knn -- number of nearest neighbors to poll.}
  \item{l}{l setting for knn -- minimum number of concordant assents.}
  \item{nnsize}{size parameter for nnet.}
  \item{nniter}{iter setting for nnet.}
}
\details{
Four options are available for classifying expression densities.
Data on each gene are shifted and scaled to have median zero and mad 1.
They are then compared to shapes of reference distributions (standard Gaussian,
chisq(1), lognorm(0,1), t(3), .75N0,1+.25N4,1, .25N0,1+.75N4,1, Beta(2,8), Beta(8,2),
U(0,1)) after each of these has been transformed to have median 0 and mad 1.
Classification proceeds by one of four methods, selected by setting of
the 'ref' argument.  Suppose there are S samples in the ExpressionSet.

multiCand -- 100 samples of size S are drawn from each reference
distribution and then scaled to med 0, mad 1.  The knn(k,l) procedure
is used to classify the genes based on proximity to representatives
in this set.

uniCand -- one representative of size S is created from each reference
distribution, using the theoretical quantiles. knn(1,0) is used
to classify genes based on proximity to these representatives.

test -- classification of each gene is based on maximum p-value
of Kolmogorov-Smirnov tests vs each reference distribution.  If
the p-value never exceeds .1, 'doubt' is declared.

nnet -- 100 samples of size S are drawn from each reference distribution
and then scaled to med 0, mad 1.  A neural net is fit to this dataset
and the associated labels.  The net is then applied to the scaled gene
expression data and the predictions are used for classification.
}
\value{
the vector of classifications, with NAs for nonclassifiable genes
}
\author{VJ Carey}
\examples{
require(Biobase)
data(sample.ExpressionSet)
print(summary(eddObsolete(sample.ExpressionSet,k=10,l=2)))

# 6 x 20 x 50 test problem
set.seed(1234)
test <- matrix(NA,nr=120,nc=50)
test[1:20,] <- rnorm(1000)
test[21:40,] <- rt(1000,3)
test[41:60,] <- rexp(1000,4)
test[61:80,] <- rmixnorm(1000,.750,0,1,4,1)
test[81:100,] <- runif(1000)
test[101:120,] <- rlnorm(1000)
labs <- c(rep("n01",20),rep("t3",20),
rep("exp",20),rep("mix1",20),rep("u01",20),rep("ln01",20))

phenoData <- new("AnnotatedDataFrame")
pData(phenoData) <- data.frame(1:50)
varLabels(phenoData) <- list("Col1")
TT <- new("ExpressionSet", exprs=test, phenoData = phenoData)

multrun <- eddObsolete(TT, k=10, l=2)
print(table(given=labs, multiCand=multrun))
netrun <- eddObsolete(TT, ref="nnet")
print(table(given=labs, netout=netrun))
newrun <- edd(TT, meth="nnet", size=10, decay=.2)
print(table(given=labs, newout=newrun))
newrun <- edd(TT, meth="test")
print(table(given=labs, newout=newrun))
}
\keyword{models}% at least one, from doc/KEYWORDS
