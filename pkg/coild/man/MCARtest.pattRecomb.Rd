\name{MCARtest.pattRecomb}
\alias{MCARtest.pattRecomb}
\title{a permutation-based test of MCAR inspired by Chen-Little recombination
estimator}
\description{
Chen-Little (1999) propose a missing pattern-based decomposition of
GEEs in incomplete datasets and provide chi-squared tests of
the hypothesis of MCAR.  This function computes scores
on elements of their decomposition and evaluates a
bootstrap test of MCAR.
}
\usage{
MCARtest.pattRecomb(formula, idname="id", family=gaussian, corstr="independence", stratDf, yname="y", xnames=c("x1", "x2"), timesname="times", stratName="stratum", sizeFun=geeCoSize, B=20)
}
\arguments{
  \item{formula}{ formula -- model formula for GEE}
  \item{idname}{ idname -- name of variable defining clusters}
  \item{family}{ family -- glm family}
  \item{corstr}{ corstr -- gee corstr}
  \item{stratDf}{ stratDf -- data frame}
  \item{yname}{ yname -- name of outcome variable}
  \item{xnames}{ xnames -- vector of names of covariates required}
  \item{timesname}{ timesname -- name of sampling time variable}
  \item{stratName}{ stratName -- names of missing pattern stratum discriminator}
  \item{sizeFun}{ sizeFun -- function on regression coefficients}
  \item{B}{ B -- number of permutations to be generated}
}
\details{
}
\value{
}
\references{  }
\author{  }
\note{  }


\seealso{  }

\examples{
data(demoRC)
require(MASS)
gsf <- coSizeMaker(2:3) # focus attention on parms 2 and 3

MCARtest.pattRecomb( edasth ~ ARM1+ARM2, idname="id", 
       yname="edasth", xnames=c("ARM1","ARM2"), timesname="YRSONST", 
       stratName="class", stratDf=demoRC, B=10, sizeFun=gsf)
}
\keyword{ models }

