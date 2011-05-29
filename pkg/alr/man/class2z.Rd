\name{class2z}
\alias{class2z}
\title{
Pairwise log OR regression design generator
}
\description{
Regression estimation for clustered binary outcomes
can be facilitated by regression analysis of pairwise
intracluster log odds ratios; such regressions form
the ``extension'' to extended GEE methods GEE2
(Liang, Zeger, Qaqish, 1992) and ALR (Carey, Zeger, Diggle, 1993).
Specification of pairwise log odds ratio regression designs
can be cumbersome; this program assists in this task.
There are assumed to be N observations in all, falling
into M = length(unique(id)) clusters; each observation
is a member of one of k classes.
}
\usage{
class2z(cvec, id, k, dmat)
~move the above line to just above the first optional argument
}
\arguments{
\item{cvec}{
class vector: an N vector with values in 1:k; cvec[i] is the
class in which the ith outcome lies
}
\item{id}{
cluster vector: an N vector with M distinct values
}
\item{k}{
the number of classes
}
\item{dmat}{
design generator.  Let q denote the number of terms
on the right hand side of the pairwise log OR regression
equation, inclusive of intercept.  dim(dmat) = c(n*(n+1)/2,q).
The rows of dmat are in lexicographic order corresponding to
non-redundant enumeration of
pairs:
(1,1), (1,2), \dots , (1,k), (2,2), \dots, (2,k), \dots ,(k,k), and
the rows of this matrix
describe how the q regression parameters are used to describe
association for pairs of elements.  Specifically, let t denote
the row of dmat corresponding to the class-pair (r,s), and let
Y1 be the outcome from a member of class r and Y2 be the outcome
from a member of class s.  Then log odds ratio(Y1,Y2)=a1*dmat[t,1]+\dots+
aq*dmat[t,q].  See example below.
}}
\value{
a self-describing list with components z, zid, dmat, and flag.
}
\details{
A simple example: a family of size 4 with 2 parents and
2 children; cvec = c(1,1,2,2), k=2, dmat = rbind(c(1,0),c(1,1),c(1,1)).
The dmat expresses the idea that first-order relations (parent-child,
sibling-sibling) have equivalent pairwise associations in the outcome,
distinct from the spouse-spouse association.  The spouse-spouse
log OR is then alpha[1], whilst the first-order relative pairwise
log OR is alpha[1]+alpha[2].
}
\references{
Heagerty, P, Ph.D. diss, Dept Biostat, Johns Hopkins University
School of Public Health.


Heagerty P, Zeger S, "Marginal regression models for clustered ordinal measurements",
JASA Sept 1996, v91, n435, 1024--1036.


Carey V, Zeger S, Diggle P, "Modelling multivariate binary data with alternating
logistic regressions", Biometrika 1993 p517.
}
\seealso{
alr
}
\examples{
# 12 observations in clusters of size 3, 5, 4
id <- c(rep(2,3),rep(4,5),rep(6,4))
# each observation is in one of two classes
k <- 2
cvec <- c(1,1,2,1,1,2,2,2,1,1,2,2)
# alpha1 is LPWOR for outcomes in the same cluster
# when both are in class 1
# alpha1+alpha2 = LPWOR for pairs with
# one element from class 1 and one from class 2
# alpha1+alpha3 = LPWOR for pairs with both elements
# from class 2
dmat <- rbind(c(1,0,0),c(1,1,0),c(1,0,1))
out <- class2z(cvec,id,k,dmat)


# Data:      
#   id 2 2 2 4 4 4 4 4 6  6  6  6
#class 1 1 2 1 1 2 2 2 1  1  2  2


out$z[1:6,]
#     [,1] [,2] [,3] 
#[1,]    1    0    0
#[2,]    1    1    0
#[3,]    1    0    0
#[4,]    1    1    0
#[5,]    1    1    0
#[6,]    1    1    0
out$zid[1:6]
# 21 22 23 24 25 26 
#  2  2  2  2  2  2
out$dmat
#      a1 a2 a3 
#(1,1)  1  0  0
#(1,2)  1  1  0
#(2,2)  1  0  1
#

}
\keyword{models}
% Converted by Sd2Rd version 1.21.
