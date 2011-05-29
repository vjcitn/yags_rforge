\name{alr}
\alias{alr}
\title{
Alternating Logistic Regression for multivariate
binary data @(#) chanlib revisions: alr.d alr.d 4.2
}
\usage{
alr(x, y, id, z=0, zmast=0, zid=0, zlocs=0, binit, ainit, tol=0.001,
       bweight="full", depmodel="general", subclust=0)
}
\arguments{
\item{x}{
matrix of predictors of logit probability of marginal outcome.
should include a column of 1 for estimation of "baseline"
probability if desired: Nxp, p includes intercept
}
\item{y}{
vector of binary outcomes corresponding to rows of x: Nx1
}
\item{id}{
cluster discriminator: N-vector
}
\item{z}{
matrix of predictors for pairwise log odds ratio regression.
this may be omitted for certain choices of "depmodel",
see below.  if used,
this matrix may take one of two forms.  it may specify
directly, for each cluster, the form of the dependency
design.  in this case, it has dimension (sum_i(n_i*(n_i - 1))) x q.
on the other hand, if the data to be analyzed are replicated
(possibly subject to missingness), then if n is the size of a complete
cluster, z has dimension (n*(n-1)) x q.  see zmast and
zlocs below
}
\item{zmast}{
if 1, then z is a "master" design for pairwise log oddsratio
regression.  each cluster is a replicate of this master
design, but clusters may have missing elements, thus
missing pairs.  the missing information is extracted from
zlocs, below.  if 0, then z has dimension (sum_i(n_i*(n_i-1))) x q
}
\item{zid}{
cluster discriminators for the z matrix, used only if zmast == 0.
if used, has dimension (sum_i(n_i *(n_i-1))) x 1.
}
\item{zlocs}{
used only if zmast == 1.  this is an N-vector of "locations" in
the dependency design.  if replication is perfect and complete,
and n is the cluster size, and there are C clusters, then
this vector should have the value rep(1:n,C). 
}
\item{binit}{
p-vector of initial values for beta -- required; should
be obtained from a logistic regression fit assuming independence
}
\item{ainit}{
q-vector of initial values for alpha -- difficult to 
recommend starting values here, but rep(.01,q) often
works
}
\item{tol}{
maximum relative change in parameter tolerated before
asserting convergence
}
\item{bweight}{
currently takes one of two values: "full" or "independence".
if "independence", then the estimating equation for
estimating regression parameters beta uses a diagonal
weighting matrix, and the estimated beta should be
identical to the GLIM estimate.  if "full", then 
the estimating equation for beta uses the estimated
covariance (nxn) matrix of the outcomes as weight.
the "independence" setting may be useful for very
large clusters (n>50?) with sparse outcomes, in which the
weighting matrix can tend to singularity.
}
\item{depmodel}{
currently takes one of three values: "general", "exchangeable",
"1-nested".  The "general" setting uses a Z matrix (see above)
to specify the dependency model.  The "exchangeable" setting
uses no Z matrix, but assumes a common pairwise log-odds ratio
for all cluster elements.   The "1-nested" setting estimates
a log-odds ratio regression with first parameter the log
pairwise odds ratio for any two elements in a cluster,
the second parameter is the incremental dependency between
members of the same subcluster.  subclusters are identified
through the "subclust" vector; see below.
}
\item{subclust}{
N-vector discriminating subclusters within clusters.  Thus
two 8-clusters might have id=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
subclust=c(10,10,11,11,11,11,11,11,20,20,20,20,21,21,21,21).
the first cluster falls into 2 subclusters, the first is size
2, the second size 6; the second has 2 subclusters each of
size 4.  the 1-nested model says that the log-pairwise
odds ratios between elements of the same cluster is a0+a1
if the two elements are ALSO of the same subcluster, and
a0 if they are from different subclusters.
}}
\value{
a structure including estimates of beta, alpha, and the
parameter covariance matrix
}
\note{
The most complicated idea in using this program concerns
the structure and modelling of the pairwise odds-ratio
data.  We deal with the "general" case
first.
If n (clustersize) = 4, then there are 6 = (n choose 2)
possible pairwise odds ratios: O12, O13, O14, O23, O24, O34.
However, to estimate these, we use the 12 (n*n-1) redundant 
pairings (Y1,Y2), (Y1,Y3), (Y1,Y4), (Y2,Y1), (Y2,Y3), \dots, 
(Y4,Y3).  It is in this order that we index the covariate vectors
for regression analysis : Zij, ordered (i = 1:n, j = 1:n, j!=i ).
The zloc argument tells uses the indices 1:n to tell which
pairing element each outcome corresponds to.  Thus suppose
the design is a replication with n = 4.  Complete data has
zloc = c(1,2,3,4); if a cluster is missing, say, the outcome
for the "second" subject, zloc = c(1,3,4), indicating that
the cluster contributes to the estimation only of O13, O14, O34.
Of course, the odds-ratio parameterization may be more
parsimonious; if the model O12=O13=\dots=O34 is adopted, the specific
values of zloc are not relevant to estimation of the common
pairwise log odds ratio.  But in the event that q>1, zloc
tells which subjects are to be used for estimation of which
combination of odds ratio parameters.


The simplest model to fit with alr is the "exchangeable" model,
in that no information is needed other than x,y, and id.


The "1-nested" model requires a subclust vector in addition
to x, y, and id.
}
\examples{
data(alrset)
a1 <- alr(alr.y ~ alr.x - 1, id=alr.id, depm="exchangeable", ainit=0.01)
summary(a1)

#using a master z matrix for a balanced design
ZMAST <- rep(1,12)
ZMAST <- cbind(ZMAST,c(0,0,0,0,1,1,0,1,1,0,1,1))
ZIND <- rep(1:4,125)
Y <- as.vector(alr.y)
X <- as.vector(alr.x[,2])
NY <- split(alr.y,alr.id)
NY <- unlist(lapply(NY,function(x)rev(x)))
NX <- split(X,alr.id)
NX <- unlist(lapply(NX,function(x)rev(x)))
Y <- as.vector(NY)
X <- as.vector(NX)
mast.out <- alr(Y ~ X, id=alr.id, depmod="general", z=ZMAST,
zmast=1, zloc=ZIND,
               ainit = c(0.01,0.01) )
summary(mast.out)


}
\keyword{cluster}
%/usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.alr.d
%alr.d 4.2 98/03/10
% Converted by Sd2Rd version 1.21.
