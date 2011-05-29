alr <- function(formula = formula(data), id = id, weights = NULL, data = sys.parent(), 
	subset, na.action, contrasts = NULL, z = 0, zmast = 0, zid = 0, zlocs
	 = 0, binit = NULL, ainit, bweight = "full", depmodel = "general", 
	subclust = 0, clnames = NULL, tol = 0.001, maxiter = 25, silent = TRUE)
{
#/* revision for symmetry and weighting: */
#/* /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.alr.q */
#/* 4.4 98/02/24 */

# genesis alternating logistic regression, Splus @(#) alr.q 3.2 97/03/23
       print("alternating logistic regression - Splus, @(#) alr.q 4.4 98/02/24"
		)	#@(#) chanlib revisions: alr.q 3.2 97/03/23
#
# genesis: /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.alr.q
#
#  Create model frame m
#  --------------------
	call <- match.call()
	m <- match.call(expand = FALSE)
	m$binit <- m$ainit <- m$tol <- m$bweight <- m$depmodel <- m$clnames <- 
		m$subclust <- m$contrasts <- m$silent <- m$z <- m$zmast <- m$
		zlocs <- m$zid <- m$maxiter <- NULL
	if(is.null(m$id)) {
		m$id <- as.name("id")
	}
	if(is.null(m$na.action)) {
	}
	else {
		if(m$na.action != "na.omit") {
			print("Warning: Only na.omit is implemented for alr")
			print("         continuing with na.action=na.omit")
			m$na.action <- as.name("na.omit")
		}
	}
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())	#
#
#  Extract x,y etc from model frame
#  --------------------------------
	Terms <- attr(m, "terms")
	y <- as.matrix(model.extract(m, response))
	x <- model.matrix(Terms, m, contrasts)
	id <- model.extract(m, id)
	if(is.null(id)) {
		stop("Id variable not found")
	}
	nobs <- nrow(x)
	nc <- ncol(x)
	xnames <- dimnames(x)[[2]]
	if(is.null(xnames)) {
		xnames <- paste("x", 1:nc, sep = "")
		dimnames(x) <- list(NULL, xnames)
	}
#
#
#  Deal with binit 
#  ---------------
	family <- binomial
	if(!is.null(binit)) {
		binit <- as.double(binit)
		if(nrow(binit) != nc) {
			stop("Size of binit is not the same as number of independents"
				)
		}
	}
	else {
		print("Running glm to get initial estimates")
		binit <- as.numeric(glm(m, family = family)$coef)
		print(binit)
	}
	p <- length(binit)	#
#
#  .... and ainit
#  --------------
	if(!is.null(ainit)) {
		ainit <- matrix(as.double(ainit), ncol = 1)
	}
	else {
		print("ainit must be specified")
	}
	q <- nrow(ainit)
	if(is.null(clnames))
		clnames <- paste("a", 1:q, sep = "")
	cov <- matrix(0, nrow = p + q, ncol = p + q)
	if(!is.matrix(z))
		z <- as.matrix(z)
	nzrows <- nrow(z)
	if(length(subclust) == nobs) subclust <- as.double(subclust)	#
#
#  pick weighting for beta
#  -----------------------
	b.indep <- 0
	if(bweight == "independence")
		b.indep <- 1
	else {
		if(bweight != "full")
			stop(paste("unknown bweight option:", bweight))	#
	}
	modeltypes <- c("general", "exchangeable", "1-nested")
	modty <- match(depmodel, modeltypes, nomatch = NA)
	if(is.na(modty))
		stop(paste("model type must be one of (", modeltypes, ")"))
	exch <- as.integer(modty == 2)
	is.twocl <- as.integer(modty == 3)
	print("nobs")
	print(nobs)	#
	if(is.null(weights))
		weights <- rep(1, nobs)
	else if(length(weights) != nobs)
		stop("weight vector wrong length")
	else weights <- 1/sqrt(weights)
	out <- .C("alr",
		as.double(x),
		as.double(y),
		as.double(z),
		as.integer(zmast),
		as.double(id),
		as.double(zid),
		as.double(zlocs),
		as.integer(nobs),
		as.integer(nzrows),
		as.integer(p),
		as.integer(q),
		beta = as.double(binit),
		alpha = as.double(ainit),
		cov = as.double(cov),
		as.double(tol),
		iter = as.integer(maxiter),
		as.integer(exch),
		as.integer(is.twocl),
		as.integer(b.indep),
		as.double(subclust),
		as.integer(silent),
		as.double(weights), PACKAGE="alr")	#
#
#  Put results together
#  --------------------
	fit <- list()
	attr(fit, "class") <- c("alr", "glm")
	fit$title <- "ALR:  ALTERNATING LOGISTIC REGRESSION"
	fit$version <- "alr S-function, version 4.4 98/02/24"
	fit$call <- call
	fit$terms <- Terms
	fit$formula <- as.vector(attr(Terms, "formula"))
	fit$contrasts <- attr(x, "contrasts")
	fit$nobs <- nobs
	fit$iterations <- out$iter
	fit$coefficients <- as.vector(out$beta)
	fit$alpha <- as.vector(out$alpha)
	fit$variance <- matrix(out$cov, nrow = p + q)
	dimnames(fit$variance) <- list(c(xnames, clnames), c(xnames, clnames))
	fit$nas <- is.na(fit$coefficients)
	names(fit$coefficients) <- xnames
	names(fit$alpha) <- clnames
	eta <- as.vector(x %*% fit$coefficients)
	fit$linear.predictors <- eta
	mu <- as.vector(family()$linkinv(eta))
	fit$fitted.values <- mu
	fit$residuals <- y - mu
	fit$family <- family
	fit$y <- as.vector(y)
	fit$id <- as.vector(id)
	fit
}


class2z <- function(cvec, id, k, dmat)
{
# (v)class2z : adaptation of Pat Heagerty's class2z function
# to produce "permutation invariant" (redundant) Z matrix for alr
# vclass2z 4.4 98/02/24
# based on Heagerty: c2z 1.3 96/07/23
	if(max(cvec) != k) {
		print(paste("k = ", k, " is not equal to the cvec max."))
		print(paste("cvec max = ", max(cvec)))
		stop("exit class2z()")
	}
	if(min(cvec) != 1) {
		print(paste("cvec min = ", min(cvec), ", must be 1."))
		stop("exit class2z()")
	}
	if(length(cvec) != length(id)) {
		print("cvec and id are unequal length.")
		stop("exit class2z()")
	}
	if(nrow(dmat) != (k * (k + 1))/2) {
		print("dmat must have k*(k+1)/2 rows.")
		stop("exit class2z()")
	}
#
	x <- table(id)
	nclust <- length(x)
	n <- length(cvec)
	total <- 0
	x <- table(table(id))
	y <- as.integer(names(x))
	for(i in 1:length(x)) {
		if(y[i] > 1) {
			total <- total + x[i] * ((y[i] * (y[i] - 1))/2)
		}
		else {
			total <- total + x[i]
		}
	}
	z <- matrix(0, total, ncol(dmat))
	zid <- matrix(0, total, 1)
	flag <- 0
	out <- .C("class2z",
		class = as.double(cvec),
		id = as.double(id),
		dmat,
		z = z,
		zid = as.double(zid),
		k = as.integer(k),
		n = as.integer(n),
		nrz = as.integer(total),
		p = as.integer(nrow(dmat)),
		q = as.integer(ncol(dmat)),
		nclust = as.integer(nclust),
		flag = as.integer(flag), PACKAGE="alr")
	dmat.names <- rep("", (k * (k + 1))/2)
	count <- 0
	for(i in 1:k) {
		for(j in i:k) {
			count <- count + 1
			dmat.names[count] <- paste("(", i, ",", j, ")", sep = 
				"")
		}
	}
	dimnames(dmat) <- list(dmat.names, paste("a", c(1:ncol(dmat)), sep = ""
		))
	x <- list(z = out$z, zid = out$zid, model = dmat, flag = out$flag)	
	#/*    61  void refab(  z, zid, nrz, q , newzp)  */
	AAA <- .C("refab",
		as.double(out$z),
		as.double(out$zid),
		as.integer(total),
		as.integer(ncol(out$z)),
		ans = double(2 * total * ncol(out$z)), PACKAGE="alr")
	OLDZID <- out$zid
	soz <- split(OLDZID, OLDZID)
	newzid <- unlist(lapply(soz, function(x)
	rep(x, 2)))
	list(z = matrix(AAA$ans, nc = ncol(out$z)), zid = newzid, dmat=dmat, flag = 
		out$flag)
}
