summary.alr <- function(object)
{
# summary function for alternating logistic regression, Splus @(#) summary.alr.q 3.1 97/03/23
# /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.summary.alr.q
#
	coef <- object$coefficients
	alpha <- object$alpha
	resid <- object$residuals
	p <- length(coef)
	q <- length(alpha)
	n <- length(resid)
	nasb <- is.na(coef)
	nasa <- is.na(alpha)
	cnames <- names(coef[!nasb])
	clnames <- names(alpha[!nasa])
	coef <- matrix(rep(coef[!nasb], 3), ncol = 3)
	alpha <- matrix(rep(alpha[!nasa], 3), ncol = 3)
	dimnames(coef) <- list(cnames, c("Estimate", "Robust S.E.", "Robust z")
		)
	dimnames(alpha) <- list(clnames, c("Estimate", "Robust S.E.", 
		"Robust z"))
	sebet <- sqrt(diag(object$variance[1:p, 1:p, drop = FALSE]))
	sealp <- sqrt(diag(object$variance[(p + 1):(p + q), (p + 1):(p + q), 
		drop = FALSE]))
	coef[, 2] <- sebet
	coef[, 3] <- coef[, 1]/coef[, 2]
	alpha[, 2] <- sealp
	alpha[, 3] <- alpha[, 1]/alpha[, 2]
	summary <- list()
	summary$call <- object$call
	summary$version <- object$version
	summary$nobs <- object$nobs
	summary$residual.summary <- quantile(as.vector(object$residuals))
	names(summary$residual.summary) <- c("Min", "1Q", "Median", "3Q", "Max"
		)
	summary$title <- object$title
	summary$coefficients <- coef
	summary$alpha <- alpha
	summary$iterations <- object$iterations
	summary$nas <- c(nasb, nasa)
	attr(summary, "class") <- "summary.alr"
	summary
}
