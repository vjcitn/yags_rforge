print.alr <- function(x, digits = NULL, quote = FALSE, prefix = "")
{
# print function for alternating logistic regression, Splus @(#) print.alr.q 3.1 97/03/23
# /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.print.alr.q
	if(is.null(digits)) digits <- options()$digits else options(digits = 
			digits)
# @(#) chanlib revisions: print.alr.q 3.1 97/03/23
	cat("\n")
	cat(x$title)
	cat("\n")
	cat(x$version, "\n")
	cat("\nCall:\n")
	dput(x$call)	#
#	cat("\nFormula:\n")
#	print(x$formula)
#	cat("\nContrasts:\n")
#	print(x$contrasts)
#	ys <- matrix(rep(as.matrix(x$id, ncol = 1), 5), ncol = 5)
#	ys[, 2] <- x$y
#	ys[, 3] <- x$linear.predictors
#	ys[, 4] <- x$fitted.values
#	ys[, 5] <- x$residuals
#	dimnames(ys) <- list(1:length(x$y), c("ID", "Y", "LP", "fitted", 
#		"Residual"))
#	cat("\nFitted Values:\n")
#	print(ys, digits = digits)
	nas <- x$nas
	if(any(nas))
		cat("\n\nCoefficients: (", sum(nas), 
			" not defined because of singularities)\n", sep = "")
	else cat("\n\nCoefficients:\n")
	print(x$coefficients, digits = digits)
	cat("\nAlpha:\n")
	print(x$alpha, digits = digits)	#	p <- length(x$coefficients)
#	q <- length(x$alpha)
#	if(!is.null(x$variance)) {
#		cat("\n\nVariance/Covariance of Estimates:\n\n")
#		print(x$variance)
#	}
#
	cat("\nNumber of observations : ", x$nobs)
	cat("\nNumber of Iterations   : ", x$iterations, "\n")
	invisible(x)
}
