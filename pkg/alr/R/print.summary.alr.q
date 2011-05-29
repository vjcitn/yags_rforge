print.summary.alr <- function(x, digits = NULL, quote = FALSE, prefix = "")
{
# print function for alternating logistic regression, Splus @(#) print.summary.alr.q 3.1 97/03/23
# /usr16/stdevs/stdev0f/SLIBS/alr.dev/SCCS/s.print.summary.alr.q
#
	if(is.null(digits)) digits <- options()$digits else options(digits = 
			digits)
	cat("\n")
	cat(x$title)
	cat("\n")
	cat(x$version, "\n")
	cat("\nCall:\n")
	dput(x$call)
	cat("\nSummary of Residuals:\n")
	print(x$residual.summary, digits = digits)
	nas <- x$nas[1]
	if(any(nas))
		cat("\n\nCoefficients: (", sum(nas), 
			" not defined because of singularities)\n", sep = "")
	else cat("\n\nCoefficients:\n")
	print(x$coefficients, digits = digits)
	nas <- x$nas[2]
	if(any(nas))
		cat("\n\nAlpha: (", sum(nas), 
			" not defined because of singularities)\n", sep = "")
	else cat("\n\nAlpha:\n")
	print(x$alpha, digits = digits)
	cat("\n\nNumber of observations : ", x$nobs)
	cat("\nNumber of Iterations   : ", x$iterations, "\n")
	invisible(x)
}
