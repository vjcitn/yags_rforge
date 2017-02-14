#' plot a qqplot relative to a weibull distribution
#' @import survival
#' @importFrom graphics abline
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom stats pweibull
#' @param x vector of numbers
#' @param conf.int numeric confidence coefficient for EDF CIs
#' @param \dots not used
#' @details The weibull parameters are estimated using survival::survreg
#' @export
weiqqci <- function(x, conf.int=0.95, ...)
{
# x = any vector of numbers
# plots a qqplot relative to a weibull distribution
# enhanced by CI for the EDF at sample points
        FF <- weifit2(x)
        EE <- nedf(x, conf.int=conf.int, ...)
        xseq <- EE$time
        msurv <- 1 - pweibull(xseq, shape = FF[1], scale = FF[2])
        esurv <- EE$surv
        plot(esurv, msurv, xlab = "EDF", ylab = "Weibull Model")
        segments(EE$lower, msurv, EE$upper, msurv)
        abline(0, 1, lty = 2)
}

weifit2 <- function(x)
{
# extracts necessary parm ests from Weibull fit
        fit <- survreg(Surv(x, rep(1, length(x))) ~ 1)
        #sha <- 1/exp(summary(fit)$scale)
        sha <- 1/(summary(fit)$scale)
        scale <- exp(fit$coef)
        c(shape = sha, scale = scale)
}

nedf <- function(x, conf.int=0.95, ...)
survfit(Surv(x, rep(1, length(x))) ~ 1, conf.int=conf.int)

