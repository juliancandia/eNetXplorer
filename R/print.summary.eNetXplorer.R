print.summary.eNetXplorer <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}
