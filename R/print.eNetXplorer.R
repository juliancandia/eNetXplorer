print.eNetXplorer <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    x$alpha_best_lambda
}
