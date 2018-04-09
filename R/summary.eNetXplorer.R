summary.eNetXplorer <- function(x, ...)
{
    coeff_mat = cbind(x$alpha,x$best_lambda,x$model_QF_est,x$QF_model_vs_null_pval)
    colnames(coeff_mat) = c("alpha","lambda.max","QF.est","model.vs.null.pval")
    res <- list(call=x$call,
    coefficients=coeff_mat
    )
    class(res) <- "summary.eNetXplorer"
    res
}
