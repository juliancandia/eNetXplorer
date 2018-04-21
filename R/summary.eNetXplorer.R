summary.eNetXplorer <- function(object, ...)
{
    coeff_mat = cbind(object$alpha,object$best_lambda,object$model_QF_est,object$QF_model_vs_null_pval)
    rownames(coeff_mat) = rep("",nrow(coeff_mat))
    colnames(coeff_mat) = c("alpha","lambda.max","QF.est","model.vs.null.pval")
    res <- list(call=object$call,
    coefficients=coeff_mat
    )
    class(res) <- "summary.eNetXplorer"
    res
}
