plot.lambda.vs.QF.eNetXplorer <- function (
x, alpha.index, xlab="lambda", ylab="QF (response vs out-of-bag predicted)", cex.lab=0.95,
main=NULL, col.main="black", cex.main=0.95, log="x", type="b",
...)
{
    for (i_alpha in alpha.index) {
        plot(x$lambda_values[[i_alpha]],x$lambda_QF_est[[i_alpha]],
        log=log,type=type,xlab=xlab,ylab=ylab,cex.lab=cex.lab, ...)
        if (is.null(main)) {
            main.title = paste0("alpha=",x$alpha[i_alpha])
            if (x$QF_label!="QF") {
                main.title = paste0(main.title," ; QF=",x$QF_label)
            }
        } else {
            main.title = main
        }
        title(main=main.title, col.main=col.main, cex.main=cex.main)
        abline(v=x$best_lambda[i_alpha],lty=3)
    }
}
