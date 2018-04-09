plot.measured.vs.oob.numer.eNetXplorer <- function (
x, alpha.index, xlab, ylab, cex.lab,
main, col.main, cex.main, instance.label,
instance.label.cex, instance.label.offset,
instance.label.added.margin, col,
...)
{
    if (is.null(xlab)) {
        xlab="response"
    }
    if (is.null(ylab)) {
        ylab="out-of-bag predicted (mean, sd)"
    }
    for (i_alpha in alpha.index) {
        
        pred_mean = x$predicted_values[[i_alpha]][,1]
        pred_sd = x$predicted_values[[i_alpha]][,2]
        
        x_range = range(x$response)
        if (!is.null(instance.label)&&is.logical(instance.label)&&instance.label==T) {
            x_range = x_range + c(-1,1)*(x_range[2]-x_range[1])*instance.label.added.margin
        }
        y_range = range(c(pred_mean-pred_sd,pred_mean+pred_sd),na.rm=T)
        
        plot(x_range,y_range,type="n",xlab=xlab,ylab=ylab,cex.lab=cex.lab,...)
        if (is.null(main)) {
            main.title = paste0("alpha=",x$alpha[i_alpha],
            "; lambda=",signif(x$best_lambda[i_alpha],digits=3),
            "; QF=",signif(x$model_QF_est[i_alpha],digits=2))
        } else {
            main.title = main
        }
        title(main=main.title, col.main=col.main, cex.main=cex.main)
        
        df = data.frame(x=x$response,y=pred_mean)
        mod <- lm(y ~ x, data = df)
        newx_range = x_range + c(-1,1)*(x_range[2]-x_range[1])*0.3
        newx <- seq(newx_range[1], newx_range[2], length.out=200)
        preds <- predict(mod, newdata = data.frame(x=newx),
        interval = 'confidence')
        polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'lightblue1', border = NA)
        abline(mod,lty=1,col="blue") # model
        lines(newx, preds[,3], lty=3, col="blue") # intervals
        lines(newx, preds[,2], lty=3, col="blue")
        
        points(x$response,pred_mean,pch=16,cex=1,col=col,lty=3,type="p")
        arrows(x$response,pred_mean-pred_sd,x$response,pred_mean+pred_sd,col=col,length=0.025,angle=90,code=3)
        if (instance.label) {
            textxy(x$response,pred_mean,x$instance,cex=instance.label.cex,offset=instance.label.offset)
        }
    }
}
