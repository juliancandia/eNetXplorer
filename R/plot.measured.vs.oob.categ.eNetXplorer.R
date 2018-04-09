plot.measured.vs.oob.categ.eNetXplorer <- function (
x, alpha.index, xlab, ylab, cex.lab,
main, col.main, cex.main, instance.label,
instance.label.cex, instance.label.offset,
transparency, jitter, cex.pt, class.color,
...)
{
    n_instance = length(x$instance)
    if (is.null(xlab)) {
        xlab="response"
    }
    if (is.null(ylab)) {
        ylab="out-of-bag predicted (accuracy)"
    }
    for (i_alpha in alpha.index) {
        accuracy = rep(NA,n_instance)
        for (i_instance in 1:n_instance) {
            accuracy[i_instance] = x$predicted_values[[i_alpha]][i_instance,as.numeric(x$response)[i_instance]]
        }
        
        boxplot(as.formula("accuracy ~ x.response"),data=data.frame(x$response,accuracy),
        xlab=xlab,ylab=ylab,cex.lab=cex.lab,outline=F,boxwex=0.5,range=0,...)
        if (is.null(main)) {
            main.title = paste0("alpha=",x$alpha[i_alpha],
            "; lambda=",signif(x$best_lambda[i_alpha],digits=3),
            "; QF=",signif(x$model_QF_est[i_alpha],digits=2))
        } else {
            main.title = main
        }
        title(main=main.title, col.main=col.main, cex.main=cex.main)
        
        class = levels(x$response)
        n_class = length(class)
        if (is.null(class.color)) {
            class.color = 1:n_class
        }
        for (i_class in 1:n_class) {
            instance_select = x$response==class[i_class]
            scatter = runif(sum(instance_select),-jitter,jitter)
            x_plot = i_class+scatter-mean(scatter)
            y_plot = accuracy[instance_select]
            col_par = as.numeric(col2rgb(class.color[i_class]))
            points(x_plot,y_plot,cex=cex.pt,col=rgb(col_par[1],col_par[2],col_par[3],
            transparency,maxColorValue=255),pch=16)
            if (instance.label) {
                textxy(x_plot,y_plot,x$instance[instance_select],
                cex=instance.label.cex,offset=instance.label.offset)
            }
        }
    }
}
