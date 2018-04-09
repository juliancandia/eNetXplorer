plot.eNetXplorer <- function(x, plot.type, alpha.index=NULL, stat=NULL, ...)
#plot.type=c("lambda.vs.QF","measured.vs.oob","feature.heatmap","feature.caterpillar")
#stat=c("freq","coef") if plot.type in c("feature.heatmap","feature.caterpillar")
{
    if (is.null(plot.type)|!plot.type%in%c("summary","lambda.vs.QF","measured.vs.oob","contingency",
    "feature.heatmap","feature.caterpillar")) {
        stop("parameter \'plot.type\' must be in c(\"summary\",\"lambda.vs.QF\",\"measured.vs.oob\"
        ,\"contingency\",\"feature.heatmap\",\"feature.caterpillar\")\n")
    }
    n_alpha = length(x$alpha)
    if (!is.null(alpha.index)) {
        if (!(is.numeric(alpha.index)&&all(floor(alpha.index)==alpha.index,na.rm=T) # check for integer
        &&min(alpha.index)>0&&max(alpha.index)<=n_alpha)) {
            stop("parameter \'alpha.index\' must be an integer in the range [1,length(alpha)]\n")
        }
    } else {
        alpha.index = 1:n_alpha
    }
    
    if (plot.type=="summary") { # QF and model.vs.null.pval vs alpha
        plot.summary.eNetXplorer(x=x, ...)
    }
    
    if (plot.type=="lambda.vs.QF") { # QF distribution vs lambda
        plot.lambda.vs.QF.eNetXplorer(x=x, alpha.index=alpha.index, ...)
    }
    
    if (plot.type=="measured.vs.oob") { # true vs predicted response
        plot.measured.vs.oob.eNetXplorer(x=x, alpha.index=alpha.index, ...)
    }
    
    if (plot.type=="contingency") { # contingency matrix for categorical models
        if (x$family%in%c("binomial","multinomial")) {
            plot.contingency.eNetXplorer(x=x, alpha.index=alpha.index, ...)
        } else {
            warning("Contingency plots only applicable to categorical models")
        }
    }
    
    if (plot.type=="feature.heatmap") { # feature heatmaps
        if (!is.null(stat)&&stat%in%c("freq","coef")) {
            plot.feature.heatmap.eNetXplorer (x=x, alpha.index=alpha.index, stat=stat, ...)
        } else {
            stop("parameter \'stat\' must be in c(\"freq\",\"coef\")\n")
        }
    }
    
    if (plot.type=="feature.caterpillar") { # feature caterpillar plots
        if (!is.null(stat)&&stat%in%c("freq","coef")) {
            plot.feature.caterpillar.eNetXplorer(x=x, alpha.index=alpha.index, stat=stat,...)
        } else {
            stop("parameter \'stat\' must be in c(\"freq\",\"coef\")")
        }
    }
}
