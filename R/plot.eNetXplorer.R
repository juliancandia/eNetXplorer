plot.eNetXplorer <- function(x, plot.type=c("summary","lambdaVsQF","measuredVsOOB","contingency",
"featureCaterpillar","featureHeatmap"), alpha.index=NULL, stat=c("freq","coef"), ...)
{
    plot.type = match.arg(plot.type)
    
    if (plot.type=="summary") { # QF and model.vs.null.pval vs alpha
        plotSummary(x=x, ...)
    }
    
    if (plot.type=="lambdaVsQF") { # QF distribution vs lambda
        plotLambdaVsQF(x=x, alpha.index=alpha.index, ...)
    }
    
    if (plot.type=="measuredVsOOB") { # true vs predicted response
        plotMeasuredVsOOB(x=x, alpha.index=alpha.index, ...)
    }
    
    if (plot.type=="contingency") { # contingency matrix for categorical models
        if (x$family%in%c("binomial","multinomial")) {
            plotContingency(x=x, alpha.index=alpha.index, ...)
        } else {
            warning("Contingency plots only applicable to categorical models")
        }
    }
    
    if (plot.type=="featureHeatmap") { # feature heatmaps
        stat = match.arg(stat)
        plotFeatureHeatmap(x=x, alpha.index=alpha.index, stat=stat, ...)
    }
    
    if (plot.type=="featureCaterpillar") { # feature caterpillar plots
        stat = match.arg(stat)
        plotFeatureCaterpillar(x=x, alpha.index=alpha.index, stat=stat,...)
    }
}
