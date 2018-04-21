summaryPDF <- function (x, path=getwd(), filename="eNetXplorerSummary.pdf")
{
    pdf(file.path(path,filename),width=7,height=5)
    plot(x=x,plot.type="summary")
    for (i_alpha in 1:length(x$alpha)) {
        plot(x=x,plot.type="lambdaVsQF",alpha.index=i_alpha)
        plot(x=x,plot.type="measuredVsOOB",alpha.index=i_alpha)
        if (x$family%in%c("binomial","multinomial")) {
            plot(x=x,plot.type="contingency",alpha.index=i_alpha)
        }
        for (stat in c("coef","freq")) {
            plot(x=x,plot.type="featureCaterpillar",alpha.index=i_alpha,stat=stat)
            plot(x=x,plot.type="featureHeatmap",alpha.index=i_alpha,stat=stat)
        }
    }
    dev.off()
}
