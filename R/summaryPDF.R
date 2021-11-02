summaryPDF <- function (x, dest_dir=getwd(), dest_dir_create=TRUE, dest_file="eNetSummary.pdf")
{
    if ((!dir.exists(dest_dir))&&(dest_dir_create)) {
        dir.create(dest_dir)
    }
    pdf(file.path(dest_dir,dest_file),width=7,height=5)
    plot(x=x,plot.type="summary")
    for (i_alpha in 1:length(x$alpha)) {
        plot(x=x,plot.type="lambdaVsQF",alpha.index=i_alpha)
        if (x$family%in%c("gaussian","binomial","multinomial")) {
            plot(x=x,plot.type="measuredVsOOB",alpha.index=i_alpha)
        }
        if (x$family%in%c("binomial","multinomial")) {
            plot(x=x,plot.type="contingency",alpha.index=i_alpha)
        }
        if (x$family=="cox") {
            plot(x=x,plot.type="KaplanMeier",alpha.index=i_alpha)
            if (x$survAUC) {
                plot(x=x,plot.type="survROC",alpha.index=i_alpha,survAUC_time=x$survAUC_time)
            }
        }
        for (stat in c("coef","freq")) {
            plot(x=x,plot.type="featureCaterpillar",alpha.index=i_alpha,stat=stat)
            plot(x=x,plot.type="featureHeatmap",alpha.index=i_alpha,stat=stat)
        }
    }
    dev.off()
}
