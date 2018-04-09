pdf.summary <- function (x, path=getwd(), filename="eNetXplorer.summary.pdf", ...)
{
    pdf(file.path(path,filename),width=7,height=5)
    plot(x=x,plot.type="summary")
    for (i_alpha in 1:length(x$alpha)) {
        plot(x=x,plot.type="lambda.vs.QF",alpha.index=i_alpha)
        plot(x=x,plot.type="measured.vs.oob",alpha.index=i_alpha)
        if (x$family%in%c("binomial","multinomial")) {
            plot(x=x,plot.type="contingency",alpha.index=i_alpha)
        }
        for (stat in c("coef","freq")) {
            plot(x=x,plot.type="feature.caterpillar",alpha.index=i_alpha,stat=stat)
            plot(x=x,plot.type="feature.heatmap",alpha.index=i_alpha,stat=stat)
        }
    }
    dev.off()
}
