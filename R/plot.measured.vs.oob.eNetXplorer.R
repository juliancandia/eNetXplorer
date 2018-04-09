plot.measured.vs.oob.eNetXplorer <- function (
x, alpha.index=alpha.index,
xlab=NULL, ylab=NULL, cex.lab=0.95, main=NULL, col.main = "black", cex.main=0.85,
instance.label=T,
instance.label.cex=NULL, instance.label.offset=NULL,
instance.label.added.margin=NULL, col=NULL, # numerical only
transparency=NULL, jitter=NULL, cex.pt=NULL, class.color=NULL, # categorical only
...)
{
    if (x$family=="gaussian") {
        if (is.null(instance.label.cex)) {
            instance.label.cex = 0.5
        }
        if (is.null(instance.label.offset)) {
            instance.label.offset = 0.75
        }
        if (is.null(instance.label.added.margin)) {
            instance.label.added.margin = 0.1
        }
        if (is.null(col)) {
            col="red"
        }
        plot.measured.vs.oob.numer.eNetXplorer(x=x, alpha.index=alpha.index, xlab=xlab, ylab=ylab, cex.lab=cex.lab, main=main, col.main=col.main, cex.main=cex.main, instance.label=instance.label, instance.label.cex=instance.label.cex, instance.label.offset=instance.label.offset, instance.label.added.margin=instance.label.added.margin, col=col,
         ...)
    } else if (x$family%in%c("binomial","multinomial")) {
        if (is.null(instance.label.cex)) {
            instance.label.cex = 0.45
        }
        if (is.null(instance.label.offset)) {
            instance.label.offset = 0.6
        }
        if (is.null(transparency)) {
            transparency = 70
        }
        if (is.null(jitter)) {
            jitter = 0.25
        }
        if (is.null(cex.pt)) {
            cex.pt = 1.7
        }
        plot.measured.vs.oob.categ.eNetXplorer(x=x, alpha.index=alpha.index, xlab=xlab, ylab=ylab, cex.lab=cex.lab, main=main, col.main=col.main, cex.main=cex.main, instance.label=instance.label, instance.label.cex=instance.label.cex, instance.label.offset=instance.label.offset,
            transparency=transparency, jitter=jitter, cex.pt=cex.pt, class.color=class.color, ...)
    }
}

