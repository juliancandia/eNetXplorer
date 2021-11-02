plotMeasuredVsOOB <- function (
x, alpha.index=NULL,
xlab=NULL, ylab=NULL, cex.lab=0.95, main=NULL, col.main = "black", cex.main=0.85,
instance.label=TRUE,
instance.label.cex=NULL, instance.label.offset=NULL,
instance.label.added.margin=NULL, col=NULL, # numerical only
box.wex=NULL, box.range=NULL, box.col=NULL,
transparency=NULL, jitter=NULL, cex.pt=NULL, class.color=NULL, # categorical only
...)
{
    if (is.null(alpha.index)) {
        alpha.index = 1:length(x$alpha)
    }
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
        plotMeasuredVsOOBNumer(x=x, alpha.index=alpha.index, xlab=xlab, ylab=ylab, cex.lab=cex.lab, main=main, col.main=col.main, cex.main=cex.main, instance.label=instance.label, instance.label.cex=instance.label.cex, instance.label.offset=instance.label.offset, instance.label.added.margin=instance.label.added.margin, col=col,
         ...)
    } else if (x$family%in%c("binomial","multinomial")) {
        if (is.null(instance.label.cex)) {
            instance.label.cex = 0.45
        }
        if (is.null(instance.label.offset)) {
            instance.label.offset = 0.6
        }
        if (is.null(box.wex)) {
            box.wex = 0.5
        }
        if (is.null(box.range)) {
            box.range = 0
        }
        if (is.null(box.col)) {
            box.col = "white"
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
        plotMeasuredVsOOBCateg(x=x, alpha.index=alpha.index, xlab=xlab, ylab=ylab, cex.lab=cex.lab, main=main, col.main=col.main, cex.main=cex.main, instance.label=instance.label, instance.label.cex=instance.label.cex, instance.label.offset=instance.label.offset,
            box.wex=box.wex, box.range=box.range, box.col=box.col,
            transparency=transparency, jitter=jitter, cex.pt=cex.pt, class.color=class.color, ...)
    }
}

