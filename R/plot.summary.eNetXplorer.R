plot.summary.eNetXplorer <- function (
x, show.pval.ref=T, main=NULL, col.main="black", cex.main=0.95, line=1,
...)
{
    y1 = x$model_QF_est
    y2 = -log10(x$QF_model_vs_null_pval)
    
    opar <- par() # to make a copy of current settings
    par(mar = c(5,4,4,4)) # default is mar = c(lower=5,left=4,top=4,right=2) + 0.1.
    
    plot(x$alpha,y1,type="b",yaxt="n",col="red",xlab="alpha",ylab="")
    axis(side=2, las=3, cex.axis=0.85, col.axis="red")
    mtext(text="QF (response vs out-of-bag predicted)", side=2, col="red", line=2.1)
    
    if (max(y2)==min(y2)) {
        delta = 1
        baseline = min(y2)-0.5
    } else {
        delta = max(y2)-min(y2)
        baseline = min(y2)
    }
    rescale = (max(y1)-min(y1))/delta
    shift = min(y1)-baseline*rescale
    y2n = y2*rescale+shift
    
    lines(x$alpha,y2n,type="b",col="blue")
    axis(side=4, las=3, labels=signif((axTicks(2)-shift)/rescale,digits=3), at=axTicks(2),cex.axis=0.85, col.axis="blue")
    mtext(text="-log10 pval (model vs null)", side=4, col="blue", line=2.1)
    
    if (show.pval.ref==T) {
        pval_ref = c(0.001,0.01,0.05,0.1)
        pval_ref_transf = -log10(pval_ref)*rescale+shift
        par(mgp = c(0,-2.8,0)) # sets location of axis labels, tick mark labels, and tick marks.
        for (i in 1:length(pval_ref)) {
            abline(h=pval_ref_transf[i],lty=3,col="blue")
            axis(side=4, las=1, labels=paste0("pval=",pval_ref), at=pval_ref_transf+0.01*(max(y1)-min(y1)),cex.axis=0.6, col.axis="blue", tcl=0.) # tcl to set tick length and orientation
        }
        par(mgp=c(3,1,0)) # resets defaults
    }

    if (is.null(main)) {
        main = "model performance vs alpha values"
    }
    title(main, col.main=col.main, line=line, cex.main=cex.main)
    
    par(opar) # restores par settings
}
