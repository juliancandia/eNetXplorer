plotSummary <- function (
x, show.pval.ref=TRUE, main=NULL, col.main="black", cex.main=0.95, line=1,
...)
{
    y1 = x$model_QF_est
    y2 = -log10(x$QF_model_vs_null_pval)
    
    par_mar = par()$mar # to save current settings
    par(mar = c(5,4,4,4))
    
    plot(x$alpha,y1,type="b",yaxt="n",col="red",xlab="alpha",ylab="",...)
    if (max(y1)>min(y1)) {
        axis(side=2, las=3, cex.axis=0.85, col.axis="red")
    } else {
        axis(side=2, las=3, cex.axis=0.85, col.axis="red", labels=signif(max(y1),digits=2), at=max(y1))
    }
    mtext(text="QF (response vs out-of-bag predicted)", side=2, col="red", line=2.1)
    
    if (max(y2,na.rm=T)>min(y2,na.rm=T)) {
        delta1 = max(y1)-min(y1)
        delta2 = max(y2,na.rm=T)-min(y2,na.rm=T)
        baseline2 = min(y2,na.rm=T)
        rescale = delta1/delta2
        shift = min(y1)-baseline2*rescale
        y2n = y2*rescale+shift
        axis(side=4, las=3, labels=signif((axTicks(2)-shift)/rescale,digits=3), at=axTicks(2),cex.axis=0.85, col.axis="blue")
        if (show.pval.ref==T) {
            pval_ref = c(0.001,0.01,0.05,0.1)
            pval_ref_transf = -log10(pval_ref)*rescale+shift
            par_mgp = par()$mgp # to save current settings
            par(mgp = c(0,-2.8,0)) # sets location of axis labels, tick mark labels, and tick marks.
            for (i in 1:length(pval_ref)) {
                abline(h=pval_ref_transf[i],lty=3,col="blue")
                axis(side=4, las=1, labels=paste0("pval=",pval_ref), at=pval_ref_transf+0.01*(delta1),cex.axis=0.6, col.axis="blue", tcl=0.) # tcl to set tick length and orientation
            }
            par(mgp=par_mgp) # restores par settings
        }
    } else {
        y2n = rep(mean(y1),length(x$alpha))
        axis(side=4, las=3, cex.axis=0.85, col.axis="blue", labels=signif(max(y2,na.rm=T),digits=2), at=mean(y1))
    }
    lines(x$alpha,y2n,type="b",col="blue")
    mtext(text="-log10 pval (model vs null)", side=4, col="blue", line=2.1)
    
    if (is.null(main)) {
        main = "model performance vs alpha values"
    }
    title(main, col.main=col.main, line=line, cex.main=cex.main)
    
    par(mar=par_mar) # restores par settings
}
