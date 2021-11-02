plotKaplanMeier <- function (
x, alpha.index=NULL, xlab="Time", ylab="Probability of Survival", cex.lab=1,
main=NULL, col.main="black", cex.main=0.95, conf.int=TRUE,
breaks_ptiles=NULL, risk.col=NULL, legend=TRUE, legend.cex=0.75,
...)
{
    if (is.null(alpha.index)) {
        alpha.index = 1:length(x$alpha)
    }
    y = x$response
    if (is.null(breaks_ptiles)) {
        breaks_ptiles = 0.5
    } else {
        breaks_ptiles = breaks_ptiles[order(breaks_ptiles)] # to enforce increasing order
    }
    n_breaks = length(breaks_ptiles)
    if (is.null(risk.col)) {
        if (n_breaks==1) {
            risk.col = c("blue","red")
        } else if (n_breaks==2) {
            risk.col = c("blue","orange","red")
        } else {
            risk.col = 1:(n_breaks+1)
        }
    }
    risk_label = NULL
    breaks_ranges = signif(c(0,breaks_ptiles*100,100),digits=2)
    for (i_breaks in 1:(n_breaks+1)) {
        risk_label = c(risk_label,paste0("Risk ",breaks_ranges[i_breaks],"-",breaks_ranges[i_breaks+1],"%"))
    }
    for (i_alpha in alpha.index) {
        pred_mean = x$predicted_values[[i_alpha]][,1]
        breaks = quantile(pred_mean,probs=breaks_ptiles)
        score = rep(1,length(x$instance))
        for (i_breaks in 1:n_breaks) {
            score[pred_mean>breaks[i_breaks]]=i_breaks+1
        }
        KM_data = data.frame(y,as.factor(score))
        colnames(KM_data)[3] = "x"
        survfit_obj <- survfit(Surv(time,status)~x, data=KM_data)
        if (conf.int) {
            lty = c(1,3,3)
            lwd = c(2,1,1)
        } else {
            lty = 1
            lwd = 2
        }
        plot(survfit_obj,col=risk.col,xlab=xlab,ylab=ylab,conf.int=conf.int,lty=rep(lty,n_breaks+1),
        lwd=rep(lwd,n_breaks+1))
        if (is.null(main)) {
            if ((x$logrank)&(isTRUE(all.equal(breaks_ptiles,0.5)))) {
                main.title = paste0("alpha=",x$alpha[i_alpha]," ; logrank p-val=",signif(x$logrank_pval[i_alpha],digits=2)) # we print the cross-validated logrank p-value
            } else {
                main.title = paste0("alpha=",x$alpha[i_alpha])
            }
        } else {
            main.title = main
        }
        title(main=main.title,cex.main=cex.main,col.main=col.main)
        if (legend) {
            legend("topright",legend=risk_label,col=risk.col,lty=1,lwd=2,cex=legend.cex)
        }
    }
}
