plotSurvROC <- function (
x, alpha.index=NULL, survAUC_time,
xlab="False positive rate (1 - Specificity)", ylab="True positive rate (Sensitivity)", cex.lab=1,
main=NULL, col.main="black", cex.main=0.95, status0="censored", status1="events", ...)
{
    if (!"package:survival"%in%search()) {
        attachNamespace("survival") # to solve some issues with timeROC
    }
    
    if (is.null(alpha.index)) {
        alpha.index = 1:length(x$alpha)
    }
    if (is.null(survAUC_time)) {
        stop("Error: survAUC_time must be provided")
    }
    y = x$response
    for (i_alpha in alpha.index) {
        pred_mean = x$predicted_values[[i_alpha]][,1]
        for (pred_time in survAUC_time) {
            survROC = timeROC(T=y[,"time"],delta=y[,"status"],cause=1,marker=pred_mean,times=pred_time,...)
            plot(survROC$FP, survROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),
            xlab=xlab, ylab=ylab)
            abline(0,1,lty=3)
            if (is.null(main)) {
                n_status0 = sum((y[,"status"]==0)&(y[,"time"]<=pred_time))
                n_status1 = sum((y[,"status"]==1)&(y[,"time"]<=pred_time))
                n_surv = sum(y[,"time"]>pred_time)
                n_total = n_status0+n_status1+n_surv
                perc_status0 = round(100*n_status0/n_total)
                perc_status1 = round(100*n_status1/n_total)
                perc_surv = round(100*n_surv/n_total)
                main_top = paste0("alpha=",x$alpha[i_alpha]," ; timepoint=",pred_time)
                main_center = paste0(status0,": n=",n_status0," (",perc_status0,"%) ; ",status1,": n=",n_status1," (",perc_status1,"%) ; survivors: n=",n_surv," (",perc_surv,"%)")
                if ((x$survAUC)&(identical(x$survAUC_time,survAUC_time))) {
                    AUC_025 = round(x$AUC_perc025[which(survAUC_time==pred_time),i_alpha],2)
                    AUC_500 = round(x$AUC_perc500[which(survAUC_time==pred_time),i_alpha],2)
                    AUC_975 = round(x$AUC_perc975[which(survAUC_time==pred_time),i_alpha],2)
                    main_bottom = paste0("AUC=",AUC_500," (95% CI: ",AUC_025,"-",AUC_975,")")
                } else {
                    main_bottom = paste0("AUC=",round(survROC$AUC,2))
                }
                main.title = paste0(main_top,"\n",main_center,"\n",main_bottom)
            } else {
                main.title = main
            }
            title(main=main.title,cex.main=cex.main,col.main=col.main)
        }
    }
}
