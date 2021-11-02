## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("eNetXplorer")

## ---- message=FALSE-----------------------------------------------------------
data_gen_pop_covmat <- function(n_inst, covmat, seed=123) {
    library (expm);
    set.seed(seed)
    data <- matrix(rnorm(n_inst*ncol(covmat)),ncol=ncol(covmat))%*%sqrtm(covmat)
    predictor=data[,-1,drop=F]
    rownames(predictor) = paste0("Inst.",1:n_inst)
    colnames(predictor) = paste0("Feat.",1:(ncol(covmat)-1))
    list(response=data[,1],predictor=predictor)
}

## -----------------------------------------------------------------------------
covmat_gen <- function(n_pred, block_size, r_resp, r_block) {
    covmat = matrix(rep(1.e-3,(n_pred+1)**2),ncol=(n_pred+1))
    for (i_pred in 1:block_size) {
        for (j_pred in (i_pred+1):(block_size+1)) {
            if (i_pred==1) {
                covmat[i_pred,j_pred] = r_resp
            } else {
                covmat[i_pred,j_pred] = r_block
            }
            covmat[j_pred,i_pred] = covmat[i_pred,j_pred]
        }
    }
    for (i_pred in 1:n_pred) {
        covmat[i_pred,i_pred] = 1
    }
    covmat
}

## ---- message=FALSE-----------------------------------------------------------
data = data_gen_pop_covmat(n_inst=50, covmat_gen(n_pred=60, block_size=5, r_resp=0.5, r_block=0.35))

## ---- echo=FALSE, message=FALSE, fig.height = 4.5, fig.width = 6.5, fig.align = "left"----
library(gplots)
library(RColorBrewer)
cor_mat = cor(data$predictor)
cor_max = max(abs(range(cor_mat[-seq(1,ncol(cor_mat)^2,ncol(cor_mat)+1)])))
n_breaks = 10
breaks = seq(-cor_max,cor_max,2*cor_max/(n_breaks-1))
heatmap.2(cor_mat,col=redgreen,breaks=breaks,scale="none",dendrogram="none",Rowv=F,Colv=F,margins=c(5,8),cexRow=0.5,cexCol=0.5,trace="none")

## ---- echo=FALSE, fig.height = 3.6, fig.width = 5.5, fig.align = "left"-------
plot(cor(data$predictor,data$response),xaxt="n",xlab="predictor",ylab="correlation to response")
axis(side=1, at=1:ncol(data$predictor),labels=colnames(data$predictor), las=2, cex.axis=0.5)

## -----------------------------------------------------------------------------
library(eNetXplorer)

## ---- warning=FALSE, eval=FALSE, tidy=TRUE------------------------------------
#  fit_def = eNetXplorer(x=data$predictor,y=data$response,family="gaussian")

## ---- warning=FALSE, eval=FALSE, tidy=TRUE------------------------------------
#  fit = eNetXplorer(x=data$predictor,y=data$response,family="gaussian",alpha=seq(0,1,by=0.1),n_run=1000,n_perm_null=250,seed=123)

## ---- echo=FALSE--------------------------------------------------------------
# OR, for our purposes, we upload the object previously generated:
load("Case1_r3_fit.Robj")

## -----------------------------------------------------------------------------
summary(fit)

## ---- warning=FALSE, fig.height = 4, fig.width = 5.5, fig.align = "left"------
plot(fit, plot.type="summary") 

## ---- warning=FALSE, tidy=TRUE, fig.height = 4.2, fig.width = 5.5, fig.align = "left"----
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type="featureCaterpillar", stat=c("coef"))

## ---- warning=FALSE, tidy=TRUE, fig.height = 4.5, fig.width = 6.5, fig.align = "left"----
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type="featureHeatmap", stat=c("coef"),notecex=1.5) 

## ---- warning=FALSE, tidy=TRUE, fig.height = 4, fig.width = 5.5, fig.align = "left"----
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type="lambdaVsQF") 

## ---- warning=FALSE, tidy=TRUE, fig.height = 4, fig.width = 5.5, fig.align = "left"----
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type="measuredVsOOB") 

## -----------------------------------------------------------------------------
data = data_gen_pop_covmat(n_inst=50, covmat_gen(n_pred=60, block_size=1, r_resp=0.7, r_block=0.35))

## ---- echo=FALSE, fig.height = 4.5, fig.width = 6.5, fig.align = "left"-------
library(gplots)
library(RColorBrewer)
cor_mat = cor(data$predictor)
cor_max = max(abs(range(cor_mat[-seq(1,ncol(cor_mat)^2,ncol(cor_mat)+1)])))
n_breaks = 10
breaks = seq(-cor_max,cor_max,2*cor_max/(n_breaks-1))
heatmap.2(cor_mat,col=redgreen,breaks=breaks,scale="none",dendrogram="none",Rowv=F,Colv=F,margins=c(5,8),cexRow=0.5,cexCol=0.5,trace="none")

## ---- echo=FALSE, fig.height = 3.6, fig.width = 5.5, fig.align ="left"--------
plot(cor(data$predictor,data$response),xaxt="n",xlab="predictor",ylab="correlation to response")
axis(side=1, at=1:ncol(data$predictor),labels=colnames(data$predictor), las=2, cex.axis=0.5)

## ---- warning=FALSE, eval=FALSE, tidy=TRUE------------------------------------
#  fit = eNetXplorer(x=data$predictor,y=data$response,family="gaussian",alpha=seq(0,1,by=0.1),n_run=1000,n_perm_null=250,seed=123)

## ---- echo=FALSE--------------------------------------------------------------
# OR, for our purposes, we upload the object previously generated:
load("Case2_r3_fit.Robj")

## -----------------------------------------------------------------------------
summary(fit)

## ---- warning=FALSE, fig.height = 4, fig.width = 5.5, fig.align ="left"-------
plot(fit, plot.type="summary") 

## ---- warning=FALSE, tidy=TRUE, fig.height = 4.2, fig.width = 5.5, fig.align ="left"----
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type="featureCaterpillar", stat=c("coef"))

## ---- warning=FALSE, tidy=TRUE, fig.height = 4.5, fig.width = 6.5, fig.align ="left"----
plot(fit, alpha.index = which.max(fit$model_QF_est), plot.type="featureHeatmap", stat=c("coef"),notecex=1.5) 

## ---- eval=FALSE--------------------------------------------------------------
#  data(H1N1_Flow)

## ---- eval=FALSE--------------------------------------------------------------
#  data(Leukemia_miR)

## ---- eval=FALSE--------------------------------------------------------------
#  expr_full = Leuk_miR_full$expression_matrix
#  miR_filter = rep(F,nrow(Leuk_miR_full$miR_metadata))
#  miR_filter[apply(expr_full,2,mean)>1.2] = T
#  sample_filter = rep(T,nrow(Leuk_miR_full$sample_metadata))
#  sample_filter[Leuk_miR_full$sample_metadata$sample_class=="Normal"] = F
#  expr_filtered = expr_full[sample_filter,miR_filter]
#  miR_filtered = Leuk_miR_full$miR_metadata[miR_filter,]
#  sample_filtered = Leuk_miR_full$sample_metadata[sample_filter,]

## ---- eval=FALSE--------------------------------------------------------------
#  data(breastCancerSurv)

## ---- echo=FALSE--------------------------------------------------------------
# OR, for our purposes, we upload the object previously generated:
load("eNet_gaussian_covsam_null.Robj")

## ---- warning=FALSE, tidy=TRUE, fig.height = 4, fig.width = 5.5, fig.align = "left"----
plot(eNet, alpha.index = which.max(eNet$model_QF_est), plot.type="measuredVsOOB") 

## ---- message=FALSE-----------------------------------------------------------
data_gen_sampl_covmat <- function(n_inst, covmat) {
    library (expm)
    mat <- matrix(rnorm(n_inst*ncol(covmat)),ncol=ncol(covmat))
    data = mat%*%sqrtm(solve(cov(mat)))%*%sqrtm(covmat)
    predictor=data[,-1,drop=F]
    rownames(predictor) = paste0("Inst.",1:n_inst)
    colnames(predictor) = paste0("Feat.",1:(ncol(covmat)-1))
    list(response=data[,1],predictor=predictor)
}

## ---- message=FALSE-----------------------------------------------------------
n_inst=50
n_fold = 5
foldid = NULL
fold_size = floor(n_inst/n_fold)
for (i_fold in 1:n_fold) {
    foldid = c(foldid,rep(i_fold,fold_size))
}
fold_rest = n_inst%%n_fold
if (fold_rest>0) {
    for (i_fold in 1:fold_rest) {
        foldid = c(foldid,i_fold)
    }
}

## ---- warning=FALSE, message=FALSE--------------------------------------------
n_run = 20
n_rdm = 10
set.seed(123)
cor_data = rep(NA,n_run)
cor_pred_full = rep(NA,n_run)
cor_pred_OOB_mean = rep(NA,n_run)
cor_pred_OOB_sd = rep(NA,n_run)
covmat=matrix(c(1,0,0,1),ncol=2)
for (i_run in 1:n_run) {
    data_gen = data_gen_sampl_covmat(n_inst,covmat)
    data = data.frame(x=data_gen$predictor[,1],y=data_gen$response)
    cor_data[i_run] = cor(data$x,data$y)
    fit = lm("y~x",data)
    y_pred = predict(fit,data)
    cor_pred_full[i_run] = cor(data$y,y_pred)
    # we generate random folds
    cor_pred_OOB = rep(NA,n_rdm)
    for (i_rdm in 1:n_rdm) {
        foldid_rdm = sample(foldid)
        y_pred_rdm = rep(NA,n_inst)
        for (i_fold in 1:n_fold) {
            IB = foldid_rdm!=i_fold
            fit = lm("y~x",data[IB,])
            y_pred_rdm[!IB] = predict(fit,data[!IB,])
        }
        cor_pred_OOB[i_rdm] = cor(data$y,y_pred_rdm)
    }
    cor_pred_OOB_mean[i_run] = mean(cor_pred_OOB)
    cor_pred_OOB_sd[i_run] = sd(cor_pred_OOB)
}

## ---- warning=FALSE, tidy=TRUE, fig.height = 4, fig.width = 5.5, fig.align = "left"----
color = c("black","blue","red")
x = 1:n_run
range_y = range(cor_data,cor_pred_full,cor_pred_OOB_mean+cor_pred_OOB_sd,
                cor_pred_OOB_mean-cor_pred_OOB_sd,na.rm=T)
plot(range(x),range_y,type="n",xlab="run",ylab="Pearson's correlation",cex.lab=1.)
lines(x,cor_data,type="b",pch=16,col=color[1],lty=3)
lines(x,cor_pred_full,type="b",pch=16,col=color[2],lty=3)
lines(x,cor_pred_OOB_mean,type="b",pch=16,col=color[3],lty=3)
arrows(x,cor_pred_OOB_mean-cor_pred_OOB_sd,x,cor_pred_OOB_mean+cor_pred_OOB_sd,col=color[3],length=0.025,angle=90,code=3)
legend_txt = c("cor(y,x)","cor(y,y_pred)","cor(y,y_pred_OOB)")
legend(11,-0.075,legend_txt,pch=16,lwd=1,col=color,lty=3,cex=0.7,pt.cex=0.7) 

## ---- eval=FALSE--------------------------------------------------------------
#  library(furrr)
#  data(QuickStartEx)
#  alpha_values = seq(0,1,by=0.2)
#  future::plan(multiprocess) # Set up parallel processing
#  eNet = alpha_values %>% furrr::future_map(~ eNetXplorer(
#    x=QuickStartEx$predictor, y=QuickStartEx$response,family="gaussian",
#    n_run=20,n_perm_null=10,save_obj=T,dest_obj=paste0("eNet_a",.x,".Robj"),alpha=.x))
#  mergeObj(paste0("eNet_a",alpha_values,".Robj"))

