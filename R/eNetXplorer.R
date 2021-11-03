eNetXplorer <- function(x, y, family=c("gaussian","binomial","multinomial","cox"),
    alpha=seq(0,1,by=0.2), nlambda=100, nlambda.ext=NULL,
    seed=NULL, scaled=TRUE, n_fold=5, n_run=100, n_perm_null=25,
    save_obj=FALSE, dest_dir=getwd(), dest_dir_create=TRUE, dest_dir_create_recur=FALSE, dest_obj="eNet.Robj",
    save_lambda_QF_full=FALSE,
    QF.FUN=NULL, QF_label=NULL, QF_gaussian=c("cor.pearson","cor.spearman","cor.kendall","mse"),
    binom_method=c("accuracy","precision","recall","Fscore","specificity","auc"),
    multinom_method=c("avg accuracy", "avg precision","avg recall","avg Fscore"),
    binom_pos=NULL, fscore_beta=NULL, fold_distrib_fail.max=100, cox_index=c("concordance","D-index"),
    logrank=FALSE, survAUC=FALSE, survAUC_time=NULL,
    ...)
{
    if (save_obj) {
        if ((!dir.exists(dest_dir))&&(dest_dir_create)) {
            dir.create(dest_dir, showWarnings=FALSE, recursive=dest_dir_create_recur)
        }
    }
    family = match.arg(family)
    if (family=="gaussian") {
        eNet <- eNetXplorerGaussian(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext,seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, save_lambda_QF_full=save_lambda_QF_full, QF.FUN=QF.FUN,
        QF_label=QF_label, QF_gaussian=match.arg(QF_gaussian), ...)
    } else if (family=="binomial") {
        eNet <- eNetXplorerBinomial(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext, seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, save_lambda_QF_full=save_lambda_QF_full, QF.FUN=QF.FUN, QF_label=QF_label,
        binom_method=match.arg(binom_method), binom_pos=binom_pos, fscore_beta=fscore_beta,
        fold_distrib_fail.max=fold_distrib_fail.max, ...)
    } else if (family=="multinomial") {
        eNet <- eNetXplorerMultinomial(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext, seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, save_lambda_QF_full=save_lambda_QF_full, QF.FUN=QF.FUN, QF_label=QF_label,
        multinom_method=match.arg(multinom_method), fscore_beta=fscore_beta,
        fold_distrib_fail.max=fold_distrib_fail.max, ...)
    } else if (family=="cox") {
        eNet <- eNetXplorerCox(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext,seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, save_lambda_QF_full=save_lambda_QF_full, QF.FUN=QF.FUN, QF_label=QF_label, cox_index=match.arg(cox_index), logrank=logrank, survAUC=survAUC, survAUC_time=survAUC_time, ...)
    }
    eNet$call <- match.call()
    class(eNet) <- "eNetXplorer"
    if (save_obj) {
        save(eNet,file=file.path(dest_dir,dest_obj))
    }
    eNet
}
