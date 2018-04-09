#eNetXplorer <- function(x, ...) UseMethod("eNetXplorer")

#eNetXplorer.default <- function(x, y, family, alpha=seq(0,1,by=0.2), nlambda=100, nlambda.ext=NULL,
eNetXplorer <- function(x, y, family, alpha=seq(0,1,by=0.2), nlambda=100, nlambda.ext=NULL,
seed=NULL, scaled=T, n_fold=5, n_run=100, n_perm_null=25, QF.FUN=NULL, QF_label=NULL, cor_method="pearson", fold_distrib_fail.max=100, ...)
{
    if (family=="gaussian") {
        res <- eNetXplorerGaussian(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext,seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,n_perm_null=n_perm_null, QF.FUN=QF.FUN, QF_label=QF_label, cor_method=cor_method, ...)
    } else if (family=="binomial") {
        res <- eNetXplorerBinomial(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext, seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, QF.FUN=QF.FUN, QF_label=QF_label, fold_distrib_fail.max=fold_distrib_fail.max, ...)
    } else if (family=="multinomial") {
        res <- eNetXplorerMultinomial(x=x, y=y, family=family, alpha=alpha, nlambda=nlambda, nlambda.ext=nlambda.ext, seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, QF.FUN=QF.FUN, QF_label=QF_label, fold_distrib_fail.max=fold_distrib_fail.max, ...)
    } else {
        stop("Family type not supported\n")
    }
    res$call <- match.call()
    class(res) <- "eNetXplorer"
    res
}
