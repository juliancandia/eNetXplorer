\name{eNetXplorer-internal}
\title{Internal eNetXplorer functions}
\alias{eNetXplorerBinomial}
\alias{eNetXplorerGaussian}
\alias{eNetXplorerMultinomial}
\alias{plotFeatureCaterpillarMultinom}
\alias{plotFeatureHeatmapMultinom}
\alias{plotMeasuredVsOOBCateg}
\alias{plotMeasuredVsOOBNumer}
\alias{print.eNetXplorer}
\alias{print.summary.eNetXplorer}
\description{Internal eNetXplorer functions}
\usage{
eNetXplorerGaussian(x, y, family, alpha, nlambda, nlambda.ext, seed, scaled, n_fold, 
n_run, n_perm_null, save_lambda_QF_full, QF.FUN, QF_label, QF_gaussian, \dots)

eNetXplorerBinomial(x, y, family, alpha, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, save_lambda_QF_full, QF.FUN, QF_label, binom_method,
binom_pos, fscore_beta, fold_distrib_fail.max, \dots)

eNetXplorerMultinomial(x, y, family, alpha, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, save_lambda_QF_full, QF.FUN, QF_label, multinom_method,
fscore_beta, fold_distrib_fail.max, \dots)

eNetXplorerCox <- function(x, y, family, alpha, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, save_lambda_QF_full, QF.FUN, QF_label, cox_index, \dots)

plotMeasuredVsOOBNumer <- function (x, alpha.index, xlab, ylab, cex.lab, main, 
col.main, cex.main, instance.label, instance.label.cex, instance.label.offset,
instance.label.added.margin, col, \dots)

plotMeasuredVsOOBCateg <- function (x, alpha.index, xlab, ylab, cex.lab, main, 
col.main, cex.main, instance.label, instance.label.cex, instance.label.offset,
box.wex, box.range, box.col,
transparency, jitter, cex.pt, class.color, \dots)

plotFeatureCaterpillarMultinom(x, alpha.index, stat, feature.all, 
feature.pval.thres, feature.set, feature.top.n, signif.code, xlab, ylab, main, 
col.main, cex.main, line, subtitle, col.subtitle, line.subtitle, cex.subtitle, 
cexRow, cex.lab, legend, \dots)

plotFeatureHeatmapMultinom(x, alpha.index, stat, feature.all, feature.pval.thres,
feature.set, feature.top.n, signif.code, xlab, ylab, main, col.main, cex.main, 
line, col, breaks, scale, Rowv, Colv, na.color, cexRow, srtRow, cexCol, srtCol,
margins, key, key.title, dendogram, trace, notecol.freq, notecol.coef, notecex,
subtitle1, col.subtitle1, line.subtitle1, cex.subtitle1, subtitle2, 
col.subtitle2, line.subtitle2, cex.subtitle2, \dots)
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\details{These functions are not part of the package application programming interface (API) and are not intended for use by end users. Functions 
\code{eNetXplorerGaussian}, \code{eNetXplorerBinomial}, \code{eNetXplorerMultinomial} and \code{eNetXplorerCox} are called by function \code{eNetXplorer} to run models of 
type linear, binomial, multinomial and Cox, respectively. Functions \code{plotMeasuredVsOOBNumer} and \code{plotMeasuredVsOOBCateg} are called by function \code{plotMeasuredVsOOB} to generate plots for linear and categorical models, respectively. For multinomial models, functions \code{plotFeatureCaterpillarMultinom} and \code{plotFeatureHeatmapMultinom} are called 
by functions \code{plotFeatureCaterpillar} and \code{plotFeatureHeatmap}, respectively.
}
\keyword{internal}
