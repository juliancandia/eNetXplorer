\name{plot}
\alias{plot.eNetXplorer}
\alias{plot}
\title{generates plots from eNetXplorer object}
\description{
This function is a wrapper for a variety of plots, namely:

\code{summary}:
model performance across \code{alpha} (to assess the relative performance among 
different member models in the elastic net family, as well as in relation to permutation null models);

\code{lambdaVsQF}:
given \code{alpha}, quality function across \code{lambda} (to examine the selection of the optimal penalty parameter);

\code{measuredVsOOB}:
(for gaussian and categorical models) given \code{alpha}, response vs out-of-bag predictions across instances (to assess 
individual instances, examine outliers, etc);

\code{contingency}:
(for categorical models) given \code{alpha}, response vs out-of-bag predictions across classes;

\code{featureCaterpillar}: 
given \code{alpha}, caterpillar plot of feature statistics compared to permutation 
null models (with statistical significance annotations for individual features);

\code{featureHeatmap}:
heatmap of feature statistics across \code{alpha} (including statistical significance 
annotations for individual features);

\code{KaplanMeier}: 
(for Cox regression models) given \code{alpha}, Kaplan-Meier plot of survival probability as a function of time (where the cohort is partitioned in two or more groups based on predicted risk); and

\code{survROC}: 
(for Cox regression models) given \code{alpha}, time-dependent ROC plot(s) based on predicted risk at the specified timepoints of interest.
}
\usage{
\method{plot}{eNetXplorer}(x, plot.type=c("summary","lambdaVsQF","measuredVsOOB","contingency",
"featureCaterpillar","featureHeatmap","KaplanMeier","survROC"), alpha.index, 
stat=c("freq","coef"), \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{plot.type}{Type of plot to be produced. Available plots 
  are \code{"summary"}, \code{"lambdaVsQF"}, \code{"measuredVsOOB"} (gaussian and categorical models only), \code{"contingency"} (categorical models only), 
    \code{"featureCaterpillar"}, \code{"featureHeatmap"},
    \code{"KaplanMeier"} (Cox models only) and \code{"survROC"} (Cox models only).}
  \item{alpha.index}{Integer indices to select \code{alpha} values. Default is \code{1:length(alpha)}}
  \item{stat}{
  Feature statistic: \code{"freq"} for nonzero frequency, \code{"coef"} for mean nonzero coefficient. 
  Used for plot types \code{"featureHeatmap"} and \code{"featureCaterpillar"}, ignored otherwise.}
  \item{\dots}{Additional plotting parameters.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{eNetXplorer}}, \code{\link{plotSummary}}, \code{\link{plotLambdaVsQF}}, \code{\link{plotMeasuredVsOOB}}, \code{\link{plotContingency}},

\code{\link{plotFeatureCaterpillar}}, \code{\link{plotFeatureHeatmap}}, \code{\link{plotKaplanMeier}}, \code{\link{plotSurvROC}}
}
\examples{
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor, y=QuickStartEx$response,
family="gaussian", n_run=20, n_perm_null=10, seed=111)
plot(x=fit,plot.type="summary")
plot(x=fit,plot.type="lambdaVsQF",alpha.index=2)
plot(x=fit,plot.type="measuredVsOOB",alpha.index=c(1,3,5))
plot(x=fit,plot.type="featureCaterpillar",stat="coef")
plot(x=fit,plot.type="featureHeatmap",stat="freq")
}
}

