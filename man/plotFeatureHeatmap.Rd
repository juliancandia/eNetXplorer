\name{plotFeatureHeatmap}
\alias{plotFeatureHeatmap}
\title{generates heatmap plot of feature statistics}
\description{This function generates a heatmap plot of feature statistics across \code{alpha}, which includes statistical significance annotations for individual features. By default, features are selected (and ordered top-down) by statistical significance based on a given value of \code{alpha}; options are provided to customize feature selection and display.}
\usage{
plotFeatureHeatmap(x, alpha.index=NULL, stat=c("freq","coef"), feature.all, feature.pval.thres, 
feature.set, feature.top.n, signif.code, xlab=NULL, ylab=NULL, main=NULL, 
col.main="black", cex.main=0.95, line=1, col=NULL, breaks=NULL, scale="none", 
Rowv=F, Colv=F, na.color=NULL, cexRow=NULL, srtRow=0, cexCol=0.75, srtCol=45, 
margins=c(5, 5), key=T, key.title=NA, dendogram="none", trace="none", 
notecol.freq="black", notecol.coef="white", notecex=1, subtitle1=NULL, 
col.subtitle1="black", line.subtitle1=-1, cex.subtitle1=0.65, subtitle2=NULL, 
col.subtitle2="darkgray", line.subtitle2=-2, cex.subtitle2=0.55, \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{alpha.index}{Integer indices to select \code{alpha} values. Default is \code{1:length(alpha)}}
  \item{stat}{
  Feature statistic: \code{"freq"} for nonzero frequency, \code{"coef"} for mean nonzero coefficient.}
  \item{feature.all}{(Feature selection argument 1) Logical to show all features. Default is FALSE.}
  \item{feature.pval.thres}{(Feature selection argument 2) P-value threshold to select features. Default is NULL.}
  \item{feature.set}{(Feature selection argument 3) Character vector of feature names to display. Default is NULL.}
  \item{feature.top.n}{(Feature selection argument 4) Number of top features (ordered by p-value) to display. Default is 25.}
  \item{signif.code}{Logical to display statistical significance annotations. Default is TRUE.}
  \item{xlab}{Label for x axis.}
  \item{ylab}{Label for y axis.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{line}{Title location.}
  \item{col}{Heatmap color vector. Length must be one less than number of breaks.}
  \item{breaks}{Color breaks vector. Default number of breaks is 10.}
  \item{scale}{Logical to scale the data for heatmap in either the row or column direction. Default is \code{"none"}.}
  \item{Rowv}{Logical to reorder rows by hierarchical clustering. Default is FALSE.}
  \item{Colv}{Logical to reorder columns by hierarchical clustering. Default is FALSE.}
  \item{na.color}{Color to use for missing values.}
  \item{cexRow}{Size of row labels.}
  \item{srtRow}{Angle of row labels, in degrees from horizontal.}
  \item{cexCol}{Size of column labels.}
  \item{srtCol}{Angle of column labels, in degrees from horizontal.}
  \item{margins}{Numeric vector of length 2 containing the margins for column and row names, respectively.}
  \item{key}{Logical to display key. Default is TRUE.}
  \item{key.title}{Main title of the color key.}
  \item{dendogram}{To draw dendograms. Default is \code{"none"}.}
  \item{trace}{To display trace lines. Default is \code{"none"}.}
  \item{notecol.freq}{Color of statistical significance annotations for feature frequency heatmaps.}
  \item{notecol.coef}{Color of statistical significance annotations for feature coefficient heatmaps.}
  \item{notecex}{Size of significance annotations.}
  \item{subtitle1}{Custom subtitle 1.}
  \item{col.subtitle1}{Color of subtitle 1.}
  \item{line.subtitle1}{Position of subtitle 1.}
  \item{cex.subtitle1}{Size of subtitle 1.}
  \item{subtitle2}{Custom subtitle 2.}
  \item{col.subtitle2}{Color of subtitle 2.}
  \item{line.subtitle2}{Position of subtitle 2.}
  \item{cex.subtitle2}{Size of subtitle 2.}
\item{\dots}{Additional plotting parameters.}
}
\details{Feature selection criteria are hierarchical based on arguments 1 through 4 (see argument description above), with argument 1 at the top of the hierarchy. E.g. if \code{feature.all} is explicitly set to TRUE, it will take precedence over any other feature selection argument. By default, the top 25 features are displayed, ordered top-down by significance based on the given value of \code{alpha}.}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{eNetXplorer}}, \code{\link{plot}}}
\examples{
\dontrun{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",n_run=20,n_perm_null=10,seed=111)
plot(x=fit,plot.type="featureHeatmap",stat="freq")
plotFeatureHeatmap(x=fit,alpha.index=3,stat="freq",main="custom title")
}
}
