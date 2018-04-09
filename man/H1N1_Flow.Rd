\name{H1N1_Flow}
\alias{H1N1}
\docType{data}
\title{Prediction of H1N1 postvaccination response based on longitudinal cell subpopulation frequencies}
\description{
  Data from a cohort of healthy subjects vaccinated against influenza virus H1N1. Cell subpopulation frequencies from deep-phenotyping flow cytometry were determined longitudinally pre- (days -7, 0) and post-vaccination (days 1, 7, 70). The response is the adjusted maximum fold change (adjMFC) of serum titers at day 70 relative to baseline, as defined in Tsang et al (2014) 
}
\usage{data(H1N1_Flow)}

\format{
For each timepoint, a numerical matrix of predictors is provided with subjects as rows and cell subpopulations as columns. The serum titer response is a named numerical vector. A metadata file with cell subpopulation annotations is also provided.}
\references{Tsang, J.S. et al. (2014)
  \emph{Global Analyses of Human Immune Variation Reveal Baseline Predictors of Postvaccination Responses}, \cr
  \emph{Cell 157, 499-513.}\cr
}

\examples{
data(H1N1_Flow)
}
%# run eNetXplorer to assess the ability to predict H1N1 serum 
%# titers at day 70 based on cell subpopulation frequencies 
%at day 7.  
%eNetXplorer(x=predictor_day7, y=response[rownames(predictor_day7)], family="gaussian", n_run=25, n_perm_null=15, seed=123)
%}

\seealso{%\code{eNetXplorer
}

\keyword{datasets}

