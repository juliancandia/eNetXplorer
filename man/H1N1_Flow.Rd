\name{H1N1_Flow}
\alias{H1N1_Flow}
\docType{data}
\title{Longitudinal cell population frequencies and titer response upon H1N1 vaccination}
\description{
Data from a cohort of healthy subjects vaccinated against influenza virus H1N1. Cell population frequencies from deep-phenotyping flow cytometry were determined longitudinally pre- (days -7, 0) and post-vaccination (days 1, 7, 70). The response is the adjusted maximum fold change (adjMFC) of serum titers at day 70 relative to baseline, as defined in Tsang et al (2014). 
}
\usage{data(H1N1_Flow)}

\format{
For each timepoint (days -7, 0, 1, 7, 70), a numerical matrix of predictors is provided with subjects as rows and cell populations as columns. Two versions of the serum titer response are given: \code{response_numer} 
as a numerical vector and \code{response_class} as a categorical vector 
discretized into low (0), 
intermediate (1) and high (2) response classes. A metadata file with cell population annotations is also provided.}

\details{
Cell populations were manually gated and expressed as percent of parent. 
Samples and cell populations were filtered independently for each timepoint. 
Samples filter: excluded if median of viable cells fraction across all 5 tubes was <0.7. 
Cell population filter: excluded if >80\% of samples had <20 cells. 
Data adjustment: data were log10-transformed and pooled across all timepoints, then adjusted for age, gender and ethnicity effects. 
For more details, see Tsang et al (2014).
}

\references{Tsang JS et al (2014)
  \emph{Global Analyses of Human Immune Variation Reveal Baseline Predictors of Postvaccination Responses},
  Cell 157: 499-513.\cr
}

\keyword{datasets}

