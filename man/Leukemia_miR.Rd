\name{Leukemia_miR}
\alias{Leukemia_miR}
\docType{data}
\title{MicroRNA expression of acute leukemia phenotypes}
\description{
Data of human microRNA (miR) expression from acute myeloid (AML) and acute lymphoblastic (ALL) leukemia cell lines and primary (patient) samples. Acute lymphoblastic leukemia samples are further classified by B-cell (B-ALL) and T-cell (T-ALL) subphenotypes.
}
\usage{data(Leukemia_miR)}

\format{
A numerical matrix of predictors is provided with samples as rows and miRs as columns. Two categorical response vectors are provided for binomial (AML, ALL) and multinomial (AML, B-ALL, T-ALL) classification.
}

\details{
An initial set of 847 human miRs was measured in multiple acute leukemia cell lines and primary samples. 
Exclusion criteria based on miR expression with mean intensity <1 or >14 (after RMA normalization and log2-transformation) 
led to a final set of 370 miRs. For more details, see Candia et al (2015).
}

\references{Candia J et al (2015)
  \emph{Uncovering low-dimensional, miR-based signatures of acute myeloid and lymphoblastic leukemias with a machine-learning-driven network approach},
  Converg Sci Phys Oncol 1: 025002.\cr
}

\keyword{datasets}

