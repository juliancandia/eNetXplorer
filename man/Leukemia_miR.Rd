\name{Leukemia_miR}
\alias{Leukemia}
\docType{data}
\title{Prediction of acute leukemia phenotypes based on microRNA expression}
\description{
   Data of human microRNA expression from acute myeloid (AML) and acute lymphoblastic (ALL) leukemia cell lines and primary (patient) samples. Acute lymphoblastic leukemia samples are further classified by B-cell and T-cell subphenotypes.
}
\usage{data(Leukemia_miR)}

\format{
A numerical matrix of predictors is provided with samples as rows and microRNAs as columns. Two categorical response vectors are provided for binomial (AML, ALL) and multinomial (AML, B-ALL, T-ALL) classification.
}
\references{Candia,J. et al. (2015)
  \emph{Uncovering low-dimensional, miR-based signatures of acute myeloid and lymphoblastic leukemias with a machine-learning-driven network approach}, \cr
  \emph{Converg. Sci. Phys. Oncol. 1, 025002.}\cr
}

\examples{
data(Leukemia_miR)
}
%# Binomial eNetXplorer models to predict acute myeloid (AML) 
%vs acute lymphoblastic (ALL) leukemias: 
%eNetXplorer(x=predictor,y=response_binomial,family="binomial",
%n_run=25, n_perm_null=15, seed=123)

%# Multinomial eNetXplorer models to predict acute myeloid (AML),
%acute B-cell lymphoblastic (B-ALL) and acute T-cell 
%lymphoblastic (T-ALL) leukemias:  
%eNetXplorer(x=predictor,y=response_multinomial,
%family="multinomial", n_run=25, n_perm_null=20, seed=111)

%# Binomial eNetXplorer models to predict B-ALL vs T-ALL: 
%eNetXplorer(x=predictor[response_multinomial!="AML",],
%y=response_multinomial[response_multinomial!="AML"],
%family="binomial", n_run=25, n_perm_null=15, seed=123)
%}

\seealso{
%\code{eNetXplorer
}

\keyword{datasets}

