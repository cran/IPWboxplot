\name{IPW.quantile}
\alias{IPW.quantile}

\title{
Computes the IPW quantiles
}
\description{
The function calculates the inverse probability weighted quantiles of a numeric vector.
}
\usage{
IPW.quantile(y, px=NULL,x=NULL,probs = seq(0, 1, 0.25))
}

\arguments{
    \item{y}{  Numerical vector of length n with possible missing values codified by NA or NAN.}
 \item{px}{Optional. Numerical vector of drop-out probabilities. If not provided a logistic fit is performed using \code{x} as predictive variable. Missing values are not admitted.}

\item{x}{  	Optional. The matrix of fully observed variables used to estimate the missing model with dimension nrows=n and ncol=p. Missing values are not admitted. One of the vectors px or x must be supplied.}

\item{probs}{	Required. Numeric vector of probabilities with values in (0,1).}
}
\details{
The function computes  inverse probability weighted (IPW) quantiles of a numeric vector \code{y} adapting for missing observations as in Zhang et al.(2012).

The user can supply a vector of drop-out probabilities \code{px} or a set of covariates \code{x} to estimate the propensity.
When both \code{px} and \code{x} are supplied, the IPW.quantile is executed using \code{px}. When \code{px} is not supplied, the happenstance probabilities are estimated assuming a logistic model depending on the covariates \code{x}.
For more details, see Bianco et al. (2018).

We adapted the function \code{weighted.fractile} from the \pkg{isotone} package to  missing values in variable \code{y}. See  \pkg{isotone} for  more details.
}
\value{
The output of the function is a list with components:\describe{
  \item{ipw.quantile }{ Numerical vector of length \code{length(probs)} containing the estimated quantiles.}
  \item{px }{Numerical vector of drop-out probabilities.}
}
}
\references{Bianco, A. M., Boente, G. and Perez-Gonzalez, A. (2018). A boxplot adapted to missing values: an R function when predictive covariates are available. Submitted.\cr\cr

Zhang, Z., Chen, Z., Troendle, J. F. and Zhang, J. (2012). Causal inference on quantiles with an obstetric application. Biometrics, 68, 697-706. \cr\cr
}
\author{
Ana Maria Bianco <abianco@dm.uba.ar>, Graciela Boente <gboente@dm.uba.ar> and Ana Perez-Gonzalez <anapg@uvigo.es>.
}
\note{
The missing values of \code{y} must be codified as  NA or NAN.

The numerical vector \code{px} and the matrix of covariates \code{x} must be fully observed.  \code{px} or  \code{x} must be supplied by the user.

The lengths of \code{y}, \code{px}, and \code{nrow(x)} must be equal.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{



## A real data example
library(mice)
data(boys)
attach(boys)
# As an illustration, we consider variable testicular volume, tv.
# To compute the inverse probability weighted (IPW) quartiles
# the covariate age is considered as covariate with predictive capability
# to estimate the vector of drop-out probabilities.

res=IPW.quantile(tv,x=age,probs=c(0.25,0.5,0.75))
res$IPW.quantile

# Compute the inverse probability weighted (IPW) quantiles
# corresponding to the fractiles  0.3, 0.8 and 0.9
# using the covariate age  to estimate the propensity.


res1=IPW.quantile(tv,x=age,probs=c(0.3,0.8,0.9))
res1$IPW.quantile

}
\keyword{ quantile }
\keyword{ missing }
\keyword{ inverse probability weighted}
