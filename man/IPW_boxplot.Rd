\name{IPW.boxplot}
\alias{IPW.boxplot}

\title{

Boxplot adapted to missing values

}
\description{
The function draws a modified boxplot adapted to missing values. The drop-out probabilities can be given by the practitioner or  fitted through a logistic model using auxiliary covariates.
The function returns the usual boxplot of the available data as well as a modified plot which takes into account the missing data model and weights the observations using the estimated/given propensity.
}
\usage{
IPW.boxplot(y,px=NULL,x=NULL,graph=c("IPW","both"),
names=c("IPW  Boxplot", "NAIVE  Boxplot"), size.letter=1.2,
lim.inf=NULL,lim.sup=NULL,main=" ",xlab = " ", ylab =" ",color="black")
}
\arguments{
  \item{y}{  Numerical vector of length n with possible missing values codified by NA or NAN.}
 \item{px}{Optional. Numerical vector of drop-out probabilities. If not provided a logistic fit is performed using \code{x} as predictive variable. Missing values are not admitted.}

\item{x}{  	Optional. The matrix of fully observed variables used to estimate the missing model with dimension nrows=n and ncol=p. Missing values are not admitted. One of the vectors px or x must be supplied.}

\item{graph}{   Optional. Character string indicating if the plot contains two boxplots ("both") or only the boxplot computed with the inverse probability weighted quantiles("IPW"). 	The default is "IPW".}
\item{names}{  	Optional. Character string to name the boxplots. The default is "IPW  Boxplot", when \code{graph="IPW"} and
                                   c("IPW  Boxplot",   "NAIVE  Boxplot") when \code{graph="both"}.}
\item{size.letter}{  	Optional. The font size of names. Default value is 1.2}
\item{lim.inf}{  	Optional. The lower limit of the plot if supplied by the user.}
	\item{lim.sup}{ Optional. The upper limit of the plot if supplied by the user.}
	\item{main}{ Optional. Character string to title the plot. By default no main title is given.}
  \item{xlab}{Optional. Character string to indicate the label of the horizontal axis.}
	\item{ylab}{ 	Optional. Character string to indicate the label of the vertical axis.}
\item{color}{  Optional. Color for the IPW Boxplot.}

}
\details{
The function draws boxplots designed to adjust for missing values. The propensity can be supplied by the user or estimated through a logistic model from given covariates.

The function plots as default a modified boxplot based on the inverse probability weighted (IPW) quantiles adapting for missing observations as in Zhang et al.(2012).

By specifying \code{graph = "both"}, the function displays two parallel boxplots. The boxplot on the left corresponds to the IPW boxplot adapted for missingness, while on the right, the naive boxplot, i.e., the usual boxplot simply computed with the observations \code{y} at hand, is displayed.


The user can supply a vector of probabilities \code{px} or a set of covariates \code{x} to estimate it.
When both \code{px} and \code{x} are supplied, the IPW.boxplot is executed using \code{px}. When \code{px} is not supplied, it is estimated assuming a logistic model depending on the covariates \code{x} .
For more details, see Bianco et al. (2018).
}
\value{
The output of the function is a list with components:
  \item{px }{Numerical vector of probabilities.}
\item{IPW.Quartiles}{Numerical vector of inverse probability weighted quartiles.}
\item{IPW.whisker}{Numerical vector of lower and upper whiskers calculated from IPW quartiles.}
	\item{out.IPW }{Numerical vector of data points detected as atypical by the IPW boxplot.}
	\item{ NAIVE.Quartiles }{Numerical vector of naive quartiles computed from the subset of non-missing values of \code{y}. Returned only when graph="both".}
\item{NAIVE.whisker }{Numerical vector of lower and upper whiskers obtained from the naive quantiles. Returned only when graph="both".}
\item{out.NAIVE }{Numerical vector of data points detected as atypical by the naive boxplot. Returned only when graph="both".}
}

\references{
Bianco, A. M., Boente, G., and Perez-Gonzalez, A. (2018). A boxplot adapted to missing values: an R function when predictive covariates are available. Submitted.\cr\cr

Zhang, Z., Chen, Z., Troendle, J. F. and Zhang, J. (2012). Causal inference on quantiles with an obstetric application. Biometrics, 68, 697-706. \cr\cr

}
\author{
Ana Maria Bianco <abianco@dm.uba.ar>, Graciela Boente <gboente@dm.uba.ar> and Ana Perez-Gonzalez <anapg@uvigo.es>.}
\note{
The missing values of \code{y} must be codified as  NA or NAN.

The numerical vector \code{px} and the matrix of covariates \code{x} must be fully observed.  \code{px} or  \code{x} must be supplied by the user.

The lengths of \code{y}, \code{px}, and \code{nrow(x)} must be equal.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{IPW.quantile, IPW.ASYM.Boxplot
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


## A real data example

library(mice)
data(boys)
attach(boys)

res1=IPW.boxplot(tv,x=age,main="IPW boxplot of the testicular volume")


# We  can compare the naive and IPW boxplots
res2=IPW.boxplot(tv,x=age,graph="both",main=" ")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ quantile }
\keyword{ boxplot }
\keyword{ missing }
\keyword{ inverse probability weighted}
