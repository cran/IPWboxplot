\name{IPW.ASYM.boxplot}
\alias{IPW.ASYM.boxplot}

\title{
Boxplot adapted to skewness and missing values
}
\description{The function draws a modified boxplot adapted to missing data and skewness. The drop-out probabilities can be given by the practitioner or  fitted through a logistic model using auxiliary covariates.  The plots are adapted to asymmetric distributions by correcting the whiskers through a measure of the skewness.
}
\usage{
IPW.ASYM.boxplot(y,px=NULL,x=NULL,graph=c("IPW","both"),names=c("IPW Asymmetric Boxplot",
                "NAIVE Asymmetric Boxplot"), size.letter=1.2,
                method=c("quartile","octile"), ctea=-4, cteb=3, lim.inf=NULL,lim.sup=NULL,
                main=" ",xlab = " ", ylab ="  ",color="black")
}

\arguments{
  \item{y}{  Numerical vector of length n with possible missing values codified by NA or NAN.}
 \item{px}{Optional. Numerical vector of drop-out probabilities. If not provided a logistic fit is performed using \code{x} as predictive variable. Missing values are not admitted.}

\item{x}{  	Optional. The matrix of fully observed variables used to estimate the missing model with dimension nrows=n and ncol=p. Missing values are not admitted. One of the vectors px or x must be supplied.}
\item{graph}{   Optional. Character string indicating if the plot contains two boxplots ("both") or only the boxplot computed with the inverse probability weighted quantiles ("IPW"). 	The default is "IPW".}
\item{names}{  	Optional. Character string to name the boxplots. The default is "IPW Asymmetric Boxplot", when \code{graph="IPW"} and
                                   c("IPW Asymmetric Boxplot",   "NAIVE Asymmetric Boxplot") when \code{graph="both"}.}
\item{size.letter}{  	Optional. The font size of names. Default value is 1.2}
\item{method}{   		Optional. Character string indicating if the skewness measure is based on the quartiles ("quartile") or the octiles ("octile").
		The default is "quartile". }
	\item{ctea}{ 	Optional. Scaling factors to compute the outlier boundary. The default is -4.}

	\item{cteb}{ Optional. Scaling factors to compute the outlier boundary. The default is 3. When ctea=cteb=0 the IPW boxplot for symmetric data is obtained.}

\item{lim.inf}{  	Optional. The lower limit of the plot if supplied by the user.}
	\item{lim.sup}{ Optional. The upper limit of the plot if supplied by the user.}
	\item{main}{ Optional. Character string to title the plot. By default no main title is given.}
\item{xlab}{Optional. Character string to indicate the label of the horizontal axis.}
	\item{ylab}{ 	Optional. Character string to indicate the label of the vertical axis.}
\item{color}{  Optional. Color for the IPW Boxplot.}

}
\details{

The function draws boxplots designed to adjust both for  skewness and missingness. The drop-out probabilities can be supplied by the user or estimated through a logistic model from given covariates.

The function plots as default a modified boxplot based on the inverse probability weighted (IPW) quantiles adapting for missing observations as in Zhang et al.(2012), but using a correction factor to adjust for skewness. For that purpose, the function incorporates a skewness measure  to compute  the whiskers and the outlier cut--off values in a similar way to that considered in Hubert and Vandervieren (2008).
The argument \code{method} selects quartiles (\code{method="quartile"}) or octiles (\code{method="octile"}) to calculate the skewness measure SKEW, respectively, as
\deqn{ SKEW=\frac{(Q_{0.75}-Q_{0.5})-(Q_{0.5}-Q_{25})}{(Q_{0.75}-Q_{0.25})},}

\deqn{ SKEW=\frac{(Q_{0.875}-Q_{0.5})-(Q_{0.5}-Q_{0.125})}{(Q_{0.875}-Q_{0.125})},}
where \eqn{Q\alpha} denotes the \eqn{\alpha-}quantile.

The whiskers and the outlier cut--off values are computed by means of an exponential model in the fashion of Hubert and Vandervieren (2008) taking into account the interval:

\deqn{(Q_{0.25}-1.5*\exp{(c_i*SKEW)}*IQR,Q_{0.75}+1.5*\exp{(c_s*SKEW)}*IQR),}

where  \eqn{IQR=Q_{0.75}-Q_{0.25}} and    \eqn{c_i}=\code{ctea} and  \eqn{c_s}=\code{cteb} if SKEW is positive, otherwise,  \eqn{c_i}=-\code{cteb} and  \eqn{c_s}=-\code{ctea}.

The default values for \code{ctea} and \code{cteb} are  \eqn{-4} and  \eqn{3}, however, the user may choose other values for these constants.

By specifying \code{graph = "both"}, the function displays two parallel modified boxplots. The boxplot on the left corresponds to the IPW boxplot adapted for missingness and skewness, while that on the right, to its naive counterpart which is simply based on the observations \code{y} at hand without any correction for missingness.


The user can supply a vector of drop-out probabilities \code{px} or a set of covariates \code{x} to estimate the propensity.
When both \code{px} and \code{x} are supplied, the IPW.ASYM.boxplot is executed using \code{px}. When \code{px} is not given, it is estimated assuming a logistic model depending on the covariates \code{x} .
For more details, see Bianco et al. (2018).



}
\value{
The output of the function is a list with components:
  \item{px }{Numerical vector of drop-out probabilities.}
\item{IPW.Quartiles}{Numerical vector of inverse probability weighted quartiles.}
\item{IPW.whisker}{Numerical vector of lower and upper whiskers calculated from IPW quantiles.}
	\item{out.IPW}{Numerical vector of data points detected as atypical by the IPW boxplot adapted to skewness.}
		\item{SKEW.IPW}{Skewness measure based on the IPW quartiles (method="quartile") or IPW octiles (method="octile"). }
\item{NAIVE.Quartiles}{Numerical vector of naive quartiles computed from the subset of non-missing values of \code{y}. Returned only when graph="both".}
	\item{NAIVE.whisker}{Numerical vector of lower and upper whiskers obtained from the naive quantiles. Returned only when graph="both".}
\item{out.NAIVE}{Numerical vector of data points detected as atypical by the naive boxplot. Returned only when graph="both".}
\item{SKEW.NAIVE}{ Skewness measure based on the naive quartiles (method="quartile") or Naive octiles (method="octile"), computed from the subset of non-missing values of \code{y}. Returned only when graph="both".}


}
\references{
Bianco, A. M., Boente, G., and Perez-Gonzalez, A. (2018). A boxplot adapted to missing values: an R function when predictive covariates are available. Submitted.\cr\cr

Hubert, M. and Vandervieren,  E. (2008).  An adjusted boxplot for skewed distributions.  Computational Statistics & Data Analysis, 52, 5186-5201. \cr\cr

Zhang, Z., Chen, Z., Troendle, J. F. and Zhang, J. (2012). Causal inference on quantiles with  an obstetric application. Biometrics, 68, 697-706. \cr\cr

}
\author{
Ana Maria Bianco <abianco@dm.uba.ar>, Graciela Boente <gboente@dm.uba.ar> and Ana Perez-Gonzalez <anapg@uvigo.es>.}
\note{
The missing values of \code{y} must be codified as  NA or NAN.

The numerical vector \code{px} and the matrix of covariates \code{x} must be fully observed.  \code{px} or  \code{x} must be supplied by the user.

The lengths of \code{y}, \code{px}, and \code{nrow(x)} must be equal.



}
\seealso{IPW.quantile,
IPW.Boxplot
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## A real data example


library(mice)
data(boys)
attach(boys)

# The function plots the IPW boxplot adapted to skewness.
# Some statistical summaries computed using the inverse probability weighting approach
# are also returned.
res1=IPW.ASYM.boxplot(hc,x=age,main="IPW boxplot adjusted for skewness of the head circumference")

# We can compare the naive and IPW approaches. We also can consider the skewness measure computed
# using the  quartiles (as default).
res2=IPW.ASYM.boxplot(hc,x=age,method="quartile",graph="both",main=" ")

# The results obtained if the skewness measure is computed with the octiles (method="octile") are:

res3=IPW.ASYM.boxplot(hc,x=age,method="octile",graph="both",main=" ")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


\keyword{ quantile }
\keyword{ boxplot }
\keyword{ missing }
\keyword{ inverse probability weighted}
