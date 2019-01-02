library(isotone)

########################################################################
# FUNCTION USED TO DRAW THE BOXPLOTS ADAPTED TO MISSING VALUES         #
########################################################################
IPW.quantile=function(y,px=NULL,x=NULL,  probs = seq(0, 1, 0.25))
{
  ############################################################################################################################################
  # Arguments
  #
  #   y		Required. Numerical vector of values with possible missing values codified NA or NAN with length n.
  #   px	Optional. Numerical vector of probabilities. If not provided a logistic fit is performed using x
  #   x 	Optional. The matrix of fully observed variables used to estimate the missing model with dimension nrows=n and ncol=dimension.
  #   probs	Required. Numeric vector of probabilities with values in [0,1]
  ############################################################################################################################################

  ########################################################################
  # Value
  #  ipw.quantile	A vector of length length(probs) is returned
  #  px 		Numerical vector of probabilities.
  ########################################################################

  #########################################################################
  #---Preliminary checks---
  #########################################################################

  dimension=NCOL(x)

  if (is.null(x)=="TRUE" & is.null(px)=="TRUE")
    stop("ERROR: It is neccesary to supply the vector of dropouts probabilities for each observation or a covariate to estimate it")

  if (is.null(px)=="FALSE")
  {
    if (min(px)<=0)
      stop("ERROR: px take positive values")
    if (NROW(y)!=NROW(px))
      stop("ERROR: 'y' and 'px' have different lengths")
    if (sum(is.na(px))+sum(is.nan(px))>0)
      stop(" ERROR: px has missing values")
  }


  if (is.null(px)=="TRUE")
  {
    if (NROW(y)!=NROW(x))
      stop("ERROR: 'y' and 'x' have different lengths")
    if (sum(is.na(x))+sum(is.nan(x))>0)
      stop(" ERROR: The  covariates matrix has missing observations")
  }

  if (dimension==1){x=as.vector(x)}

  if (sum(is.na(y))+sum(is.nan(y))==length(y))
    stop(" ERROR: All values are missing")

  nsamp=length(y)

  delta=rep(1,nsamp)
  for (i in 1: nsamp){
    if (is.na(y[i])=="TRUE"|is.nan(y[i])=="TRUE")
    {
      delta[i]<-0
      y[i]=NA
    }
  }


  if(is.null(px)=="FALSE"){METHOD="The dropout probability is given"}

  if(is.null(px)){
    ###############################################################
    #  Estimation of the dropout probability for each observation #
    ###############################################################
    a=glm(delta~x,family="binomial")
    px=a$fitted.values
    METHOD="LOGISTIC"
  }

  ###############################################################
  #  Estimation of the IPW QUANTILES 					  #
  ###############################################################

  for (i in 1:nsamp) px[i]=  replace(px[i],which(px[i]<=10^(-50)),10^(-50))

  peso=delta/px
  tau= peso/sum(peso)

  lprobs=length(probs)

  IPW.quantile=rep(NA,length=lprobs)

  for (i in 1:lprobs){
    IPW.quantile[i]=weighted.fractile(y, tau, p=probs[i])
  }

  res=list(px=px, IPW.quantile=IPW.quantile)
  return(res)
}



