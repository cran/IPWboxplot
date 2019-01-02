library(isotone)

################################################################################################
# FUNCTION USED TO DRAW THE BOXPLOTS ADAPTED TO MISSING VALUES   AND SKEWED DISTRIBUTIONS      #
################################################################################################
IPW.ASYM.boxplot=function(y,px=NULL,x=NULL,graph=c("IPW","both"),names=c("IPW Asymmetric Boxplot", "NAIVE Asymmetric Boxplot"), size.letter=1.2,
                          method=c("quartile","octile"), ctea=-4, cteb=3, lim.inf=NULL,lim.sup=NULL,main=" ",xlab = " ", ylab ="  ",color="black")
{
  ############################################################################################################################################
  # Arguments
  #
  #   y			Required. Numerical vector of values with possible missing values codified NA or NAN with length n.
  #   px		Optional. Numerical vector of probabilities. If not provided a logistic fit is performed using x
  #   x 		Optional. The matrix of fully observed variables used to estimate the missing model with dimension nrows=n and ncol=dimension.
  #   graph   	Optional. Character string indicating if the plot contains two boxplots ("both") or
  #			only the boxplot computed with the inversely probability weighted quantiles("IPW").
  #			The default is "IPW".
  #   names		Optional. Character string to name the boxplots.
  #			The default is "IPW Asymmetric Boxplot" when graph="IPW" and
  #			c("IPW Asymmetric Boxplot", "NAIVE Asymmetric Boxplot") when graph= "both"
  #   size.letter	Optional. The font size of names.
  #   method		Optional. Character string indicating if the measure of asymmetry is based on the quartiles ("quartile") or the octiles ("octile").
  #			The default is "quartile".
  #   ctea, cteb	Optional. Scaling factors multiplied by the asymmetry measure to determine outlier boundary.
  #			When ctea=cteb=0 the IPW boxplot is obtained. ctea is a negative value while cteb is positive.
  #			The default ones correspond to the choices in Hubert and Vandervieren (2008) when using the medcouple.
  #   lim.inf		Optional. The lower limit of the plot if supplied by the user.
  #   lim.sup		Optional. The upper limit of the plot if supplied by the user.
  #   main		Optional. Character string to title the plot. The default is "IPW Boxplot".
  #   xlab		Optional. Character string to indicate the label of the horizontal axis.
  #   ylab		Optional. Character string to indicate the label of the vertical axis.
  #   color   	Optional. Color for the IPW Boxplot.
  ############################################################################################################################################

  ########################################################################
  # Value
  #
  ############################################################################################################################################
  # Value
  #
  #  px 			Numerical vector of probabilities.
  #  IPW.Quartiles  	Numerical vector of inversely probability weighted quartiles.
  #  IPW.whisker    	Numerical vector of lower and upper whisker calculated from IPW quantiles.
  #  out.IPW        	Numerical vector of data points detected as atypical by the IPW boxplot.
  #  SKEW.IPW	  	Skewness measure based on the IPW quartiles (method="quartile") or IPW octiles (method="octile").
  #  NAIVE.Quartiles	Numerical vector of naive quartiles computed from the subset of non-missing values of y. Returned only when graph="both".
  #  NAIVE.whisker  	Numerical vector of lower and upper whisker obtained from the Naive quantiles. Returned only when graph="both".
  #  out.NAIVE      	Numerical vector of data points detected as atypical by the Naive boxplot. Returned only when graph="both".
  #  SKEW.NAIVE	  	Skewness measure based on the Naive quartiles (method="quartile").
  #				or IPW octiles (method="octile"), computed from the subset of non-missing values of y.
  #				Returned only when graph="both".
  ########################################################################

  #########################################################################
  #---Preliminary checks---
  #########################################################################

  dimension=NCOL(x)

  if (is.null(x)=="TRUE" & is.null(px)=="TRUE")
    stop("ERROR: It is neccesary to supply the vector of dropout probabilities for each observation or a covariate to estimate it")

  if (is.null(px)=="FALSE")
  {
    if (min(px)<=0)
      stop("ERROR: px should take positive values")
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

  GRAPH=graph[1]
  COLOR=color[1]
  nsamp=length(y)
  METHOD=method[1]

  delta=rep(1,nsamp)
  for (i in 1: nsamp){
    if (is.na(y[i])=="TRUE"|is.nan(y[i])=="TRUE")
    {
      delta[i]<-0
      y[i]=NA
    }
  }


  if(is.null(px)=="FALSE"){PROBS="The dropout probability is given"}

  if(is.null(px)){
    ###############################################################
    #  Estimation of the dropout probability for each observation #
    ###############################################################
    a=glm(delta~x,family="binomial")
    px=a$fitted.values
    PROBS="LOGISTIC"
  }

  ###############################################################
  #  Estimation of the IPW QUANTILES AND OUTLIERS DETECTION     #
  ###############################################################

  for (i in 1:nsamp) px[i]=  replace(px[i],which(px[i]<=10^(-50)),10^(-50))

  peso=delta/px
  tau= peso/sum(peso)
  yp <- y[ tmp <- (!is.na(y))]

  mediana.IPW=weighted.fractile(y, tau, p=0.5)
  cuantil025.IPW= weighted.fractile(y, tau, 0.25)
  cuantil075.IPW= weighted.fractile(y, tau, 0.75)
  theta.IPW=c(cuantil025.IPW,mediana.IPW,cuantil075.IPW)
  IQR.IPW=theta.IPW[3]-theta.IPW[1]

  if(METHOD=="quartile"){
    deno=cuantil075.IPW-cuantil025.IPW
    if(deno==0)
      print("WARNING: Too many ties. The quartiles are equal.")
    deno=  replace(deno,which(deno<=10^(-50)),10^(-50))
    num=(cuantil075.IPW-mediana.IPW)-(mediana.IPW-cuantil025.IPW)
    SKEW.IPW=num/deno
  }else
  {
    cuantil0125.IPW= weighted.fractile(y, tau, 0.125)
    cuantil0875.IPW= weighted.fractile(y, tau, 0.875)
    deno=cuantil0875.IPW-cuantil0125.IPW
    if(deno==0)
      print("WARNING: Too many ties. The octiles are equal.")
    deno=  replace(deno,which(deno<=10^(-50)),10^(-50))
    num=(cuantil0875.IPW-mediana.IPW)-(mediana.IPW-cuantil0125.IPW)
    SKEW.IPW=num/deno
  }

  cteinf.IPW=ctea*(SKEW.IPW>0) - cteb*(SKEW.IPW<0)
  ctesup.IPW= cteb*(SKEW.IPW>0) - ctea*(SKEW.IPW<0)

  bigote.IPW=c(theta.IPW[1]-1.5*exp(cteinf.IPW*SKEW.IPW)*IQR.IPW,theta.IPW[3]+1.5*exp(ctesup.IPW*SKEW.IPW)*IQR.IPW)

  bigo.IPW=bigote.IPW
  bigo.IPW[1]=min(yp[yp>=bigote.IPW[1]])
  bigo.IPW[2]=max(yp[yp<=bigote.IPW[2]])

  outsup.IPW=yp[yp>bigote.IPW[2]]
  outinf.IPW=yp[yp<bigote.IPW[1]]
  out.IPW=c(outinf.IPW,outsup.IPW)

  ###############################################################
  #  Estimation of the NAIVE QUANTILES AND OUTLIERS DETECTION   #
  ###############################################################

  largo=length(yp)
  tau.NAIVE=rep(1,length=largo)/largo
  mediana.NAIVE=weighted.fractile(yp, tau.NAIVE, p=0.5)
  cuantil025.NAIVE= weighted.fractile(yp, tau.NAIVE, 0.25)
  cuantil075.NAIVE= weighted.fractile(yp, tau.NAIVE, 0.75)
  theta.NAIVE=c(cuantil025.NAIVE,mediana.NAIVE,cuantil075.NAIVE)
  IQR.NAIVE=theta.NAIVE[3]-theta.NAIVE[1]

  if(METHOD=="quartile"){
    deno=cuantil075.NAIVE-cuantil025.NAIVE
    if(deno==0)
      print("WARNING: Too many ties. The quartiles are equal.")
    deno=  replace(deno,which(deno<=10^(-50)),10^(-50))
    num=(cuantil075.NAIVE-mediana.NAIVE)-(mediana.NAIVE-cuantil025.NAIVE)
    SKEW.NAIVE=num/deno
  }else
  {
    cuantil0125.NAIVE= weighted.fractile(yp, tau.NAIVE, 0.125)
    cuantil0875.NAIVE= weighted.fractile(yp, tau.NAIVE, 0.875)
    deno=cuantil0875.NAIVE-cuantil0125.NAIVE
    if(deno==0)
      print("WARNING: Too many ties. The octiles are equal.")
    deno=  replace(deno,which(deno<=10^(-50)),10^(-50))
    num=(cuantil0875.NAIVE-mediana.NAIVE)-(mediana.NAIVE-cuantil0125.NAIVE)
    SKEW.NAIVE=num/deno
  }

  cteinf.NAIVE=ctea*(SKEW.NAIVE>0) - cteb*(SKEW.NAIVE<0)
  ctesup.NAIVE=cteb*(SKEW.NAIVE>0) - ctea*(SKEW.NAIVE<0)

  bigote.NAIVE=c(theta.NAIVE[1]-1.5*exp(cteinf.NAIVE*SKEW.NAIVE)*IQR.NAIVE,theta.NAIVE[3]+1.5*exp(ctesup.NAIVE*SKEW.NAIVE)*IQR.NAIVE)

  bigo.NAIVE=bigote.NAIVE
  bigo.NAIVE[1]=min(yp[yp>=bigote.NAIVE[1]])
  bigo.NAIVE[2]=max(yp[yp<=bigote.NAIVE[2]])

  outsup.NAIVE=yp[yp>bigote.NAIVE[2]]
  outinf.NAIVE=yp[yp<bigote.NAIVE[1]]
  out.NAIVE=c(outinf.NAIVE,outsup.NAIVE)

  ####################################################################################
  #---Now we draw the boxes---
  ####################################################################################


  if(is.null(lim.inf)=="TRUE") lim.inf=min(yp)
  if(is.null(lim.sup)=="TRUE") lim.sup=max(yp)
  if(min(yp)<lim.inf|lim.sup<max(yp))
    print("WARNING: The box is truncated")

  if (GRAPH=="both")
  {
    ###################################################################################
    # We draw the NAIVE and IPW boxplots in the same graph
    ###################################################################################

    plot(c(-0.1,3.9),c(lim.inf,lim.sup), main = main, xlab = xlab,ylab = ylab, xaxt = "n" ,pch=21,col="white",bg="white")

    axis(1, at=c(1,3), labels=names,las=1,cex.axis=size.letter)

    rango=0.5
    punto=c(1-rango,1+rango)
    lines( punto,c(theta.IPW[3],theta.IPW[3]),col=COLOR,lwd=1)
    lines( punto,c(theta.IPW[1],theta.IPW[1]),col=COLOR,lwd=1)
    lines( punto,c(theta.IPW[2],theta.IPW[2]),col=COLOR,lwd=3)
    lines( c(punto[1],punto[1]),c(theta.IPW[1],theta.IPW[3]),col=COLOR,lwd=1)
    lines( c(punto[2],punto[2]),c(theta.IPW[1],theta.IPW[3]),col=COLOR,lwd=1)

    rango2=0.3
    punto2=c(1-rango2,1+rango2)
    lines( punto2,c(bigo.IPW[1],bigo.IPW[1]),col=COLOR,lwd=1)
    lines(punto2,c(bigo.IPW[2],bigo.IPW[2]),col=COLOR,lwd=1)
    lines(c(1,1),c(theta.IPW[3],bigo.IPW[2]),col=COLOR,lwd=1)
    lines(c(1,1),c(theta.IPW[1],bigo.IPW[1]),col=COLOR,lwd=1)

    eje=rep(1,length=length(outsup.IPW))
    points(eje,outsup.IPW,col=COLOR,pch=1)
    eje=rep(1,length=length(outinf.IPW))
    points(eje,outinf.IPW,col=COLOR,pch=1)

    punto.NAIVE=c(3-rango,3+rango)
    lines( punto.NAIVE,c(theta.NAIVE[3],theta.NAIVE[3]),col=COLOR,lwd=1)
    lines( punto.NAIVE,c(theta.NAIVE[1],theta.NAIVE[1]),col=COLOR,lwd=1)
    lines( punto.NAIVE,c(theta.NAIVE[2],theta.NAIVE[2]),col=COLOR,lwd=3)
    lines( c(punto.NAIVE[1],punto.NAIVE[1]),c(theta.NAIVE[1],theta.NAIVE[3]),col=COLOR,lwd=1)
    lines( c(punto.NAIVE[2],punto.NAIVE[2]),c(theta.NAIVE[1],theta.NAIVE[3]),col=COLOR,lwd=1)

    punto.NAIVE2=c(3-rango2,3+rango2)
    lines( punto.NAIVE2,c(bigo.NAIVE[1],bigo.NAIVE[1]),col=COLOR,lwd=1)
    lines(punto.NAIVE2,c(bigo.NAIVE[2],bigo.NAIVE[2]),col=COLOR,lwd=1)
    lines(c(3,3),c(theta.NAIVE[3],bigo.NAIVE[2]),col=COLOR,lwd=1)
    lines(c(3,3),c(theta.NAIVE[1],bigo.NAIVE[1]),col=COLOR,lwd=1)

    eje=rep(3,length=length(outsup.NAIVE))
    points(eje,outsup.NAIVE,col=COLOR,pch=1)
    eje=rep(3,length=length(outinf.NAIVE))
    points(eje,outinf.NAIVE,col=COLOR,pch=1)

    cat("The method used to estimate the dropout probability is ",PROBS,"\n")

    cat("IPW Quartiles","\n","25%    50%   75%","\n",round(theta.IPW,4),"\n")
    cat("Naive Quartiles","\n","25%    50%   75%","\n",round(theta.NAIVE,4),"\n")

    cat("Lower and upper whiskers of the IPW Boxplot","\n","Lower   Upper","\n",round(bigo.IPW,4),"\n")
    cat("Lower and upper whiskers of the Naive Boxplot","\n","Lower   Upper","\n",round(bigo.NAIVE,4),"\n")

    cat("Skewness measure computed from the IPW ",METHOD,"\n",round(SKEW.IPW,4),"\n")
    cat("Skewness measure computed from the NAIVE ",METHOD,"\n",round(SKEW.NAIVE,4),"\n")

    res=list(px=px, IPW.Quartiles=theta.IPW,IPW.whisker=bigo.IPW,out.IPW=out.IPW,SKEW.IPW=SKEW.IPW,NAIVE.Quartiles=theta.NAIVE,NAIVE.whisker=bigo.NAIVE,out.NAIVE=out.NAIVE,SKEW.NAIVE=SKEW.NAIVE)

  }else{

    ###################################################################################
    # We draw the  IPW boxplot
    ###################################################################################

    plot(c(0,2),c(lim.inf,lim.sup), main = main, xlab = xlab,ylab = ylab, xaxt = "n" ,pch=21,col="white",bg="white")

    axis(1, at=c(1), labels=names[1],las=1,cex.axis=size.letter)

    rango=0.5
    punto=c(1-rango,1+rango)
    lines( punto,c(theta.IPW[3],theta.IPW[3]),col=COLOR,lwd=1)
    lines( punto,c(theta.IPW[1],theta.IPW[1]),col=COLOR,lwd=1)
    lines( punto,c(theta.IPW[2],theta.IPW[2]),col=COLOR,lwd=3)
    lines( c(punto[1],punto[1]),c(theta.IPW[1],theta.IPW[3]),col=COLOR,lwd=1)
    lines( c(punto[2],punto[2]),c(theta.IPW[1],theta.IPW[3]),col=COLOR,lwd=1)

    rango2=0.3
    punto2=c(1-rango2,1+rango2)
    lines( punto2,c(bigo.IPW[1],bigo.IPW[1]),col=COLOR,lwd=1)
    lines(punto2,c(bigo.IPW[2],bigo.IPW[2]),col=COLOR,lwd=1)
    lines(c(1,1),c(theta.IPW[3],bigo.IPW[2]),col=COLOR,lwd=1)
    lines(c(1,1),c(theta.IPW[1],bigo.IPW[1]),col=COLOR,lwd=1)

    eje=rep(1,length=length(outsup.IPW))
    points(eje,outsup.IPW,col=COLOR,pch=1)
    eje=rep(1,length=length(outinf.IPW))
    points(eje,outinf.IPW,col=COLOR,pch=1)

    cat("The method used to estimate the dropout probability is ",PROBS,"\n")

    cat("IPW Quartiles","\n","25%    50%   75%","\n",round(theta.IPW,4),"\n")

    cat("Lower and upper whiskers of the IPW Boxplot","\n","Lower   Upper","\n",round(bigo.IPW,4),"\n")

    cat("Skewness measure computed from the IPW ",METHOD,"\n",round(SKEW.IPW,4),"\n")

    res=list(px=px, IPW.Quartiles=theta.IPW,IPW.whisker=bigo.IPW,out.IPW=out.IPW,SKEW.IPW=SKEW.IPW)
  }

  return(res)
}
