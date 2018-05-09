#' Time-uncertain time-persistence AR(1) model
#'
#' \code{tuar1redf(y,ti.mu,ti.sd,n.sim,CV=FALSE,...)} estimates parameters of AR(1) model using specification described in
#'  \emph{"Climate Time Series Analysis"} by M.Mudelsee. We modify the model to account for time-uncertainty.
#'
#'  Model estimates the autocorrelation parameter with a certain bias, which is still present after correction.
#'  The model is also not suitable for estimation of negative values of autocorrelation.
#'
#'  We include this function for further validation of applicability of the REDFIT model, which is often used
#'  in spectral analysis of unevenly spaced data.
#'
#'  Note: the function tuar1 generates unbiased estimates, and is not limited to positive values
#'   unlike the tuar1redf.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates of timing of observations.
#' @param ti.sd A vector of standard deviations of timing.
#' @param n.sim A number of simulations.
#' @param n.chains A number of chains.
#' @param CV cross-validation indicator.
#' @param ... A list of optional parameters. The list contains thinning parameter,
#' with the default value set to Thin=4, and the number of mcmc chains with the default value set
#'  to n.chains=2.
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Run the tuar1redf model with cross validation:
#' n.sim=1000; n.chains=2
#' AR1REDF=tuar1redf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=TRUE)
#'
#' # Generate summary results (optional parameters are listed in brackets):
#' summary(AR1REDF)                                # Summary statistics (urn, CI).
#'
#' # Plots and diagnostics (optional parameters are listed in brackets):
#' plot(AR1REDF,type='predTUTS',burn=0.2,CI=0.99)  # One step out of salmple predictions of the model (CI, burn).
#' plot(AR1REDF,type='par', burn=0.4)              # Distributions of parameters; (burn).
#' plot(AR1REDF,type='mcmc')                       # mcmc diagnostics.
#' plot(AR1REDF,type='cv', burn=0.4)               # 5 fold cross validation plot (CI, burn).
#' plot(AR1REDF,type='GR', burn=0.4)               # Gelman-Rubin diagnostic (CI, burn).
#' plot(AR1REDF,type='tau')                        # Volatility realizaitons.
#'
#' @export
tuar1redf=function(y,ti.mu,ti.sd,n.sim,CV=FALSE,...){
# Data checking and basic operations ---------------------------------------------------------------------
if (length(y)*2!=length(ti.mu)+length(ti.sd)){stop("Verify the input data.")}
if(is.numeric(y)==FALSE ){stop("y must be a vector of rational numbers.")}
if(is.numeric(ti.mu)==FALSE | sum((ti.mu)<0)>0 ){
  stop("ti.mu must be a vector of positive rational numbers.")}
if(is.numeric(ti.sd)==FALSE | sum((ti.sd)<0)>0 ){
  stop("ti.sd must be a vector of positive rational numbers.")}
if (sum(is.na(c(y,ti.mu,ti.sd)))>0){stop("Remove NAs.")}
if (n.sim!=abs(round(n.sim))){stop("n.sim must be a positive integer.")}
if (!is.logical(CV)){stop("CV must be a logical value.")}
  # Optional parameters
  dots = list(...)
  if(missing(...)){Thin=4; n.chains=2}
  if(!is.numeric(dots$Thin)){
    Thin=4
  } else{
    Thin=round(abs(dots$Thin))
  }

  if(!is.numeric(dots$n.chains)){
    n.chains=2
  } else{
    n.chains=round(abs(dots$n.chains))
  }
y=y[order(ti.mu,decreasing = FALSE)]; ti.sd=ti.sd[order(ti.mu,decreasing = FALSE)]
ti.mu=ti.mu[order(ti.mu,decreasing = FALSE)]
# JAGS model --------------------------------------------------------------------------------------------
modelstring="model {
for (i in 2:n) {
  y[i]~ dnorm(y[i-1] * exp(-(ti.sim[i]-ti.sim[i-1])/tau),1/(1-exp(-2*(ti.sim[i]-ti.sim[i-1])/tau )))
}
for (i in 1:n) {
  ti.sim.tmp[i]~ dnorm(ti.mu[i],ti.prec[i])
}
ti.sim<-sort(ti.sim.tmp)

for (i in 1:(n-1)){
  dt.sim[i]<- ti.sim[i+1]-ti.sim[i]
}
ar1<-exp(-mean(dt.sim)/tau)
ar1adj<-ar1+(1+3*ar1)/(n-1)
tau~dlnorm(0, 0.001)
}"

# R2Jags Main Sim  ---------------------------------------------------------------------------------------
data=list(y=y,ti.mu=ti.mu,ti.prec=1/ti.sd^2,n=length(ti.mu))
for(k in (1:n.chains)){
  inits = parallel.seeds("base::BaseRNG", n.chains)
}
model=jags.model(textConnection(modelstring), data=data,inits=inits, n.chains=n.chains)
update(model,n.iter=n.sim,thin=Thin)
output=coda.samples(model=model,variable.names=c("ar1","ar1adj","tau","ti.sim"), n.iter=n.sim, thin=Thin)
DIC = dic.samples(model=model,n.iter=n.sim,thin=Thin)

Sim.Objects=JAGS.objects(output)
Sim.Objects$JAGS=output
Sim.Objects$DIC=DIC

# Cross Validation -------------------------------------------------------------------------------------------------
if(CV==TRUE){
  print(noquote('Cross-validation of the model....'))
  folds = 5
  fold= sample(rep(1:folds,length=length(y)))
  for (i in 2:length(fold)){
    if (fold[i-1]==fold[i]){
      Sample=c(1:5)
      Sample=Sample[Sample !=fold[i]]
      fold[i]=sample(Sample,size=1)
    }
  }
  TI.SIM=apply(Sim.Objects$ti.sim,2,'quantile',0.5)
  Cores=min(c(detectCores()-1,folds));  cl = makeCluster(Cores); registerDoParallel(cl)

  BSFCV=function(i,y,ti.mu,ti.sd,modelstring,n.sim,fold){
    Y=y; Y[fold==i]=NA; Y[1]=c(y[1])
    data=list(y=Y, ti.mu=TI.SIM,ti.prec=1/ti.sd^2, n=length(ti.mu))
    model=jags.model(textConnection(modelstring), data=data,n.chains=1)
    update(model,n.iter=n.sim,thin=Thin)
    output=coda.samples(model=model,variable.names=c("y"), n.iter=n.sim, thin=Thin)
    return(output)
  }

  CVRES=foreach(i=1:folds,.export=c('jags.model','coda.samples')) %dopar%
    BSFCV(i,fold=fold,y=y,ti.mu=TI.SIM,ti.sd=ti.sd, modelstring=modelstring, n.sim=n.sim)

  stopCluster(cl)

  pred_y = array(NA, dim=c(dim(JAGS.objects(CVRES[[1]])$y)[1], length(ti.mu)))
  colnames(pred_y)= paste("y[",1:length(y),"]",sep="")

  for (i in 1:folds){
    pred_y[,fold==i] = JAGS.objects(CVRES[[i]])$y[,fold==(i)]
  }

  Sim.Objects$CVpred=pred_y[,2:dim(pred_y)[2]]
}
Sim.Objects$y=y
Sim.Objects$ti.mu=ti.mu
class(Sim.Objects)='tuts_ar1redf'
return(Sim.Objects)
}
#' Summaries of tuts_ar1redf objects
#'
#' \code{summary.tuts_ar1redf(x,...)} Prints summary of a tuts_ar1redf object.
#'
#' @param object A tuts_ar1redf object.
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Run the tuar1redf model with cross validation:
#' n.sim=1000; n.chains=2
#' AR1REDF=tuar1redf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=TRUE)
#'
#' # Generate summary results (optional parameters are listed in brackets):
#' summary(AR1REDF)                                # Summary statistics (urn, CI).
#' @export
#'
summary.tuts_ar1redf = function(x, ...) {
  dots = list(...)
  if(missing(...)){burn=0; CI=0.99}

  if(!is.numeric(dots$CI)){
    CI=0.99
  } else{
    if(dots$CI<=0.5 | dots$CI> 1){stop('Credible interval is bounded between 0.5 and 1')}
    CI=dots$CI
  }

  if(!is.numeric(dots$burn)){
    burn=0
  } else{
    burn=dots$burn
    if(burn<0 | burn> 0.5){stop('burn is bounded between 0 and 0.5')
    }
  }
  n.sim=dim(x$ar1)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}
  # ----------------------------------------------------------------------------
  cat('\n')
  cat('Estimates of parameters of interest and timing:\n')
  cat('-----------------------------------------------\n')
  #"ar1","ar1adj","tau","ti.sim"
  ar1=x$ar1[BURN:dim(x$ar1)[1]]
  ar1.lwr=quantile(ar1,(1-CI)/2)
  ar1.med=quantile(ar1,0.5)
  ar1.upr=quantile(ar1,1-(1-CI)/2)

  ar1adj=x$ar1adj[BURN:dim(x$ar1adj)[1]]
  ar1adj.lwr=quantile(ar1adj,(1-CI)/2)
  ar1adj.med=quantile(ar1adj,0.5)
  ar1adj.upr=quantile(ar1adj,1-(1-CI)/2)

  tau=x$tau[BURN:dim(x$tau)[1]]
  tau.lwr=quantile(tau,(1-CI)/2)
  tau.med=quantile(tau,0.5)
  tau.upr=quantile(tau,1-(1-CI)/2)

  ti=x$ti.sim[BURN:dim(x$ti.sim)[1],]
  ti.lwr=apply(ti,2,'quantile',(1-CI)/2)
  ti.med=apply(ti,2,'quantile',0.5)
  ti.upr=apply(ti,2,'quantile',1-(1-CI)/2)
  tiNames=names(ti.med)

  LWR=c(ar1.lwr,ar1adj.lwr,tau.lwr,ti.lwr)
  MED=c(ar1.med,ar1adj.med,tau.med,ti.med)
  UPR=c(ar1.upr,ar1adj.upr,tau.upr,ti.upr)
  TABLE2=data.frame(LWR,MED,UPR)
  row.names(TABLE2)=c('ar1','ar1adj','tau',tiNames)

  colnames(TABLE2)=c(paste(round((1-CI)/2,3)*100,"%",sep=""),'50%',paste(round(1-(1-CI)/2,3)*100,"%",sep=""))
  print(TABLE2)
  # ----------------------------------------------------------------------------
  cat('\n')
  cat('Deviance information criterion:\n')
  cat('-------------------------------\n')
  print(x$DIC)
  cat('-------------------------------\n')
}
####################################################################################################
#' Graphical summaries and diagnostics of tuts_ar1redf objects
#'
#' \code{plot.tuts_ar1redf(x,type,...)} plots summaries and diagnostics of tuts_ar1redf objects.
#'
#' @param x A tuts_tuar1 objects.
#'
#' @param type plot/disgnostics type (options:
#' 'predTUTS' plots one step preditions of the model, 'GR' plots Gelman-Rubin diagnostics,
#' 'cv' plots 5-fold cross validation, 'mcmc' plots diagnostics of mcmc objects)
#' @param ... list of optional parameters: 'burn' (burn-in parameter ranging from 0 to 0.7 with
#'  default value set to 0), and CI (credible interval ranging from 0.3 to 1 with
#'  default value set to 0.95)
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Run the tuar1redf model with cross validation:
#' n.sim=1000; n.chains=2
#' AR1REDF=tuar1redf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=TRUE)
#'
#' # Plots and diagnostics (optional parameters are listed in brackets):
#' plot(AR1REDF,type='predTUTS',burn=0.2,CI=0.99)  # One step out of salmple predictions of the model (CI, burn).
#' plot(AR1REDF,type='par', burn=0.4)              # Distributions of parameters; (burn).
#' plot(AR1REDF,type='mcmc')                       # mcmc diagnostics.
#' plot(AR1REDF,type='cv', burn=0.4)               # 5 fold cross validation plot (CI, burn).
#' plot(AR1REDF,type='GR', burn=0.4)               # Gelman-Rubin diagnostic (CI, burn).
#' plot(AR1REDF,type='tau')                        # Volatility realizaitons.
#'
#' @export
plot.tuts_ar1redf = function(x, type, ...) {
  if (sum(type==c('predTUTS','GR','cv','mcmc','par','tau'))==0){
    stop('type should be set as either par, predTUTS, GR, cv, mcmc or tau')
  }

  dots = list(...)
  if(missing(...)){burn=0; CI=0.99}

  if(!is.numeric(dots$CI)){
    CI=0.99
  } else{
    if(dots$CI<=0.5 | dots$CI> 1){stop('Credible interval is bounded between 0.5 and 1')}
    CI=dots$CI
  }

  if(!is.numeric(dots$burn)){
    burn=0
  } else{
    burn=dots$burn
    if(burn<0 | burn>0.7){stop('burn is bounded between 0 and 0.7')
    }
  }
  n.sim=dim(x$ar1)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}


  par(mfrow=c(1,1))
  #############################################################
  if(type=='par') {

    par(mfrow=c(1,3))
    plot(density(x$ar1[BURN:dim(x$ar1)[1]]),main="AR(1)")
    plot(density(x$ar1adj[BURN:dim(x$ar1adj)[1]]),main="AR(1) adjusted")
    plot(density(x$tau[BURN:dim(x$tau)[1]]),main="tau")
    par(mfrow=c(1,1))
  }
  ##########################################################################
  if(type=='predTUTS') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED.LWR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',(1-CI)/2)
    PRED.MED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    PRED.UPR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',1-(1-CI)/2)

    ti.sim=apply(x$ti.sim[BURN:dim(x$ti.sim)[1],],2,'quantile',0.5)

    MAIN=paste("One step out of sample predictions at CI= ", CI*100,"%",sep='')
    plot(y=x$y,x=x$ti.mu,type='l',main=MAIN,ylab="Observations",xlab='time',
         ylim=c(min(x$CVpred),1.2*max(x$CVpred)), xlim=c(min(x$ti.mu,ti.sim),
                                                         max(x$ti.mu,ti.sim)),lwd=2)
    lines(y=PRED.LWR,x=ti.sim[2:length(ti.sim)],type='l',col='blue',lwd=1,lty=2)
    lines(y=PRED.MED,x=ti.sim[2:length(ti.sim)],type='l',col='blue',lwd=1,lty=1)
    lines(y=PRED.UPR,x=ti.sim[2:length(ti.sim)],type='l',col='blue',lwd=1,lty=2)

    legend("topright",legend = c("Observed","Upper CI","Medium","Lower CI"),
           col=c("black","blue","blue","blue"),lwd=c(2,1,1,1),lty=c(1,2,1,2))
  }
  #################################################################################
  if(type=='cv') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    MAIN="Cross-Validation: One step out of sample predictions"

    par(mfrow=c(1,1))
    plot(x=x$y[2:length(x$y)],y=PRED,xlab="Actual",ylab="Predicted", main=MAIN, pch=18)
    abline(0,1,col='blue')

    RSQ=cor(x$y[2:length(x$y)],PRED)^2* 100

    LAB = bquote(italic(R)^2 == .(paste(format(RSQ, digits = 0),"%",sep="")))
    text(x=(min(x$y)+0.9*(max(x$y)-min(x$y))),y=(min(PRED)+0.1*(max(PRED)-min(PRED))),LAB)
  }

  #################################################################################
  if(type=='GR') {
    if(BURN>1){ABURNIN=TRUE} else{ABURNIN=FALSE}

    GELMAN=gelman.diag(x$JAGS, confidence =CI,autoburnin=ABURNIN,multivariate = FALSE)
    TIOBJ=rownames(GELMAN$psrf)%in% c(paste("ti.sim[",1:length(x$ti.mu),"]",sep=""))

    REG=GELMAN$psrf[!TIOBJ,]
    labels =rownames(REG)
    TI=GELMAN$psrf[TIOBJ,]
    par(mfrow=c(2,1))
    plot(REG[,1],ylim=c(0,max(REG)+1.5),xaxt='n',ylab="Factor", xlab='Parameters',
         main=paste("Gelman-Rubin diagnostics: Potential scale reduction factors \n with the upper confidence bounds at ",
                    CI*100,"%",sep="")
    )
    text(seq(1,length(labels),by=1), par("usr")[3]-max(REG)*0.2, srt = 00, adj= 0.5, xpd = TRUE, labels =labels, cex=0.65)
    arrows(x0=1:dim(REG)[1],y0=REG[,1],x1=1:dim(REG)[1],y1=REG[,2],code=3,length=0.04,angle=90,col='darkgray')
    abline(h=1);
    legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")

    plot(TI[,1],ylim=c(0,(max(TI)*1.8)),ylab="Factor", xlab='Simulated timings of observations')
    arrows(x0=1:dim(TI)[1],y0=TI[,1],x1=1:dim(TI)[1],y1=TI[,2],code=3,length=0.04,angle=90,col='darkgray')
    abline(h=1);
    legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")

    par(mfrow=c(1,1))
  }
  if(type=='mcmc') {
    mcmcplot(x$JAGS, parms=c('ar1','ar1adj','tau','ti.sim'))
  }
  if(type=='tau') {
    plot(sqrt(sqrt(1/sqrt(x$tau[BURN:dim(x$tau)[1],]))),type='l', xlab='Sim ID',ylab='tau',main='Persistence of time')
  }
}

