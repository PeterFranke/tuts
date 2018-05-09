#' Time-uncertain AR(1) model
#'
#' \code{tuar1} estimates unbiased parameters of AR(1) model.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates of timing of observations.
#' @param ti.sd A vector of standard deviations of timing.
#' @param n.sim A number of simulations.
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
#' # Run the tuar1 model fitting with cross validation:
#' n.sim=1000
#' TUAR1=tuar1(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim, CV=TRUE)
#'
#' # Generate summary results (optional parameters are listed in brackets):
#' summary(TUAR1)                                # Summary tables (CI, burn).
#'
#' # Plots and diagnostics (optional parameters are listed in brackets):
#' plot(TUAR1,type='predTUTS')                   # One step out of salmple predictions of the model (CI, burn).
#' plot(TUAR1,type='par', burn=0.4)              # Distributions of parameters of the AR(1) model (burn).
#' plot(TUAR1,type='mcmc')                       # mcmc diagnostics.
#' plot(TUAR1,type='cv', burn=0.4, CI=0.9)       # 5 fold cross validation plot(CI, burn).
#' plot(TUAR1,type='GR')                         # Gelman-Rubin diagnostic (CI, burn).
#' plot(TUAR1,type='volatility')                 # Volatility realizaitons.
#'
#' @export
#'
tuar1=function(y,ti.mu,ti.sd,n.sim,CV=FALSE, ...){
# Data checking and basic operations ---------------------------------------------------------------------
if (length(y)*2!=length(ti.mu)+length(ti.sd)){stop("Verify the input data.")}
if(is.numeric(y)==FALSE ){stop("y must be a vector of rational numbers.")}
#if(is.numeric(ti.mu)==FALSE | sum((ti.mu)<0)>0 ){
 # stop("ti.mu must be a vector of positive rational numbers.")}
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
modelstring="model {
  # Likelihood  ----------------------------------------------------------------------------------------
  for (i in 2:n) {
    y[i]~ dnorm(const+ alpha1*y[i-1]*(ti.sim[i]-ti.sim[i-1]),precision/sqrt(ti.sim[i]-ti.sim[i-1]))
  }
  for (i in 1:n) {
    ti.sim.tmp[i]~ dnorm(ti.mu[i],ti.prec[i])
  }
  ti.sim<-sort(ti.sim.tmp)

  # Priors   ----------------------------------------------------------------------------------------------
  const~ dnorm(0,1)
  alpha1~ dnorm(0,1)
  precision~dgamma(1.0E-3, 1.0E-3)
}"

# R2Jags Main Sim  ---------------------------------------------------------------------------------------
data=list(y=y,ti.mu=ti.mu,ti.prec=1/ti.sd^2,n=length(ti.mu))
for(k in (1:n.chains)){
  inits = parallel.seeds("base::BaseRNG", n.chains)
}
model=jags.model(textConnection(modelstring), data=data,inits=inits, n.chains=n.chains)
update(model,n.iter=n.sim,thin=Thin)
output=coda.samples(model=model,variable.names=c("const","alpha1","precision","ti.sim"), n.iter=n.sim, thin=Thin)
DIC = dic.samples(model=model,n.iter=n.sim,thin=Thin)
Sim.Objects=JAGS.objects(output)
Sim.Objects$JAGS=output
Sim.Objects$DIC=DIC
Sim.Objects$y=y
Sim.Objects$ti.mu=ti.mu
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
class(Sim.Objects)='tuts_ar1'
return(Sim.Objects)
}

#' Prints summary of tuts_ar1 objects
#'
#' \code{summary.tuts_ar1(x,...)} Prints summary of a tuts_ar1 object.
#'
#' @param x A tuts_ar1 object.
#' @param ... list of optional parameters. The list contains burn-in parameter
#' ranging from 0 to 0.5, with the default value burn=0, and the credible interval parameter
#' ranging between 0.5 and 1, with the default CI=0.99.
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Run the tuar1 model fitting with cross validation:
#' n.sim=1000
#' TUAR1=tuar1(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim, CV=TRUE)
#'
#' # Generate summary results (optional parameters are listed in brackets):
#' summary(TUAR1)                                # Summary tables (CI, burn).
#'
#' @export
summary.tuts_ar1 = function(x, ...) {
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
  n.sim=dim(x$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}
  # ----------------------------------------------------------------------------
  cat('\n')
  cat('Estimates of parameters of interest and timing:\n')
  cat('-----------------------------------------------\n')

  const=x$const[BURN:dim(x$const)[1]]
  const.lwr=quantile(const,(1-CI)/2)
  const.med=quantile(const,0.5)
  const.upr=quantile(const,1-(1-CI)/2)

  alpha1=x$alpha1[BURN:dim(x$alpha1)[1]]
  alpha1.lwr=quantile(alpha1,(1-CI)/2)
  alpha1.med=quantile(alpha1,0.5)
  alpha1.upr=quantile(alpha1,1-(1-CI)/2)

  precision=x$precision[BURN:dim(x$precision)[1]]
  precision.lwr=quantile(precision,(1-CI)/2)
  precision.med=quantile(precision,0.5)
  precision.upr=quantile(precision,1-(1-CI)/2)

  ti=x$ti.sim[BURN:dim(x$ti.sim)[1],]
  ti.lwr=apply(ti,2,'quantile',(1-CI)/2)
  ti.med=apply(ti,2,'quantile',0.5)
  ti.upr=apply(ti,2,'quantile',1-(1-CI)/2)
  tiNames=names(ti.med)

  LWR=c(const.lwr,alpha1.lwr,precision.lwr,ti.lwr)
  MED=c(const.med,alpha1.med,precision.med,ti.med)
  UPR=c(const.upr,alpha1.upr,precision.upr,ti.upr)
  TABLE2=data.frame(LWR,MED,UPR)
  row.names(TABLE2)=c('const','alpha1','precision',tiNames)

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
#' Graphical summaries and diagnostics of tuts_ar1' objects
#'
#' \code{plot.tuts_ar1(x,type,...)} plots summaries and diagnostics of tuts_ar1 objects.
#'
#' @param x A tuts_ar1 objects.
#'
#' @param type plot/diagnostic type (options:
#' 'predTUTS' plots one step predictions of the model, 'GR' plots Gelman-Rubin diagnostics,
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
#' # Run the tuar1 model fitting with cross validation:
#' n.sim=1000
#' TUAR1=tuar1(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim, CV=TRUE)
#'
#' # Plots and diagnostics (optional parameters are listed in brackets):
#' plot(TUAR1,type='predTUTS')                   # One step out of salmple predictions of the model (CI, burn).
#' plot(TUAR1,type='par', burn=0.4)              # Distributions of parameters of the AR(1) model (burn).
#' plot(TUAR1,type='mcmc')                       # mcmc diagnostics.
#' plot(TUAR1,type='cv', burn=0.4, CI=0.9)       # 5 fold cross validation plot(CI, burn).
#' plot(TUAR1,type='GR')                         # Gelman-Rubin diagnostic (CI, burn).
#' plot(TUAR1,type='volatility')                 # Volatility realizaitons.
#'
#' @export
plot.tuts_ar1 = function(x, type, ...) {
  if (sum(type==c('predTUTS','GR','cv','mcmc','par','volatility'))==0){
    stop('type should be set as either par, predTUTS, GR, cv, mcmc or volatility')
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
  n.sim=dim(x$alpha1)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}

  #############################################################
  if(type=='par') {
    par(mfrow=c(1,3))
    plot(density(x$const[BURN:dim(x$const)[1]]),main="const")
    plot(density(x$alpha1[BURN:dim(x$alpha1)[1]]),main="alpha1")
    plot(density(x$precision[BURN:dim(x$precision)[1]]),main="precision")
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

  ##########################################################################
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
    #text(x=max(x$y[2:(length(x$y)-2)])-2,y=min(PRED[2:(length(x$y)-2)])+1,LAB)
  }

  #################################################################################
  if(type=='GR') {
    if(burn>0){ABURNIN=TRUE} else{ABURNIN=FALSE}

    GELMAN=gelman.diag(x$JAGS,confidence =CI,autoburnin=ABURNIN)
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
    mcmcplot(x$JAGS, parms=c('const','alpha1','precision','ti.sim'))
  }
  if(type=='volatility') {
    plot(sqrt(sqrt(1/sqrt(x$precision[BURN:dim(x$precision)[1],]))),type='l', xlab='Sim ID',ylab='Std Deviation',main='Standard Deviaiton')
  }
}

