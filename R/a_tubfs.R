#' Bayesian Frequency Selection of time-uncertain data sets
#'
#' \code{tubfs}
#' spectral analysis of time-uncertain time series using the Bayesian Frequency Selection method
#'  described in the paper \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2492}{Frequency selection in paleoclimate time series:
#'  A model-based approach incorporating possible time uncertainty} by P. Franke, Prof B. Huntley, Dr A. Parnell.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates/observed timings of observations.
#' @param ti.sd A vector of standard deviations of timings.
#' @param n.sim A number of simulations.
#' @param CV TRUE/FALSE cross-validation indicator.
#' @param ... optional arguments: \cr
#' - n.chains: number of MCMC chains, the default number of chains is set to 2.\cr
#' - Thin: thinning factor, the default values is set to 4.\cr
#' - m: maximum number of significant frequencies in the data, the default value is set to 5. \cr
#' - polyorder: the polynomial regression component, the default odrer is set to 3. \cr
#' - freqs: set to a positive integer k returns a vector of k equally spaced frequencies in the Nyquist
#'   range. freqs can be provided as a vector of custom frequencies of interest. Set to 'internal'
#'   (the default value) generates a vector of equally spaced frequencies in the Nyquist range.
#'
#' @examples
#' \donttest{
#' #1. Import or simulate the data (a simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' #2. Fit the model:
#' n.sim=1000
#' BFS=tubfs(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=2, CV=TRUE)
#'
#' #3. Generate summary results (optional parameters are listed in brackets):
#' summary(BFS)                               # Summary results (CI, burn).
#' summary(BFS,burn=0.2)                      # Results after 20% of burn-in (CI).
#'
#' #4. Generate plots and diagnostics of the model (optional parameters are listed in brackets):
#' plot(BFS,type='periodogram')               # spectral analysis (CI, burn).
#' plot(BFS,type='predTUTS', CI=0.99)         # One step predictions (CI, burn).
#' plot(BFS,type='cv')                        # 5 fold cross validation plot (CI, burn).
#' plot(BFS,type='GR')                        # Gelman-Rubin diagnostics (CI, burn).
#' plot(BFS,type='mcmc')                      # mcmc diagnostics.
#' plot(BFS,type='volatility')                # Volatility realizaitons.
#'}
#' @export
tubfs=function(y,ti.mu,ti.sd,n.sim,CV=FALSE,...){
# Data checking and basic operations
if (length(y)*4!=length(ti.mu)*2+length(ti.sd)*2){stop("Vectors y, ti.mu and ti.sd should be of equal lengths.")}
if(is.numeric(y)==FALSE ){stop("y must be a vector of rational numbers.")}
if(is.numeric(ti.sd)==FALSE | sum((ti.sd)<0)>0 ){
  stop("ti.sd must be a vector of positive rational numbers.")}
if (sum(is.na(c(y,ti.mu,ti.sd)))>0){stop("Remove NAs.")}
if (n.sim!=abs(round(n.sim))){stop("n.sim must be a positive integer.")}
if (!is.logical(CV)){stop("CV must be a logical value.")}

dots = list(...)
if(missing(...)){Thin=4; n.chains=2; polyorder=3 }
if(!is.numeric(dots$Thin)){
  Thin=4
} else{
  Thin=round(abs(dots$Thin))
}

if(!is.numeric(dots$m)){
  m=2
} else{
  m=round(abs(dots$m))
}

if(!is.numeric(dots$n.chains)){
  n.chains=2
} else{
  n.chains=round(abs(dots$n.chains))
}

if(!is.numeric(dots$freqs)){
  freqs='internal'
} else{
  freqs=dots$freqs
}

if(!is.numeric(dots$polyorder)){
  polyorder=3
} else{
  if (dots$polyorder!=abs(round(dots$polyorder))){stop("polyorder must be an integer >0.")}
  if (dots$polyorder < 0){stop("polyorder must be an integer >=0.")}
  polyorder=round(abs(dots$polyorder))
}

y=y[order(ti.mu,decreasing = FALSE)]; ti.sd=ti.sd[order(ti.mu,decreasing = FALSE)]
ti.mu=ti.mu[order(ti.mu,decreasing = FALSE)]
# Freq vector
if (is.character(freqs)){
  freqs=seq(1/(max(ti.mu)-min(ti.mu)),floor(0.5 * length(y)) /(max(ti.mu)-min(ti.mu)),
            by = 1/(max(ti.mu)-min(ti.mu)))
} else if (is.numeric(freqs)){
   if(sum(freqs<0)>0){stop("review the frequnecy vector")}

}
if(length(freqs)>50 &length(freqs)>1){
  cat("Frequency vector has more that 50 frequencies","\n",
      "Enter 'ok' to proceed or provide the maximum number of frequencies")
  MaxFreq=readline(prompt="Max frequency: ")
  if(MaxFreq!="ok"){
    MaxF=round(abs(as.numeric(MaxFreq)))
    seq(from=min(freqs), to=max(freqs), length.out=MaxF)
  }
}

# JAGS model
STR1="model {
for(i in 1:n) {
  y[i]~dnorm(mu[i], precision)
  #y_PRED[i]~dnorm(mu[i], precision)
  mu[i]<-const"
STR2="+"
if (polyorder>0){
  for(i in 1:polyorder){
    K=paste("alpha",i,'*ti.sim.centr[i]^',i,sep='')
    STR2=paste(STR2,K,'+')
  }
}
STR3="inprod(X[i,],IBeta)
  }
for(j in 1:(2*n.freqs)) {
  beta[j] ~ dnorm(0,0.01)
}
for(j in 1:n.freqs) {
  M[j]~dunif(0,m)
  Ind[j]~dbern(M[j]/n.freqs)
  IBeta[j]<-Ind[j]*beta[j]
  Ind[j+n.freqs] <-Ind[j]
  IBeta[j+n.freqs]<-Ind[j+n.freqs]*beta[j+n.freqs]
}
for(i in 1:n) {
  for(j in 1:n.freqs) {
    X[i,j] <- sin(2*pi*ti.sim[i]*freqs[j])
    }
  for(j in (n.freqs+1):(2*n.freqs)) {
    X[i,j] <- cos(2*pi*ti.sim[i]*freqs[j-n.freqs])
    }
}
for (i in 1:n) {
  ti.sim.tmp[i]~dnorm(ti.mu[i],1/ti.sd[i]^2)
  }
ti.sim<-sort(ti.sim.tmp)
ti.sim.centr<-(ti.sim-ti.mu.mu)/ti.sd.mu
Spectrum<-(IBeta[1:n.freqs]^2+IBeta[(n.freqs+1):(2*n.freqs)]^2)/2
const~dnorm(const.mean,const.prec)
precision~dgamma(0.01,0.01)"

STR4=""
if (polyorder>0){
  for(i in 1:polyorder){
    M=paste("alpha",i,"~dnorm(0,",'0.01)',sep='')
    STR4=paste(STR4,'\n',M)
  }
}
STR5="}"

ti.mu.mu=mean(ti.mu)
ti.sd.mu=mean(ti.sd)

modelstring=paste(STR1,STR2,STR3,STR4,STR5)
data=list(y=y, ti.mu=ti.mu,ti.sd=ti.sd, n=length(ti.mu), n.freqs=length(freqs), freqs=freqs,
          pi=pi,const.mean=0, const.prec=0.01,m=m,ti.mu.mu=ti.mu.mu,ti.sd.mu=ti.sd.mu) #
inits = parallel.seeds("base::BaseRNG", n.chains)
for(k in (1:n.chains)){
  inits[[k]]$const = 0
  inits[[k]]$beta = rep(0,2*length(freqs))
  inits[[k]]$ti.sim.tmp=ti.mu
  if (polyorder>0){
    for(i in 1:polyorder){
      eval(parse(text = paste("inits[[k]]$alpha",i,'=0',sep='')))
    }
  }
}

model=jags.model(textConnection(modelstring), data=data, inits=inits,n.chains=n.chains)
update(model,n.iter=n.sim,thin=Thin)

variable.names="const"
if (polyorder>0){
  for (i in 1:polyorder){
    variable.names=c(variable.names,paste("alpha",i,sep=''))
  }
}
output=coda.samples(model=model,variable.names=c("Spectrum","Ind", variable.names,"ti.sim","beta","precision"),
                    n.iter=n.sim, thin=Thin,n.chains=n.chains)
DIC = dic.samples(model=model,n.iter=n.sim,thin=Thin)

Sim.Objects=JAGS.objects(output)
Sim.Objects$freqs=freqs
Sim.Objects$JAGS=output
Sim.Objects$DIC=DIC
# Cross Validation
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
  Cores=min(c(parallel::detectCores()-1,folds));  cl = parallel::makeCluster(Cores);
  doParallel::registerDoParallel(cl)

  BSFCV=function(i,fold,y,ti.mu,ti.sd,freqs,modelstring,n.sim,polyorder,m){
    Y=y; Y[fold==i]=NA;
    data=list(y=Y, ti.mu=TI.SIM,ti.sd=ti.sd,n=length(y), n.freqs=length(freqs),
              freqs=freqs, pi=pi,const.mean=0, const.prec=0.01, polyorder=polyorder,m=m,ti.mu.mu=ti.mu.mu,ti.sd.mu=ti.sd.mu) #
    model=jags.model(textConnection(modelstring), data=data,n.chains=1)
    update(model,n.iter=n.sim,thin=Thin)
    output=coda.samples(model=model,variable.names=c("y"), n.iter=n.sim, thin=Thin)
    return(output)
  }


  CVRES=foreach(i=1:folds,.export=c('jags.model','coda.samples')) %dopar%
    BSFCV(i,fold=fold,y=y,ti.mu=ti.mu,ti.sd=ti.sd, freqs=freqs,
          modelstring=modelstring, n.sim=n.sim, polyorder=polyorder,m=m)

  parallel::stopCluster(cl)

  pred_y = array(NA, dim=c(dim(JAGS.objects(CVRES[[1]])$y[,fold==(1)])[1],length(ti.mu)))
  colnames(pred_y)= names(y)

  for (i in 1:folds){
    pred_y[,fold==(i)] = JAGS.objects(CVRES[[i]])$y[,fold==(i)]
  }

  Sim.Objects$CVpred=pred_y

}
Sim.Objects$y=y
Sim.Objects$polyorder=polyorder
Sim.Objects$ti.mu=ti.mu
class(Sim.Objects)='tuts_BFS'
return(Sim.Objects)
}

#' Summary of tuts_BFS objects.
#'
#' \code{summary.tuts_BFS} prints summaries of tuts_BFS objects.
#'
#' @param object A tuts_BFS object.
#' @param ... A list of optional parameters: \cr
#'  - burn: burn-in parameter ranging from 0 to 0.7, the default value is 0.\cr
#'  - CI: confidence interval, the default value is set to 0.99.
#' @examples
#' \donttest{
#' #1. Import or simulate the data (a simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' #2. Fit the model:
#' n.sim=1000
#' BFS=tubfs(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=2, CV=TRUE)
#'
#' #3. Generate summary results (optional parameters are listed in brackets):
#' summary(BFS)                               # Summary results (CI, burn).
#' summary(BFS,burn=0.2)                      # Results after 20% of burn-in (CI).
#'}
#'
#' @export
summary.tuts_BFS = function(object, ...) {
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
  n.sim=dim(object$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}
  #
  cat('\n')
  cat('Bayesian Frequency Selection:\n')
  cat('-----------------------------\n')

  IND=object$Ind[BURN:dim(object$Ind)[1],]
  IND=IND[,1:(dim(IND)[2]/2)]+IND[,(dim(IND)[2]/2+1):dim(IND)[2]]
  IND[IND==2]=1

  PWR=object$Spectrum[BURN:dim(object$Spectrum)[1],]
  PWR.lwr=apply(PWR,2,'quantile',0.005)
  PWR.med=apply(PWR,2,'quantile',0.5)
  PWR.upr=apply(PWR,2,'quantile',0.995)

  Frequency=object$freqs
  Period=1/Frequency
  Probability=apply(IND,2,sum)/(dim(IND)[1])

  TABLE=data.frame(Frequency,Period,PWR.med,PWR.upr,PWR.lwr,Probability)
  rownames(TABLE)=paste("Pwr",1:length(Frequency),sep=" ")
  print(TABLE)
  #
  cat('\n')
  cat('Regression Parameters and estimates of timing:\n')
  cat('----------------------------------------------\n')

  const=object$const[BURN:length(object$const)]
  const.lwr=quantile(const,(1-CI)/2)
  const.med=quantile(const,0.5)
  const.upr=quantile(const,1-(1-CI)/2)
  constName="const"

  lwr=med=upr=NAMES=NA
  if(object$polyorder>0){
    alpha.lwr=NA
    alpha.med=NA
    alpha.upr=NA
    alphaNames=NA
    for (i in 1:object$polyorder){
      eval(parse(text=paste("alpha",i,"=object$alpha",i,"[BURN:dim(object$alpha",i,")[1],]",sep="")))
      alpha.lwr[i]=eval(parse(text=paste("quantile(alpha",i,",(1-CI)/2)",sep="")))
      alpha.med[i]=eval(parse(text=paste("quantile(alpha",i,",0.5)",sep="")))
      alpha.upr[i]=eval(parse(text=paste("quantile(alpha",i,",1-(1-CI)/2)",sep="")))
      alphaNames[i]=c(paste("alpha",i,sep=""))
    }
  }

  Spectrum=object$Spectrum[BURN:dim(object$Spectrum)[1],]
  Spectrum.lwr=apply(Spectrum,2,'quantile',(1-CI)/2)
  Spectrum.med=apply(Spectrum,2,'quantile',0.5)
  Spectrum.upr=apply(Spectrum,2,'quantile',1-(1-CI)/2)
  SpectrumNames=names(Spectrum.med)


  beta=object$beta[BURN:dim(object$beta)[1],]
  beta.lwr=apply(beta,2,'quantile',(1-CI)/2)
  beta.med=apply(beta,2,'quantile',0.5)
  beta.upr=apply(beta,2,'quantile',1-(1-CI)/2)
  betaNames=names(beta.med)

  precision=object$precision[BURN:length(object$precision)]
  precision.lwr=quantile(precision,(1-CI)/2)
  precision.med=quantile(precision,0.5)
  precision.upr=quantile(precision,1-(1-CI)/2)
  precisionName="precision"

  ti=object$ti.sim[BURN:dim(object$ti.sim)[1],]
  ti.lwr=apply(ti,2,'quantile',(1-CI)/2)
  ti.med=apply(ti,2,'quantile',0.5)
  ti.upr=apply(ti,2,'quantile',1-(1-CI)/2)
  tiNames=names(ti.med)

  if(!exists("alphaNames")){
    LWR=c(const.lwr,Spectrum.lwr, beta.lwr,precision.lwr,ti.lwr)
    MED=c(const.med,Spectrum.med, beta.med,precision.med,ti.med)
    UPR=c(const.upr,Spectrum.upr, beta.upr,precision.upr,ti.upr)
    TABLE2=data.frame(LWR,MED,UPR)
    row.names(TABLE2)=c(constName,SpectrumNames,betaNames,precisionName,tiNames)

  }else{
    LWR=c(const.lwr,alpha.lwr,Spectrum.lwr, beta.lwr,precision.lwr,ti.lwr)
    MED=c(const.med,alpha.med,Spectrum.med, beta.med,precision.med,ti.med)
    UPR=c(const.upr,alpha.upr,Spectrum.upr, beta.upr,precision.upr,ti.upr)
    TABLE2=data.frame(LWR,MED,UPR)
    row.names(TABLE2)=c(constName,alphaNames,SpectrumNames,betaNames,precisionName,tiNames)
  }


  colnames(TABLE2)=c(paste(round((1-CI)/2,3)*100,"%",sep=""),'50%',paste(round(1-(1-CI)/2,3)*100,"%",sep=""))
  print(TABLE2)

  #
  cat('\n')
  cat('Deviance information criterion:\n')
  cat('-------------------------------\n')
  print(object$DIC)
  cat('-------------------------------\n')
}
####################################################################################################
#' Graphical summaries and diagnostics of tuts_BFS objects.
#'
#' \code{plot.tuts_BFS} plots and diagnostics of tuts_BFS objects.
#'
#' @param x A tuts_BFS objects.
#' @param type plot type with the following options:\cr
#'  - 'periodogram' plots estimates of power spectrum.  \cr
#'  - 'predTUTS' plots one step predictions of the model. \cr
#'  - 'GR' plots Gelman-Rubin diagnostics. \cr
#'  - 'cv' plots 5-fold cross validation. \cr
#'  - 'mcmc' plots diagnostics of MCMC/JAGS objects. \cr
#'  - 'volatility' plots volatility realizations. \cr
#' @param ... list of optional parameters:\cr
#'  - burn: burn-in parameter ranging from 0 to 0.7 with default value set to 0. \cr
#'  - CI: credible interval ranging from 0.3 to 1 with default value set to 0.95.
#' @examples
#' \donttest{
#' #1. Import or simulate the data (a simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' #2. Fit the model:
#' n.sim=1000
#' BFS=tubfs(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=2, CV=TRUE)
#'
#' #3. Generate plots and diagnostics of the model (optional parameters are listed in brackets):
#' plot(BFS,type='periodogram')               # spectral analysis (CI, burn).
#' plot(BFS,type='predTUTS', CI=0.99)         # One step predictions (CI, burn).
#' plot(BFS,type='cv')                        # 5 fold cross validation plot (CI, burn).
#' plot(BFS,type='GR')                        # Gelman-Rubin diagnostics (CI, burn).
#' plot(BFS,type='mcmc')                      # mcmc diagnostics.
#' plot(BFS,type='volatility')                # Volatility realizaitons.
#'}
#' @export
plot.tuts_BFS = function(x, type, ...) {
  dots = list(...)
  #type=dots$type
  if (sum(type==c('periodogram','predTUTS','GR','cv','mcmc','volatility'))==0){
    stop('type should be set as either periodogram, predTUTS, GR, cv, mcmc or volatility')
  }

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
  n.sim=dim(x$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}

  #
  graphics::par(mfrow=c(1,1))

  if(type=='periodogram') {
    IND=x$Ind[BURN:dim(x$Ind)[1],]
    IND=IND[,1:(dim(IND)[2]/2)]+IND[,(dim(IND)[2]/2+1):dim(IND)[2]]
    IND[IND==2]=1

    PWR=x$Spectrum[BURN:dim(x$Spectrum)[1],]
    PWR.lwr=apply(PWR,2,'quantile',(1-CI)/2)
    PWR.med=apply(PWR,2,'quantile',0.5)
    PWR.upr=apply(PWR,2,'quantile',(1+CI)/2)

    Frequency=x$freqs
    Period=1/Frequency
    Probability=apply(IND,2,sum)/(dim(IND)[1])
    DIV=round(max(PWR.upr)/4,- floor(log10(max(PWR.upr)/4)))
    YLAB=seq(from=0,to=4*DIV, by=DIV)
    #
    graphics::par(mar=c(5,5,5,5))
    graphics::par(mfrow=c(1,1))
    graphics::plot(x=Frequency,y=PWR.med,pch=20,xlab="frequency", ylab="",ylim=c(-0.4*max(PWR.upr), 1.25*max(PWR.upr)),yaxt="n", main="BFS Spectrum")
    graphics::axis(side=2, at=YLAB)
    graphics::mtext("Power", side=2, line=2.5, at=max(PWR.upr)/2)

    graphics::lines(x=Frequency,y=PWR.med,type='l',col="gray")
    options(warn=-1)
    graphics::arrows(x0=Frequency,y0=PWR.lwr,x1=Frequency,y1=PWR.upr,code=3,length=0.04,angle=90,col='black')
    options(warn=0)
    graphics::legend("topright", paste("Power with ", CI*100, "%CI",sep=""),lty=c(NA),pch=c(20),lwd=c(1), col="black",border="white")
    graphics::par(new=TRUE)
    graphics::plot(y=Probability*max(PWR.upr)*0.04,x=Frequency,pch=20,xlab="",ylab="",yaxt="n", ylim=c(0,max(PWR.upr)*0.2))
    graphics::lines(y=Probability*max(PWR.upr)*0.04,x=Frequency,type='h',col='gray',xlab="",ylab="",yaxt="n", ylim=c(0,max(PWR.upr)*0.2))
    graphics::axis(side=4, at=c(0,max(PWR.upr)*0.04),labels=c("0%","100%"))
    graphics::mtext("Probability", side=4, line=2.5,at=max(PWR.upr)*0.04/2)
  }

  #
  if(type=='cv') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    MAIN="Cross-Validation: One step out of sample predictions"

    graphics::plot(x=x$y,y=PRED,xlab="Actual",ylab="Predicted", main=MAIN, pch=18)
    graphics::abline(0,1,col='blue')

    RSQ=cor(x$y,PRED)^2*100

    LAB = bquote(italic(R)^2 == .(paste(format(RSQ, digits = 1),"%",sep="")))
    graphics::text(x=(min(x$y)+0.9*(max(x$y)-min(x$y))),y=(min(PRED)+0.1*(max(PRED)-min(PRED))),LAB)
  }

  #
  if(type=='predTUTS') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED.LWR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',(1-CI)/2)
    PRED.MED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    PRED.UPR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',(1+CI)/2)

    ti.sim=apply(x$ti.sim[BURN:dim(x$ti.sim)[1],],2,'quantile',0.5)

    MAIN=paste("One step out of sample predictions at CI= ", CI*100,"%",sep='')
    graphics::plot(y=x$y,x=x$ti.mu,type='l',main=MAIN,ylab="Observations",xlab='time',lwd=2,
         ylim=c(min(x$CVpred),1.2*max(x$CVpred)), xlim=c(min(x$ti.mu,ti.sim),max(x$ti.mu,ti.sim)))

    graphics::lines(y=PRED.LWR,x=ti.sim,type='l',col='blue',lwd=1,lty=2)
    graphics::lines(y=PRED.MED,x=ti.sim,type='l',col='blue',lwd=1,lty=1)
    graphics::lines(y=PRED.UPR,x=ti.sim,type='l',col='blue',lwd=1,lty=2)

    graphics::legend("topright",legend = c("Observed","Upper CI","Medium","Lower CI"),
           col=c("black","blue","blue","blue"),lwd=c(2,1,1,1),lty=c(1,2,1,2))
  }
  #
  if(type=='GR') {
    if(burn>0){ABURNIN=TRUE} else{ABURNIN=FALSE}

    ALL_Objects=names(x)
    Remove=c("Ind", "Spectrum", "freqs" , "beta","const", "JAGS" , "DIC",
             "CVpred","y","polyorder" , "ti.mu" ,"ti.sim")

    for (i in 1:length(Remove)){ALL_Objects=ALL_Objects[!ALL_Objects==Remove[i]]}
    PN_Objects=c("const",ALL_Objects)

    GELMAN.PN=array(NA,dim=c(length(PN_Objects),3))
    rownames(GELMAN.PN)=PN_Objects

    for (i in 1:length(PN_Objects)){
      PRM=grep(c(PN_Objects[i]),colnames(x$JAGS[[1]]))
      GELMAN.PN[i,1:2]=gelman.diag(x$JAGS[,PRM],multivariate=FALSE,confidence =CI,autoburnin=ABURNIN)$psrf
      GELMAN.PN[i,3]=i
    }

    graphics::par(mfrow=c(3,1),oma = c(0, 0, 2, 0))
    graphics::plot(y=GELMAN.PN[1,1], x=GELMAN.PN[1,3],ylim=c(0,max(GELMAN.PN[,1:2])),xlim=c(1,(dim(GELMAN.PN)[1]))
         ,xaxt='n',ylab="Factor",  main="Parameters of polynomial regression and precision of the model",xlab="")
    for(i in 2:dim(GELMAN.PN)[1]){
      graphics::points(x=GELMAN.PN[i,3],y=GELMAN.PN[i,1])
    }

    for(i in 1:dim(GELMAN.PN)[1]){
      graphics::arrows(x0=GELMAN.PN[i,3],
             y0=GELMAN.PN[i,1],
             x1=GELMAN.PN[i,3],
             y1=GELMAN.PN[i,2],
             code=3,length=0.04,angle=90,col='darkgray')
    }
    graphics::text(x=seq(1,length(PN_Objects),by=1),y=-max(GELMAN.PN[,1:2])/7, srt = 00, adj= 0.5, xpd = TRUE, labels =PN_Objects, cex=0.65)
    graphics::legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")
    graphics::legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")
    graphics::abline(h=1);

    #
    par_to_use = grep('Spectrum',colnames(x$JAGS[[1]]))
    GELMAN.SPC <- matrix(NA, nrow=length(par_to_use), ncol=3)
    for (v in 1:length(par_to_use)) {
      GELMAN.SPC[v,1:2] <- gelman.diag(x$JAGS[,par_to_use[v]],multivariate=FALSE,confidence =CI,autoburnin=ABURNIN)$psrf
      GELMAN.SPC[v,3] <- v
    }
    GELMAN.SPC[is.na(GELMAN.SPC)]=1
    graphics::plot(GELMAN.SPC[,1], main="Power estimates", xaxt='n',ylab="Factor",xlab="",
         ylim=c(0,max(GELMAN.SPC[,1:2])),xlim=c(1,(dim(GELMAN.SPC)[1])))

    for(i in 1:dim(GELMAN.SPC)[1]){
      graphics::arrows(x0=GELMAN.SPC[i,3],
             y0=GELMAN.SPC[i,1],
             x1=GELMAN.SPC[i,3],
             y1=GELMAN.SPC[i,2],
             code=3,length=0.04,angle=90,col='darkgray')
    }
    graphics::legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")
    graphics::abline(h=1)

    if (dim(GELMAN.SPC)[1]<11){
      graphics::text(x=seq(1,dim(GELMAN.SPC)[1],by=1),y=-max(GELMAN.SPC[,1:2]) /2.5, srt = 00, adj= 0, xpd = TRUE, srt = 90,
           labels =paste("P",1:dim(GELMAN.SPC)[1]), cex=0.65)}else {
      graphics::text(x=seq(5,dim(GELMAN.SPC)[1],by=5),y=-max(GELMAN.SPC[,1:2]) /2.5, srt = 00, adj= 0, xpd = TRUE, srt = 90,
           labels =paste("P",seq(5,dim(GELMAN.SPC)[1],by=5)), cex=0.65)}

    #
    par_to_use = grep('ti.sim',colnames(x$JAGS[[1]]))
    GELMAN.TIS <- matrix(NA, nrow=length(par_to_use), ncol=3)
    for (v in 1:length(par_to_use)) {
      GELMAN.TIS[v,1:2] <- gelman.diag(x$JAGS[,par_to_use[v]],multivariate=FALSE,confidence =CI,autoburnin=ABURNIN)$psrf
      GELMAN.TIS[v,3] <- v
    }
    GELMAN.TIS[is.na(GELMAN.TIS)]=1
    graphics::plot(GELMAN.TIS[,1], main="Estimated timings of observations", xaxt='n',ylab="Factor",xlab="",
         ylim=c(0,max(GELMAN.TIS[,1:2])),xlim=c(1,(dim(GELMAN.TIS)[1])))

    for(i in 1:dim(GELMAN.TIS)[1]){
      graphics::arrows(x0=GELMAN.TIS[i,3],
             y0=GELMAN.TIS[i,1],
             x1=GELMAN.TIS[i,3],
             y1=GELMAN.TIS[i,2],
             code=3,length=0.04,angle=90,col='darkgray')
    }
    graphics::abline(h=1)

    if (dim(GELMAN.TIS)[1]<21){
      graphics::text(x=seq(1,dim(GELMAN.TIS)[1],by=1),y=-max(GELMAN.TIS[,1:2])/2.5, srt = 00, adj= 0.5, xpd = TRUE,srt = 90,
           labels =paste("t",1:dim(GELMAN.TIS)[1]), cex=0.65)}else {
      graphics::text(x=seq(5,dim(GELMAN.TIS)[1],by=5),y=-max(GELMAN.TIS[,1:2])/2.5, srt = 00, adj= 0.5, xpd = TRUE,srt = 90,
           labels =paste("t",seq(5,dim(GELMAN.TIS)[1],by=5)), cex=0.65)}

    graphics::legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")


    TITLE=paste("Gelman-Rubin diagnostics: Potential scale reduction factors with the upper confidence bounds at ",
               CI*100,"%",sep="")

    graphics::mtext(TITLE, outer = TRUE, cex = 0.9)
    #
    graphics::par(mfrow=c(1,1))
  }
  if(type=='volatility') {
    graphics::plot(sqrt(sqrt(1/sqrt(x$precision[BURN:dim(x$precision)[1],]))),type='l', xlab='Sim ID',ylab='Std Deviation',main='Standard Deviaiton')
  }
  if(type=='mcmc') {

    ALL_Objects=names(x)
    Remove=c("Ind","const","Spectrum", "beta","freqs","JAGS","DIC","y","polyorder","ti.mu","CVpred",
             "precision", "ti.sim" )

    for (i in 1:length(Remove)){ALL_Objects=ALL_Objects[!ALL_Objects==Remove[i]]}
    PN_Objects=c("const",ALL_Objects,"Spectrum","ti.sim" ,"beta","precision"  )

    mcmcplot(x$JAGS, parms=c(PN_Objects))
  }

}
