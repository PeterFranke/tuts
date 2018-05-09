#' Time-uncertain polynomial regression of counting series
#'
#' \code{tupoispn}  poisson polynomial regression of the N-th order of time-uncertain time series.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates of timing of observations.
#' @param ti.sd A vector of standard deviations of timing.
#' @param polyorder Order of the polynomial regression.
#' @param n.sim A number of simulations.
#' @param CV cross-validation indicator.
#' @param ... A list of optional parameters. The list contains polynomial order with the default value
#' set to polyorder=3, thinning parameter, with the default value set to Thin=4,
#' and the number of mcmc chains with the default value set to n.chains=2.
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0), trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' y=round(y-min(y))
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Set parameters and run the polynomial regression of counting series:
#' polyorder=2
#' n.sim=1000
#' PPN=tupoispn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder,n.sim=n.sim, CV=TRUE)
#'
#' # Generate summary results (optional parameters are listed in brackets):
#' summary(PPN)                           # Summary statistics (burn, CI).
#'
#' # Plots and diagnostics (optional parameters are listed in brackets):
#' plot(PPN,type='predTUTS',CI=0.95)      # One step out of salmple predictions of the model (CI, burn).
#' plot(PPN,type='cv',burn=0.3)           # 5 fold cross-validation (CI, burn).
#' plot(PPN,type='GR',CI=0.95)            # Gelman-Rubin diagnostic (CI).
#' plot(PPN,type='mcmc')                  # MCMC diagnostics.
#' plot(PPN,type='lambda')                # Volatility realizaitons.
#'
#' @export
#'
tupoispn=function(y,ti.mu,ti.sd, n.sim, CV=FALSE, ... ){
# Data checking and basic operations ---------------------------------------------------------------------
if (length(y)*4!=length(ti.mu)*2+length(ti.sd)*2){stop("Vectors y, ti.mu and ti.sd should be of equal lengths.")}
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
  if(missing(...)){Thin=4; n.chains=2; polyorder=3 }
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

  if(!is.numeric(dots$polyorder)){
    polyorder=3
  } else{
    if (dots$polyorder < 0){stop("polyorder must be an integer >=0.")}
    if (dots$polyorder!=abs(round(dots$polyorder))){stop("polyorder must be an integer >0.")}
    polyorder=round(abs(dots$polyorder))
  }

y=y[order(ti.mu,decreasing = FALSE)]; ti.sd=ti.sd[order(ti.mu,decreasing = FALSE)]
ti.mu=ti.mu[order(ti.mu,decreasing = FALSE)]

# JAGS model --------------------------------------------------------------------------------------------
modelstring= "model {
for(i in 1:n) {
  y[i]~dpois(lambda[i])
  log(lambda[i])<-const+inprod(ti.sim.MX[i,],alpha)
}
for (i in 1:n) {
  ti.sim.tmp[i]~dnorm(ti.mu[i],1/ti.sd[i])
}
ti.sim<-sort(ti.sim.tmp)

for (i in 1:n) {
  for(j in 1:(polyorder)){
    ti.sim.MX[i,j]<- ti.sim[i]^j
  }
}
for(i in 1:polyorder){
  alpha[i]~dnorm(0,0.01)
}
const~dnorm(0,0.01)
}"

# JAGS data ---------------------------------------------------------------------------------------------
data=list(y=y, ti.mu=ti.mu,ti.sd=ti.sd, n=length(ti.mu),polyorder=polyorder)

init=parallel.seeds("base::BaseRNG", n.chains)
for(k in (1:n.chains)){
  init[[k]]$ti.sim.tmp=ti.mu
}

model=jags.model(textConnection(modelstring), data=data, init=init,n.chains=n.chains)
update(model,n.iter=n.sim,thin=Thin)

output=coda.samples(model=model, variable.names=c("lambda","const","alpha","ti.sim")
                      ,n.iter=n.sim, thin=Thin,n.chains=n.chains)
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
    Y=y; Y[fold==i]=NA;
    data=list(y=Y, ti.mu=TI.SIM,ti.sd=ti.sd, n=length(y),polyorder=polyorder)
    model=jags.model(textConnection(modelstring), data=data,n.chains=1)
    update(model,n.iter=n.sim,thin=Thin)
    output=coda.samples(model=model,variable.names=c("y"), n.iter=n.sim, thin=Thin)
    return(output)
  }

  CVRES=foreach(i=1:folds,.export=c('jags.model','coda.samples')) %dopar%
    BSFCV(i,fold=fold,y=y,ti.mu=TI.SIM,ti.sd=ti.sd, modelstring=modelstring, n.sim=n.sim)

  stopCluster(cl)

  pred_y = array(NA, dim=c(dim(JAGS.objects(CVRES[[1]])$y[,fold==(1)])[1],length(ti.mu)))
  colnames(pred_y)= paste("y[",1:length(y),"]",sep="")

  for (i in 1:folds){
    pred_y[,fold==(i)] = JAGS.objects(CVRES[[i]])$y[,fold==(i)]
  }
  Sim.Objects$CVpred=pred_y
  }
Sim.Objects$y=y
Sim.Objects$ti.mu=ti.mu
Sim.Objects$polyorder=polyorder
class(Sim.Objects)='tuts_poisPN'
return(Sim.Objects)
}

#' Prints summaries of tuts_poisPN objects.
#'
#' \code{summary.tuts_poisPN} Prints summaries of tuts_poisPN objects.
#'
#' @param x A tuts_poisPN object.
#' @param ... list of optional parameters. The list contains burn-in parameter
#' ranging from 0 to 0.5, with the default value burn=0, and the credible interval parameter
#' ranging between 0.5 and 1, with the default CI=0.99.
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0), trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' y=round(y-min(y))
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Set parameters and run the polynomial regression of counting series:
#' polyorder=2
#' n.sim=1000
#' PPN=tupoispn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder,n.sim=n.sim, CV=TRUE)
#'
#' # Generate summary results (optional parameters are listed in brackets):
#' summary(PPN)                           # Summary statistics (burn, CI).
#'
#' @export
summary.tuts_poisPN = function(x, ...) {
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
  n.sim=dim(x$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}

  # ----------------------------------------------------------------------------
  cat('\n')
  cat('Regression Parameters and estimates of timing:\n')
  cat('----------------------------------------------\n')

  const=x$const[BURN:length(x$const)]
  const.lwr=quantile(const,(1-CI)/2)
  const.med=quantile(const,0.5)
  const.upr=quantile(const,1-(1-CI)/2)
  constName="const"

  lwr=med=upr=NAMES=NA
  polyorder=x$polyorder
  alpha=x$alpha[BURN:dim(x$alpha)[1],]
  alpha.lwr=apply(alpha,2,'quantile',(1-CI)/2)
  alpha.med=apply(alpha,2,'quantile',0.5)
  alpha.upr=apply(alpha,2,'quantile', 1-(1-CI)/2)
  alphaNames=colnames(alpha)

  lambda=x$lambda[BURN:dim(x$lambda)[1],]
  lambda.lwr=apply(lambda,2,'quantile',(1-CI)/2)
  lambda.med=apply(lambda,2,'quantile',0.5)
  lambda.upr=apply(lambda,2,'quantile',1-(1-CI)/2)
  lambdaNames=names(lambda.med)

  ti=x$ti.sim[BURN:dim(x$ti.sim)[1],]
  ti.lwr=apply(ti,2,'quantile',(1-CI)/2)
  ti.med=apply(ti,2,'quantile',0.5)
  ti.upr=apply(ti,2,'quantile',1-(1-CI)/2)
  tiNames=names(ti.med)

  LWR=c(const.lwr,alpha.lwr,lambda.lwr,ti.lwr)
  MED=c(const.med,alpha.med,lambda.med,ti.med)
  UPR=c(const.upr,alpha.upr,lambda.upr,ti.upr)
  TABLE2=data.frame(LWR,MED,UPR)
  row.names(TABLE2)=c(constName,alphaNames,lambdaNames,tiNames)

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
#' Graphical summaries and diagnostics of tuts_poisPN objects
#'
#' \code{plot.tuts_poisPN} plots summaries and diagnostics of a tuts_poisPN object.
#'
#' @param x A tuts_poisPN object.
#' @param type plot/diagnostic type (options:'predTUTS' plots one step out of sample predictions of the model,
#' 'GR' plots Gelman-Rubin diagnostics, 'cv' plots 5-fold cross validation, 'mcmc' plots diagnostics
#'  of mcmc objects, and 'lambda' plots lambda realizations).
#' @param ... list of optional parameters: 'burn' (burn-in parameter ranging from 0 to 0.7 with
#'  default value set to 0), and CI (credible interval ranging from 0.3 to 1 with
#'  default value set to 0.95).
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0), trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' y=round(y-min(y))
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' # Set parameters and run the polynomial regression of counting series:
#' polyorder=2
#' n.sim=1000
#' PPN=tupoispn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder,n.sim=n.sim, CV=TRUE)
#'
#' # Plots and diagnostics (optional parameters are listed in brackets):
#' plot(PPN,type='predTUTS',CI=0.95)      # One step out of salmple predictions of the model (CI, burn).
#' plot(PPN,type='cv',burn=0.3)           # 5 fold cross-validation (CI, burn).
#' plot(PPN,type='GR',CI=0.95)            # Gelman-Rubin diagnostic (CI).
#' plot(PPN,type='mcmc')                  # MCMC diagnostics.
#' plot(PPN,type='lambda')                # Volatility realizaitons.
#'
#' @export
plot.tuts_poisPN = function(x, type, ...) {
  if (sum(type==c('predTUTS','GR','cv','mcmc','lambda'))==0){
    stop('type should be set as either predTUTS, GR, cv, mcmc or lambda')
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
  n.sim=dim(x$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}

  par(mfrow=c(1,1))
  ##########################################################################
  if(type=='cv') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    MAIN="Cross-Validation: One step out of sample predictions"

    plot(x=x$y,y=PRED,xlab="Actual",ylab="Predicted", main=MAIN, pch=18)
    abline(0,1,col='blue')

    RSQ=cor(x$y,PRED)^2* 100

    LAB = bquote(italic(R)^2 == .(paste(format(RSQ, digits = 0),"%",sep="")))
    text(x=(min(x$y)+0.9*(max(x$y)-min(x$y))),y=(min(PRED)+0.1*(max(PRED)-min(PRED))),LAB)
  }
  ##########################################################################
  if(type=='predTUTS') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED.LWR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',(1-CI)/2)
    PRED.MED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    PRED.UPR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',1-(1-CI)/2)

    ti.sim=c(x$ti.mu[1],apply(x$ti.sim[BURN:dim(x$ti.sim)[1],],2,'quantile',0.5))

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
  if(type=='GR') {
    if(burn>0){ABURNIN=TRUE} else{ABURNIN=FALSE}

    ALL_Objects=names(x)
    Remove=c("Ind", "Spectrum", "freqs" , "beta","const", "JAGS" , "DIC",
             "CVpred","y","polyorder" , "ti.mu" ,"ti.sim","IBeta","lambda")

    for (i in 1:length(Remove)){ALL_Objects=ALL_Objects[!ALL_Objects==Remove[i]]}
    PN_Objects=c("const",paste("alpha[",1:x$polyorder,"]",sep=""))

    GELMAN.PN=array(NA,dim=c((1+x$polyorder),3))
    rownames(GELMAN.PN)=PN_Objects

    NCOL=seq(1:length(colnames(x$JAGS[[1]])))
    for (i in 1:dim(GELMAN.PN)[1]){
      PRM=NCOL[colnames(x$JAGS[[1]])==PN_Objects[i]]
      GELMAN.PN[i,1:2]=gelman.diag(x$JAGS[,PRM],multivariate=FALSE,confidence =CI,autoburnin=ABURNIN)$psrf
      GELMAN.PN[i,3]=i
    }

    par(mfrow=c(3,1))
    plot(y=GELMAN.PN[1,1], x=GELMAN.PN[1,3],ylim=c(0,max(GELMAN.PN[,1:2])),xlim=c(1,(dim(GELMAN.PN)[1]))
         ,xaxt='n',ylab="Factor",  xlab="Parameters of polynomial regression",
         main=paste("Gelman-Rubin diagnostics: Potential scale reduction factors \n with the upper confidence bounds at ",
                    CI*100,"%",sep=""))
    for(i in 2:dim(GELMAN.PN)[1]){
      points(x=GELMAN.PN[i,3],y=GELMAN.PN[i,1])
    }

    for(i in 1:dim(GELMAN.PN)[1]){
      arrows(x0=GELMAN.PN[i,3],
             y0=GELMAN.PN[i,1],
             x1=GELMAN.PN[i,3],
             y1=GELMAN.PN[i,2],
             code=3,length=0.04,angle=90,col='darkgray')
    }
    text(x=seq(1,length(PN_Objects),by=1),y=-max(GELMAN.PN[,1:2])/7, srt = 00, adj= 0.5, xpd = TRUE, labels =PN_Objects, cex=0.65)
    legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")
    abline(h=1);

    ##############
    ###############
    par_to_use = grep('ti.sim',colnames(x$JAGS[[1]]))
    GELMAN.TIS <- matrix(NA, nrow=length(par_to_use), ncol=3)
    for (v in 1:length(par_to_use)) {
      GELMAN.TIS[v,1:2] <- gelman.diag(x$JAGS[,par_to_use[v]],multivariate=FALSE,confidence =CI,autoburnin=ABURNIN)$psrf
      GELMAN.TIS[v,3] <- v
    }
    GELMAN.TIS[is.na(GELMAN.TIS)]=1
    plot(GELMAN.TIS[,1], xlab="Estimated timings of observations", xaxt='n',ylab="Factor",
         ylim=c(0,max(GELMAN.TIS[,1:2])),xlim=c(1,(dim(GELMAN.TIS)[1])))

    for(i in 1:dim(GELMAN.TIS)[1]){
      arrows(x0=GELMAN.TIS[i,3],
             y0=GELMAN.TIS[i,1],
             x1=GELMAN.TIS[i,3],
             y1=GELMAN.TIS[i,2],
             code=3,length=0.04,angle=90,col='darkgray')
    }
    abline(h=1)
    text(x=seq(1,dim(GELMAN.TIS)[1],by=1),y=-max(GELMAN.TIS[,1:2])/2, srt = 00, adj= 0.5, xpd = TRUE,srt = 90,
         labels =paste("t.sim",1:dim(GELMAN.TIS)[1]), cex=0.65)
    legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")

    ###############
    par_to_use = grep('lambda',colnames(x$JAGS[[1]]))
    GELMAN.TIS <- matrix(NA, nrow=length(par_to_use), ncol=3)
    for (v in 1:length(par_to_use)) {
      GELMAN.TIS[v,1:2] <- gelman.diag(x$JAGS[,par_to_use[v]],multivariate=FALSE,confidence =CI,autoburnin=ABURNIN)$psrf
      GELMAN.TIS[v,3] <- v
    }
    GELMAN.TIS[is.na(GELMAN.TIS)]=1
    plot(GELMAN.TIS[,1], xlab="Estimates of lambda", xaxt='n',ylab="Factor",
         ylim=c(0,max(GELMAN.TIS[,1:2])),xlim=c(1,(dim(GELMAN.TIS)[1])))

    for(i in 1:dim(GELMAN.TIS)[1]){
      arrows(x0=GELMAN.TIS[i,3],
             y0=GELMAN.TIS[i,1],
             x1=GELMAN.TIS[i,3],
             y1=GELMAN.TIS[i,2],
             code=3,length=0.04,angle=90,col='darkgray')
    }
    abline(h=1)
    text(x=seq(1,dim(GELMAN.TIS)[1],by=1),y=-max(GELMAN.TIS[,1:2])/2, srt = 00, adj= 0.5, xpd = TRUE,srt = 90,
         labels =paste("lambda",1:dim(GELMAN.TIS)[1]), cex=0.65)
    legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")

    par(mfrow=c(1,1))
  }
  if(type=='lambda') {
    lambda=x$lambda[BURN:dim(x$lambda)[1],]
    lambda.lwr=apply(lambda,2,'quantile',(1-CI)/2)
    lambda.med=apply(lambda,2,'quantile',0.5)
    lambda.upr=apply(lambda,2,'quantile',1-(1-CI)/2)
    lambdaNames=names(lambda.med)

    ti.sim=apply(x$ti.sim[BURN:dim(x$ti.sim)[1],],2,'quantile',0.5)

    #plot(x=ti.sim,y=lambda.med,type='l', xlab='Sim ID',ylab='lambda',main='lambda')

    MAIN=paste("Realizations of lambda with CI= ", CI*100,"%",sep='')
    plot(y=lambda.med,x=ti.sim,type='l',main=MAIN,ylab="Level",xlab='time',lwd=1,lty=1,
         ylim=c(min(lambda),1.2*max(lambda)), xlim=c(min(x$ti.mu,ti.sim),max(x$ti.mu,ti.sim)))

    lines(y=lambda.lwr,x=ti.sim,type='l',col='blue',lwd=1,lty=2)
    lines(y=lambda.upr,x=ti.sim,type='l',col='blue',lwd=1,lty=2)

    legend("topright",legend = c("Upper CI","Medium","Lower CI"),
           col=c("blue","black","blue"),lwd=c(1,1,1),lty=c(2,1,2))

  }
  if(type=='mcmc') {
    mcmcplot(x$JAGS)
  }
}

