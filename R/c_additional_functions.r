#' JAGS output
#'
#' \code{JAGS.objects(JAGS.output)} an internal function used within the tuts package. Additional
#'  plots and diagnostics are available from model-specific methods, for more information see the
#'  help documents of individual models.
#'
#' @param JAGS.output An JAGS MCMC Object.
JAGS.objects=function(JAGS.output){
  chain=1
  Output.JAGS=list()
  JAGS.names=varnames(JAGS.output[[chain]])
  JAGS.objects=unique(gsub("\\[.*","",JAGS.names))
  for (i in 1:length(JAGS.objects)){
    ncol=sum((gsub("\\[.*","",JAGS.names))==JAGS.objects[i])
    objects=paste("JAGS.",JAGS.objects[i],"=matrix(data=NA,nrow=dim(JAGS.output[[1]])[1],ncol=ncol)",sep="")
    eval(parse(text=objects))
    object.names=paste("JAGS.",JAGS.objects[i],".names=rep(NA,ncol)",sep="")
    eval(parse(text=object.names))
  }
  for (Obj in JAGS.objects){
    j=1
    for (i in 1:length(JAGS.names)){
      if (gsub("\\[.*","",JAGS.names[i])==Obj){
        eval(parse(text=paste("JAGS.",Obj,".names[j]=JAGS.names[i]",sep="")))
        eval(parse(text=paste("JAGS.",Obj,"[,j]=JAGS.output[[1]][,JAGS.names[i]]",sep="")))
        j=j+1
      }
    }
    eval(parse(text=paste("colnames(JAGS.",Obj,")=JAGS.",Obj,".names",sep="")))
  }
  for (Obj in JAGS.objects){
    Output.JAGS[[(length(Output.JAGS)+1)]]=eval(parse(text=paste("JAGS.",Obj,sep="")))
  }
  # add optput
  names(Output.JAGS)=JAGS.objects
  return(Output.JAGS)
}

##############################################################################################
#' Wrapper of models in the tuts package
#'
#' \code{tuwrap=function(y,ti.mu,ti.sd,n.sim, ...)} tuwrap compares results obtained from
#'  fitting multiple of time-uncertain time series models.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates/observed timings of observations.
#' @param ti.sd A vector of standard deviations of timings.
#' @param n.sim A number of simulations.
#' @param ... optional arguments. n.chains: number of MCMC chains, default is set to 2.
#' Thin: Thinning factor, default is set to 4.m: maximum number of significant frequencies
#'  in the data, the default is set to 5. polyorder: order of the polynomial regression component,
#'  the default is set to 3. freqs: set to a positive integer k returns a vector of k equally spaced frequencies in the Nyquist
#'   range. freqs can be provided as a vector of custom frequencies of interest. Set to 'internal' generates an equally
#'   spaced vector of frequencies in the Nyquist range
#'
#' @examples
#' # Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=50,Harmonics=c(10,30,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),
#' trend=0,y.sd=3, ti.sd=1)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(1, length(ti.mu))
#'
#' Run simulations:
#' n.sim=1000
#' WRAP=tuwrap(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim)
#'
#' Compare the models:
#' summary(WRAP)
#'
#' # Individual summaries, plots and diagnostics:
#' # Note: Accessing of individual summaries, diagnostics and plots is presented in manuals of all of the funcions contained in the wrapper.SOme examples:
#'
#' plot(WRAP$BFS,type='periodogram')
#' summary(WRAP$BFS,tCI=0.99)
#'
#' @export
#'

tuwrap=function(y,ti.mu,ti.sd,n.sim, ...){
  dots = list(...)

  if (is.null(dots$n.chains)) {n.chains=2}else{n.chains=dots$n.chains}
  if (is.null(dots$CV)) {CV=FALSE}else{CV=dots$CV}
  if (is.null(dots$freqs)) {freqs="internal"}else{freqs=dots$freqs}
  if (is.null(dots$polyorder)) {polyorder=c(1:3)}else{polyorder=dots$polyorder}
  if (is.null(dots$type)) {
   if(!sum(y==round(y) & y>=0)==length(y)){
     type="continuous"
    } else {
      RESP=readline(prompt="Are the observations discrete? (y/n):")
      if (RESP=='y' | RESP=='Y'){
        type="counting"
      } else{
        type="continuous"
      }
    }
  }

  output=list()
  if (type=="continuous"){
    BFS=tubfs(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs=freqs,n.sim=n.sim,n.chains=n.chains, CV=CV)
    output$BFS=BFS

    TUAR1=tuar1(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=CV)
    output$TUAR1=TUAR1

    AR1REDF=tuar1redf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=CV)
    output$AR1REDF=AR1REDF

    for (i in 1: length(polyorder)){
      eval(parse(text = paste(
        "PN",polyorder[i],"=tupolyn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder[i],
            n.sim=n.sim,n.chains=n.chains, CV=CV)"
        ,sep='')))
      eval(parse(text = paste("output$PN",polyorder[i],"=PN",polyorder[i],sep='')))
    }


  } else if(type=="counting"){
    output$TUPOISBSF=tupoisbsf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=n.chains, CV=CV)
    output$TUPOISPN3=tupoispn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=CV)
    output$TUNBINOM=tunbinombsf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=n.chains, CV=CV)
  } else{
    stop("Type must be either a 'continuous' or 'counting' process.")
  }
  class(output)="tuts_wrap"
  return(output)
}

# =========
#' Summary table comparing models based on DIC criterion.
#'
#' \code{summary.tuts_wrap(x)} Summary table of a tuts_wrap objects.
#'
#' @export
#'
summary.tuts_wrap = function(x) {
  NAMES=names(x)

  DIC=array(NA,dim=c(length(NAMES),3))
  rownames(DIC)=NAMES; colnames(DIC)=c("Deviance","Penalty", "Penalized Dev")

  for(i in 1:length(NAMES)){
    Deviance=eval(parse(text = paste( "round(sum(x$",NAMES[i],"$DIC$deviance),2)",sep='')))
    Penalty=eval(parse(text = paste( "round(sum(x$",NAMES[i],"$DIC$penalty),2)",sep='')))
    DIC[i,1]=Deviance
    DIC[i,2]=Penalty
    DIC[i,3]=Deviance+Penalty
  }
  print(DIC)
}
