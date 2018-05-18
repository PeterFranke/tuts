#' Generating time-uncertain time series
#'
#' \code{simtuts} A function generating time-uncertain time series. It returns two data frames
#' containing simulation of an actual process and its observations.\cr
#' The actual process consists of a sum of a constant, a linear trend, and three sine and three cosine functions, and its
#' observations are normally distributed \emph{y.obs~N(y.act, y.sd)}.\cr
#' Timing of simulated processes is modeled  as \emph{t.act~U(0,N)} and sorted in the ascending order.
#' Observations of timings are modeled  in two ways:
#' \enumerate{
#' \item  Normally distributed timing \emph{t.obs.norm~N(ti.act,ti.sd)},
#' sorted from the smallest to the largest value to ensure non-overlapping feature of observations,
#' \item Timing simulated with truncated normal distribution t.obs.tnorm~N(ti.act,ti.sd,....).}
#' Note: variability of timing can be substantially greater when the normal distribution is chosen,
#' the truncated distribution utilizes enforced limits applied in the midpoints of the actual timing.
#'
#' @param N A number of observations.
#' @param Harmonics A vector of three harmonics, typically integers.
#' @param sin.ampl A vector of three amplitudes of the sine terms.
#' @param cos.ampl vector of three amplitudes of the cosine terms.
#' @param trend A constant trend.
#' @param y.sd A standard deviation of observations.
#' @param ti.sd A standard deviation of estimates of timing.
#'
#' @examples
#' # 1. Generate actual and observed time series as a sum of 2 sine functions:
#' DATA=simtuts(N=50,Harmonics=c(10,20,0), sin.ampl=c(10,10, 0), cos.ampl=c(0,0,0),trend=0,
#' y.sd=2, ti.sd=0.3)
#'
#' @references  \url{https://en.wikipedia.org/wiki/Truncated_normal_distribution}
#' @import truncnorm
#' @import rjags
#' @import coda
#' @import stats
#' @import mcmcplots
#' @import doParallel
#' @import foreach
#' @import lomb
#' @import parallel
#' @export
simtuts=function(N,Harmonics,sin.ampl,cos.ampl,trend=0, y.sd,ti.sd){
  if((round(N)==N)==FALSE | N<1){stop("N must be a positive integer")}
  if(is.numeric(trend)==FALSE | length(trend)>1){stop("trend must be a real number.")}
  if(is.numeric(y.sd)==FALSE | length(y.sd)>1 | y.sd<0){stop("y.sd must be a positive rational number.")}
  if(is.numeric(ti.sd)==FALSE | length(ti.sd)>1 | ti.sd<0){stop("ti.sd must be a positive rational number.")}
# Actual process ----------------------------------------------------------
  ti.act=sort(runif(N,0,1))
  y.act=(trend*ti.act*N)
  for (i in 1:length(Harmonics)) {
    data= sin.ampl[i]*sin(2*pi*ti.act*Harmonics[i])+cos.ampl[i]*cos(2*pi*ti.act*Harmonics[i])
    y.act=y.act +data
  }
  DATA.actual=data.frame(ti.act=ti.act*N,y.act=y.act)
# Observed data -----------------------------------------------------------
  y.obs=rnorm(N,y.act,y.sd)
  # using normally distributed timing
  ti.obs.norm=rnorm(n=N,mean = ti.act*N, sd=ti.sd)
  ti.obs.norm=sort(ti.obs.norm)-(min(ti.obs.norm))
  # using truncated normal distribution  timing
  DIFFs=diff(ti.act)
  ti=ti.upper=ti.lower=rep(NA,N)
  ti.upper[1:(N-1)]=ti.act[1:(N-1)]+DIFFs/2;ti.upper[N]=ti.act[N]+DIFFs[N-1]/2
  ti.lower[2:N]=ti.act[2:N]-DIFFs/2; ti.lower[1]=ti.act[1]-DIFFs[1]/2
  ti.obs.tnorm=rtruncnorm(1,a=ti.lower,b=ti.upper,mean=ti.act,sd=ti.sd)*N
  # DATA object
  DATA.observed=data.frame(y.obs=y.obs,ti.obs.norm=ti.obs.norm,ti.obs.tnorm=ti.obs.tnorm)
# Output --------------------------------------------------------------------
  OUTPUT= list(observed=DATA.observed, actual=DATA.actual)
  return(OUTPUT)
}
