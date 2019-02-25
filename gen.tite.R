#'
#' Generate time to toxicity based on unconditional model (uniform, Weibull and log-logistic)
#'
#' @usage gen.tite(dist=1, n, p, alpha=0.5, T=1)
#'
#' @param dist the underlying distribution for the time-to-event data; \code{dist=1} is the
#'             uniform distribution, \code{dist=2} corresponds to the Weibull distribution,
#'             \code{dist=3} is the log-logistic distribution.
#' @param n the number of samples required to be simulated
#' @param pi the toxicity probability
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param alpha a number from (0,1) that controls alpha*100% events in (0, 1/2T). 
#'              The default is \code{alpha=0.5}.
#' @param T the observation period
#' @details This function is used to generate time to toxicity data for dose-finding trials with
#'          late-onset toxicities. The distributions for the time to toxicity data include uniform,
#'          Weibull, and log-logistic distributions. In addition, by spcifying \code{alpha}, users can
#'          control the degree of late-onset. 
#'
#' @return a vector of time to toxicity data
#' @examples
#' gen.tite(dist=1, 5, 0.3, alpha=0.5, T=1)

gen.tite<-function(dist=1, n, pi, alpha=0.5, Tobs=1)
{
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
        ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
        alpha = log(log(1-pi)/log(1-pihalft))/log(2);
        lambda = -log(1-pi)/(Tobs^alpha);
        t = (-log(runif(n))/lambda)^(1/alpha);
        return(t);
    }
    
    llogit<-function(n, pi, pihalft)
    {
        ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
        alpha = log((1/(1-pi)-1)/(1/(1-pihalft)-1))/log(2);
        lambda = (1/(1-pi)-1)/(Tobs^alpha);
        t = ((1/runif(n)-1)/lambda)^(1/alpha);
        return(t);
    }
    ############ end of subroutines ############
    
    
    tox = rep(0, n);
    t.tox = rep(0, n);
    
    #### uniform
    if(dist==1) {  # 50% event in (0, 1/2T)
        tox = rbinom(n, 1, pi);
        ntox.st = sum(tox);
        t.tox[tox==0]=Tobs;
        t.tox[tox==1]=runif(ntox.st, 0, Tobs);
    }
    #### Weibull
    if(dist==2)
    {
        pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
        t.tox = weib(n, pi, pihalft);
        tox[t.tox<=Tobs]=1;
        ntox.st = sum(tox);
        t.tox[tox==0]=Tobs;
    }
    #### log-logistic
    if(dist==3)
    {
        pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
        t.tox = llogit(n, pi, pihalft);
        tox[t.tox<=Tobs]=1;
        ntox.st = sum(tox);
        t.tox[tox==0]=Tobs;
    }
    return(list(tox=tox, t.tox=t.tox, ntox.st=ntox.st));
}
