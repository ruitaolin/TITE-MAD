
#' Obtain the operating characteristics of the Time-to-event model-assisted design for single agent trials with late-onset
#' toxicities by simulating trials.
#'
#'
#' @param target the target toxicity rate
#' @param p.true a vector containing the true toxicity probabilities of the
#'              investigational dose levels.
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param maxt the maximum follow-up time 
#' @param prior.p a vector of length 3, which specifies the prior probability that the time to toxicity lies
#'              inside the time interval (0,\code{maxt}/3), (\code{maxt}/3,\code{2maxt}/3), (\code{2maxt}/3,1).
#'              The default value is \code{prior.p=c(1/3,1/3,1/3)}. 
#' @param accrual the accrual rate, i.e., the number of patients accrued in 1 unit of time
#' @param dist1 the underlying distribution of the time to toxicity outcomes; \code{dist1=1} is the
#'              uniform distribution, \code{dist1=2} corresponds to the Weibull distribution,
#'              \code{dist1=3} is the log-logistic distribution.
#' @param dist2 the underlying distribution of patient arrival time; \code{dist2=1} is the 
#'              uniform distribution, \code{dist2=2} is the exponential distribution
#' @param alpha a number from (0,1) that controls alpha*100% events in (0, 1/2T). 
#'              The default is \code{alpha=0.5}.             
#' @param n.earlystop the early stopping parameter. If the number of patients
#'                    treated at the current dose reaches \code{n.earlystop},
#'                    stop the trial and select the MTD based on the observed data.
#'                    The default value \code{n.earlystop=100} essentially turns
#'                    off this type of early stopping.
#' @param startdose the starting dose level for the trial
#' @param p.saf the lower bound of the proper dosing interval, e.g., \code{p.saf=target-0.05}
#' @param p.tox the upper bound of the proper dosing interval, e.g., \code{p.saf=target-0.05}
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset a small positive number (between 0 and 0.5) to control how strict the
#'               stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#'               strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param ntrial the total number of trials to be simulated.
#' 
#' @details This function generates he operating characteristics of the Time-To-Event Bayesian Optimal 
#'          Interval (TITE-BOIN) design for trials 
#'          with delayed toxicities by simulating trials under the prespecified true toxicity 
#'          probabilities of the investigational doses. Various O/I ratios, underlying distributions 
#'          of time to toxicity and time to arrival can be specified. 
#' @return \code{get.oc()} returns the operating characteristics of the TITE-BOIN design as a data frame,
#'         including: (1) selection percentage at each dose level (\code{selpercent}),
#'         (2) the number of patients treated at each dose level (\code{nptsdose}),
#'         (3) the number of toxicities observed at each dose level (\code{ntoxdose}),
#'         (4) the average number of toxicities (\code{totaltox}),
#'         (5) the average number of patients (\code{totaln}),
#'         (6) the percentage of early stopping without selecting the MTD (\code{pctearlystop}).   
#'         (7) the average number of suspending times (\code{npend})
#'         (8) the average trial duration needed for the trial based on the TITE-BOIN design (\code{duration})    
#' @export
#' @note We should avoid setting the values of \code{p.saf} and \code{p.tox} very close to the \code{target}.
#'       This is because the small sample sizes of typical phase I trials prevent us from
#'       differentiating the target toxicity rate from the rates close to it. In addition,
#'       in most clinical applications, the target toxicity rate is often a rough guess,
#'       and finding a dose level with a toxicity rate reasonably close to the target rate
#'       will still be of interest to the investigator. In addition, we recommend setting 
#'       the value of \code{priortox} relatively small, for example, \code{priortox=target/2} to accelerate
#'       the escalation procedure. The default values provided by  \code{get.iboundary()} 
#'       are generally reasonable for most clinical applications.
#'               
#' @example 
#' target<-0.3
#' p.true<-c(0.14,0.30,0.39,0.48,0.56,0.64,0.70)
#' get.oc.tite(target,p.true,12,3,rate=6,ntrial=500)


get.oc.tite <- function(target, p.true, ncohort, cohortsize, maxt=1, prior.p=NA, accrual=3, maxpen=0.5, 
                        dist1=1, dist2=1,alpha=0.5,n.earlystop = 100, startdose = 1, 
                        p.saf = target-0.05, p.tox = target+0.05, cutoff.eli = 0.95, 
                        extrasafe = FALSE, offset = 0.05, ntrial = 1000, seed=123, design=1)
{
  
  select.mtd <- function(target, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05, verbose = TRUE) {
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function(x, wt = rep(1, length(x))) {
      n <- length(x)
      if (n <= 1)
        return(x)
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol)))
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    ## determine whether the dose has been eliminated during the trial
    y = ntox
    n = npts
    ndose = length(n)
    elimi = rep(0, ndose)
    for (i in 1:ndose) {
      if (n[i] >= 3) {
        if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) > cutoff.eli) {
          elimi[i:ndose] = 1
          break
        }
      }
    }
    
    if (extrasafe) {
      if (n[1] >= 3) {
        if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) > cutoff.eli - offset) {
          elimi[1:ndose] = 1
        }
      }
    }
    
    ## no dose should be selected (i.e., selectdose=99) if the first dose is already very toxic or all uneliminated doses are never used to treat patients
    if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
      selectdose = 99
    } else {
      adm.set = (n != 0) & (elimi == 0)
      adm.index = which(adm.set == T)
      y.adm = y[adm.set]
      n.adm = n[adm.set]
      
      ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
      phat = (y.adm + 0.05)/(n.adm + 0.1)
      phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm + 0.1)^2 * (n.adm + 0.1 + 1))
      
      ## perform the isotonic transformation using PAVA
      phat = pava(phat, wt = 1/phat.var)
      phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number
      selectd = sort(abs(phat - target), index.return = T)$ix[1]  ## select dose closest to the target as the MTD
      selectdose = adm.index[selectd]
    }
    
    if (verbose == TRUE) {
      if (selectdose == 99) {
        out = list(target = target, MTD = selectdose,
                   p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=rep("----",length(npts)),
                                            'CI'=paste("(", rep("----",length(npts)),",",rep("----",length(npts)),")",sep="")))
        )
      } else {
        trtd = (n != 0)
        poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] - y[trtd] + 0.05))
        phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] + 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 * (n[trtd] + 0.1 + 1))))
        
        A1 = A2 = NA
        A3 = NA
        A4 = A5 = NA
        ## output summary statistics
        for (i in 1:ndose) {
          if (n[i] > 0) {
            A1 = append(A1, formatC(phat.all[i], digits = 2, format = "f"))
            A2 = append(A2, formatC(qbeta(0.025, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
            A3 = append(A3, formatC(qbeta(0.975, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
            A4 = append(A4, formatC(poverdose[i], digits = 2, format = "f"))
          } else {
            # no estimate output for doses never used to treat patients
            A1 = append(A1, "----")
            A2 = append(A2, "----")
            A3 = append(A3, "----")
            A4 = append(A4, "----")
          }
        }
        p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=A1[-1], 'CI'=paste("(", A2[-1],",",A3[-1],")",sep="")))
        
        out = list(target = target, MTD = selectdose, p_est=p_est, p_overdose = A4[-1])
      }
    } else {
      out = list(target = target, MTD = selectdose)
    }
    return(out)
  }
  
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
  
  
  ### simple error checking
  if(target<0.05) {cat("Error: the target is too low! \n"); return();}
  if(target>0.6)  {cat("Error: the target is too high! \n"); return();}
  if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return();}
  if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return();}
  if(p.saf<0.05) {cat("Error: the lower interval boundary cannot be too close to 0! \n"); return();}
  if(p.tox>0.95) {cat("Error: the upper interval boundary cannot be too close to 1! \n"); return();}
  if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
  if(!is.na(prior.p[1])){if(length(prior.p)!=3){cat("Error: The length of the prior probabilities should be 3! \n"); return();}}
  if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}
  if(is.na(maxpen)){maxpen=0.5;}
  if(maxpen<0 || maxpen>0.65) {cat("Error: the value of maxpen should lie within (0,0.65]!  \n"); return();}
  epi1<-target-p.saf;
  epi2<-p.tox-target;
  
    #set.seed(seed);
    if(is.na(prior.p[1])){prior.p = rep(1/3,3)}
    prior.p = prior.p/sum(prior.p)
    ndose = length(p.true);
    Y = matrix(rep(0, ndose * ntrial), ncol = ndose);
    N = matrix(rep(0, ndose * ntrial), ncol = ndose);
    dselect = rep(0, ntrial);
    durationV = rep(0, ntrial);
  	npendV = rep(0, ntrial);
    npts = ncohort*cohortsize;
    a<-b<-1

    for(trial in 1:ntrial)
    {
        y=NULL;  #toxicity indicator for each subject
        dv=NULL;  #dose for each subject
        n.d = rep(0, ndose);  # number of toxicity at each dose
        y.d = rep(0, ndose);  # number of patient at each dose
        t.enter=NULL; # time enter the study
        t.event=NULL; # time to event
        t.decision = 0; # decision making time
        d = startdose;  # current dose level
        earlystop = 0; #indicate if trial stops early
        elimi = rep(0, ndose)
	    	npend = 0;
        
        for(i in 1:ncohort)
        {
            # generate data for the new patient
            for(j in 1:cohortsize)
            {
                if(j==1) { t.enter = c(t.enter, t.decision); }
                else { 
                  if(dist2==1){ t.enter = c(t.enter, t.enter[length(t.enter)] + runif(1, 0, 2/accrual))}
                  if(dist2==2){ t.enter = c(t.enter, t.enter[length(t.enter)] + rexp(1, rate=accrual))}
                }
            }
            obscohort = gen.tite(dist1, cohortsize, p.true[d], alpha=alpha,T=maxt);
            t.event = c(t.event, obscohort$t.tox);
            y = c(y, obscohort$tox);
            dv = c(dv, rep(d, cohortsize));
            t.decision = t.enter[length(t.enter)];
            nobs=-1; pending=1;
            d.curr=d;
            npend = npend-1;
            while(pending==1)
            {
                npend = npend+1;
			         	pending = 0;
                if(i==ncohort) { t.decision = t.decision + maxt; }
                else {
                  if(dist2==1){t.decision = t.decision + runif(1, 0, 2/accrual)}
                  if(dist2==2){t.decision = t.decision + rexp(1, rate=accrual)}
                  }
			         	
                # determine which observation are observed
                delta = ((t.enter+t.event)<=t.decision);
                t = pmin(t.event, t.decision-t.enter, maxt);  ## used for recording potential censoring time
                cset = (dv==d);
                delta.curr = delta[cset];
                t.curr = t[cset];
                ntox.curr = sum((y[cset])[delta.curr==1]);
                #totalt = sum(t.curr[delta.curr==0])/maxt;
                totalt = t.curr[delta.curr==0]
                totalt = 3*prior.p[1]*totalt*(totalt<=maxt/3)+
                  ((prior.p[1]-prior.p[2])*maxt+3*prior.p[2]*totalt)*(maxt/3<totalt & totalt<=2*maxt/3)+
                  ((prior.p[1]+prior.p[2]-2*prior.p[3])*maxt+3*prior.p[3]*totalt)*(2*maxt/3<totalt & totalt<=maxt)
                totalt = sum(totalt)/maxt
                n.curr = sum(cset);
                n.pend = sum(delta[cset]==0);
                nobs = sum(delta[cset]);
                
                # determine which dose level should be eliminated
                for(dd in 1:ndose){
                  cset1 = dv==dd;
                  delta.curr1 = delta[cset1];
                  ntox.curr1 = sum((y[cset1])[delta.curr1==1]);
                  n.curr1 = sum(cset1);
                  if (1-pbeta(target, ntox.curr1+1, n.curr1-ntox.curr1+1)>cutoff.eli && n.curr1>=3){
                    elimi[dd:ndose]=1;
                    break;
                  }
                }
                
                #check whether extra safey rule should be applied
                if(extrasafe)
                {
                  if(d==1){
                    if(1-pbeta(target, ntox.curr+1, n.curr-ntox.curr+1)>cutoff.eli-offset && n.curr>=3) {
                      earlystop = 1;   break;}
                  }
                }
                
                # check whether the current dose level should be eliminated
                if(elimi[d.curr]==1) {
                  d=which(elimi==1)[1]-1
                  if(d==0){earlystop = 1; break;}
                  next;
                }
                
                #check whether the trial should be early terminated
                if(n.curr>=n.earlystop){break;}

                #check whether the current dose is toxic based on observed ata
                q0<-pbeta(target+epi2,ntox.curr+a,n.curr-ntox.curr+b)-pbeta(target-epi1,ntox.curr+a,n.curr-ntox.curr+b)
                q1<-pbeta(target-epi1,ntox.curr+a,n.curr-ntox.curr+b)-pbeta(target-2*epi1-epi2,ntox.curr+a,n.curr-ntox.curr+b)
                q2<-pbeta(target+epi1+2*epi2,ntox.curr+a,n.curr-ntox.curr+b)-pbeta(target+epi2,ntox.curr+a,n.curr-ntox.curr+b)
                if(q2>q1 & q2>q0){
                  if(d==1){d=d; if(n.pend>0){pending=1};} else{d=d-1};
                  next;
                }
                
                #check whether the trial should be suspended
                if(n.pend<0) {pending=1;}
                else
                {
                  y0<-ntox.curr
                  n0<-n.curr-n.pend+totalt
                  if(design==2){
                    q0<-pbeta(target+epi2,y0+1,n0-y0+1)-pbeta(target-epi1,y0+1,n0-y0+1)
                    q0<-q0/(epi2+epi1)
                    q1<-pbeta(target-epi1,y0+1,n0-y0+1)-pbeta(0,y0+1,n0-y0+1)
                    q1<-q1/(target-epi1)
                    q2<-pbeta(1,y0+1,n0-y0+1)-pbeta(target+epi2,y0+1,n0-y0+1)
                    q2<-q2/(1-target-epi2)
                  } else {
                    q0<-pbeta(target+epi2,y0+a,n0-y0+b)-pbeta(target-epi1,y0+a,n0-y0+b)
                    q1<-pbeta(target-epi1,y0+a,n0-y0+b)-pbeta(target-2*epi1-epi2,y0+a,n0-y0+b)
                    q2<-pbeta(target+epi1+2*epi2,y0+a,n0-y0+b)-pbeta(target+epi2,y0+a,n0-y0+b)
                  }

                  
                  #make dose assignment decisions
                  if(q1>q0 & q1>q2){
                    if(nobs<=2){
                      if(n.curr==1){d=d;} else {if(n.pend==0){d=min(d+1,ndose)} else {pending=1;}}
                    } else {
                      #check whether the current dose is the highest
                      if(d==ndose){d=d;} else{
                        if(elimi[d+1]==1){d=d} else{ d=d+1}
                      }}} else if (q2>q0 & q2>q1){
                        if(d==1){d=d} else {d=d-1}
                      } else {d=d;
                      }
                  if(elimi[d]==1){d=which(elimi==1)[1]-1}
                }
            }
            if(earlystop==1){break;}
        }
        
	    	for(k in 1:ndose){
	    	  y.d[k] = sum(y[dv==k]);
	    	  n.d[k] = sum(dv==k);
	    	}
	    	
	    	npendV[trial]= npend;
        Y[trial, ] = y.d
        N[trial, ] = n.d
        durationV[trial] = t.decision
        if (earlystop == 1) {
            dselect[trial] = 99
        }
        else {dselect[trial] = select.mtd(target, n.d, y.d, cutoff.eli, extrasafe, offset, verbose=FALSE)$MTD}
    }

    selpercent = rep(0, ndose)
    selpercent=rep(0, ndose);
    nptsdose = apply(N,2,mean);
    ntoxdose = apply(Y,2,mean);
    
    for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100; }
    
    if(length(which(p.true==target))>0) # if MTD exists, calculate risk of overdosing
    {
      if (which(p.true==target) == ndose-1) {
        overdosing60=mean(N[,p.true>target]>0.6*npts)*100;
        overdosing80=mean(N[,p.true>target]>0.8*npts)*100;
      } else {
        overdosing60=mean(rowSums(N[,p.true>target])>0.6*npts)*100;
        overdosing80=mean(rowSums(N[,p.true>target])>0.8*npts)*100;
      }
      
      out=list(selpercent=selpercent, npatients=nptsdose, ntox=ntoxdose, totaltox=sum(Y)/ntrial, totaln=sum(N)/ntrial,
               percentstop=sum(dselect== 99)/ntrial*100, # poorallocation=mean(N[, p.true==target]<npts/ndose)*100,
               overdose60=overdosing60, overdose80=overdosing80, duration=mean(durationV),sdduration=sqrt(var(durationV)),simu.setup=data.frame(target=target, p.true=p.true, ncohort=ncohort, cohortsize = cohortsize,
                                                                                       startdose = startdose,p.saf = p.saf, p.tox = p.tox, cutoff.eli = cutoff.eli, extrasafe = extrasafe, offset = offset,
                                                                                       ntrial = ntrial, dose=1:ndose),flowchart=TRUE);
    }
    else {
      out=list(selpercent=selpercent, npatients=nptsdose, ntox=ntoxdose, totaltox=sum(Y)/ntrial, totaln=sum(N)/ntrial,
               percentstop=sum(dselect== 99)/ntrial*100, duration=mean(durationV),sdduration=sqrt(var(durationV)),simu.setup=data.frame(target=target, p.true=p.true, ncohort=ncohort, cohortsize = cohortsize,
                                                                               startdose = startdose,p.saf = p.saf, p.tox = p.tox, cutoff.eli = cutoff.eli, extrasafe = extrasafe, offset = offset, ntrial = ntrial,
                                                                               dose=1:ndose),flowchart=TRUE);
    }
    return(out);

}


target<-0.3
P<-matrix(0,6,6)

## Six fixed scenarios in the paper 
P[1,]<-c(0.13,0.28,0.41,0.50,0.60,0.70)
P[2,]<-c(0.08,0.15,0.29,0.43,0.50,0.57)
P[3,]<-c(0.28,0.42,0.49,0.61,0.76,0.87)
P[4,]<-c(0.05,0.10,0.20,0.31,0.50,0.70)
P[5,]<-c(0.06,0.08,0.12,0.18,0.30,0.41)
P[6,]<-c(0.05,0.06,0.08,0.11,0.19,0.32)

## TITE-keyboard design  
get.oc.tite(target,P[1,],ncohort=12,cohortsize=3,accrual=2,ntrial=5000,maxt=3)

get.oc.tite(target,P[5,],ncohort=12,cohortsize=3,accrual=2,ntrial=5000,maxt=3)


## TITE-mTPI design

get.oc.tite(target,P[1,],ncohort=12,cohortsize=3,accrual=2,ntrial=5000,maxt=3, design=2)

