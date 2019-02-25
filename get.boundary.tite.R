#'
#' Generate the optimal dose escalation and deescalation boundaries for conducting the trial with late-onset toxicities.
#'
#' Use this function to generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#' @param target the target toxicity rate
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param n.earlystop the early stopping parameter. If the number of patients treated at
#'                    the current dose reaches \code{n.earlystop}, stop the trial
#'                    and select the MTD based on the observed data. the default
#'                    value \code{n.earlystop=100} essentially turns off the type
#'                    of early stopping.
#'               
#' @param prior.p a vector of length 3, which specifies the prior probability that the time to toxicity lies
#'              inside the time interval (0,T/3), (T/3,2T/3), (2T/3,1).
#'              The default value is \code{prior.p=c(1/3,1/3,1/3)}.        
#' @param p.saf the lower bound of the proper dosing interval, e.g., \code{p.saf=target-0.05}
#' @param p.tox the upper bound of the proper dosing interval, e.g., \code{p.saf=target-0.05}
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                   We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more strict stopping rule for extra safety
#' @param offset a small positive number (between 0 and 0.5) to control how strict
#'               the stopping rule is when \code{extrasafe=TRUE}. A larger value leads
#'               to a more strict stopping rule. the default value \code{offset=0.05}
#'               generally works well.
#' @param print to print out the boundary results.
#' @param output to print the decision table to an excel file. 
#' @param design the design indicator: \code{design=1} means the TITE-keyboard design, \code{desing=2} means
#'               the TITE-mTPI design. 
#'
#'  Author: Ruitao Lin  ruitaolin@gmail.com 
#'          Ying Yuan yyuan@mdanderson.org
#'          
#'              
get.boundary.tite <- function(target, ncohort, cohortsize, n.earlystop=100,prior.p=NA,
                          p.saf = target-0.05, p.tox = target+0.05,
                          cutoff.eli=0.95,extrasafe=FALSE, offset=0.05, print=TRUE, output=TRUE, design=1)
{
  ### simple error checking
  if(target<0.05) {cat("Error: the target is too low! \n"); return();}
  if(target>0.6)  {cat("Error: the target is too high! \n"); return();}
  if((target-p.saf)<(0.1*target)) {cat("Error: the lower interval boundary cannot be higher than or too close to the target! \n"); return();}
  if((p.tox-target)<(0.1*target)) {cat("Error: the upper interval boundary cannot be lower than or too close to the target! \n"); return();}
  if(p.saf<0.05) {cat("Error: the lower interval boundary cannot be too close to 0! \n"); return();}
  if(p.tox>0.95) {cat("Error: the upper interval boundary cannot be too close to 1! \n"); return();}
  if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
  if(!is.na(prior.p[1])){if(length(prior.p)!=3){cat("Error: The length of the prior probabilities should be 3! \n"); return();}}
  if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}

  epi1<-target-p.saf;
  epi2<-p.tox-target;
  
  if(design==2){
    boundf<-function(n.curr,ntox.curr,n.pend,totalt,target,epi1,epi2,type=1){
      y0<-ntox.curr
      n0<-n.curr-n.pend+totalt
      if(type==1){
        q0<-(pbeta(target+epi2,y0+1,n0-y0+1)-pbeta(target-epi1,y0+1,n0-y0+1))/(epi1+epi2)
        q1<-pbeta(target-epi1,y0+1,n0-y0+1)/(target-epi1)
        return(q1-q0)
      }
      if(type==2){
        q0<-(pbeta(target+epi2,y0+1,n0-y0+1)-pbeta(target-epi1,y0+1,n0-y0+1))/(epi1+epi2)
        q2<-(1-pbeta(target+epi2,y0+1,n0-y0+1))/(1-target-epi2)
        return(q2-q0)
      }
    }    
  } else {
    boundf<-function(n.curr,ntox.curr,n.pend,totalt,target,epi1,epi2,type=1){
      y0<-ntox.curr
      n0<-n.curr-n.pend+totalt
      if(type==1){
        q0<-pbeta(target+epi2,y0+1,n0-y0+1)-pbeta(target-epi1,y0+1,n0-y0+1)
        q1<-pbeta(target-epi1,y0+1,n0-y0+1)-pbeta(max(target-2*epi1-epi2,0),y0+1,n0-y0+1)
        return(q1-q0)
      }
      if(type==2){
        q0<-pbeta(target+epi2,y0+1,n0-y0+1)-pbeta(target-epi1,y0+1,n0-y0+1)
        q2<-pbeta(target+min(epi1+2*epi2,1),y0+1,n0-y0+1)-pbeta(target+epi2,y0+1,n0-y0+1)
        return(q2-q0)
      }
    }    
  }

  
    if(is.na(prior.p[1])){prior.p = rep(1/3,3)}
    prior.p = prior.p/sum(prior.p)
    boundary = NULL
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
    pattern=0;
    decision=NULL;
    nprior=1; # the effective sample size for the Beta prior.
    priortox=target/2; # the prior mean for the Beta prior.
    for(n in (1:ncohort)*cohortsize){
        lambda1<-uniroot(boundf,c(0,n),n.curr=n,n.pend=0,totalt=0,target=target,epi1=epi1,epi2=epi2,type=1)$root/n
        lambda2<-uniroot(boundf,c(0,n),n.curr=n,n.pend=0,totalt=0,target=target,epi1=epi1,epi2=epi2,type=2)$root/n 
        ntrt = c(ntrt, n);
        cutoff1 = floor(lambda1*n);
        cutoff2 = ceiling(lambda2*n);
        elimineed=0; # indicating whether elimination is needed
        if(n<3) { elim = c(elim, NA); j=NA;}  # require treating at least 3 patients before eliminating a dose
        else
        {
          for(j in 1:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
          {
            if(1-pbeta(target, j+1, n-j+1)>cutoff.eli) {elimineed=1; break;}
          }
          if(elimineed==1) { elim = c(elim, j); if (cutoff2>j) {cutoff2=j} }
          else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
        }
        
          for(ntox in 0:n){
            te.prev=-99;td.prev=-99
            if(ntox>=cutoff2) {
              elimineed=1;
              if(is.na(j)){
                decision=rbind(decision,
                               c(n,ifelse(ntox<n,paste(">",ntox-1),paste(ntox)),
                                 ifelse(n-ntox>0,paste("<",n-ntox+1),paste(n-ntox)),
                                 "No","No","Yes" ));
              } else {
                if(j==(ntox+1)){
                  decision=rbind(decision,
                                 c(n,ntox,
                                   ifelse(n-ntox>0,paste("<",n-ntox+1),paste(n-ntox)),
                                   "No","No","Yes" ));
                }
                if(j>(ntox+1)){
                  decision=rbind(decision,
                                 c(n,paste(as.character(seq(ntox,j-1)), collapse=", "),
                                   ifelse(n-ntox>0,paste("<",n-ntox+1),paste(n-ntox)),
                                   "No","No","Yes" ));
                }
                decision=rbind(decision,
                               c(n,ifelse(j<n,paste(">",j-1),paste(j)),
                                 ifelse(n-j>0,paste("<",n-j+1),paste(n-j)),
                                 "No","No","Yes & Eliminate"))
              }
              break;}
            for(npend in 0:(n-ntox)){
              if(npend==0){ 
                if((ntox/n)<=target){
                  if(boundf(n.curr=n,ntox.curr=ntox,n.pend=0,totalt=0,target=target,epi1=epi1,epi2=epi2,type=1)>0){
                    te=-Inf;td=-Inf
                  } else {te=Inf;td=-Inf}
                }
                if((ntox/n)>target){
                  if(boundf(n.curr=n,ntox.curr=ntox,n.pend=0,totalt=0,target=target,epi1=epi1,epi2=epi2,type=2)>0){
                    te=Inf;td=Inf
                  } else {te=Inf;td=-Inf}
                }
              } else {
                #compute the boundaries for the total follow-up proportion
                a1<-boundf(n.curr=n,ntox.curr=ntox,n.pend=npend,totalt=0,target=target,epi1=epi1,epi2=epi2,type=1)
                a2<-boundf(n.curr=n,ntox.curr=ntox,n.pend=npend,totalt=npend,target=target,epi1=epi1,epi2=epi2,type=1)
                a3<-boundf(n.curr=n,ntox.curr=ntox,n.pend=npend,totalt=0,target=target,epi1=epi1,epi2=epi2,type=2)
                a4<-boundf(n.curr=n,ntox.curr=ntox,n.pend=npend,totalt=npend,target=target,epi1=epi1,epi2=epi2,type=2)
                if(a1<0 & a2<0){te=Inf} else if (a1>0 & a2>0){te=-Inf} else{
                  te<-round(uniroot(boundf,c(0,npend),n.curr=n,ntox.curr=ntox,n.pend=npend,target=target,epi1=epi1,epi2=epi2,type=1)$root,2)
                }
                if(a3<0 & a4<0){td=-Inf} else if (a3>0 & a4>0) {td=-Inf} else{
                  td<-round(uniroot(boundf,c(0,npend),n.curr=n,ntox.curr=ntox,n.pend=npend,target=target,epi1=epi1,epi2=epi2,type=2)$root,2)
                }
                if(te>=npend) te=Inf;
                if(te<=0) te=-Inf;
                if(td<=0) td=-Inf;
                
                te <- te + (n-npend)
                td <- td + (n-npend)
              }
              
              
              if(npend>=2 & n<=3 & ntox==0){
                decision=rbind(decision,c(n,ntox,ifelse(npend<n,paste(">",npend-1),paste(npend)),"Suspend","Suspend","Suspend"))
                break;
              } else{
                if (te==-Inf && td==-Inf){
                  decision=rbind(decision,c(n,ntox,npend,"Yes","No","No"))
                } else if (te==Inf && td==-Inf){
                  decision=rbind(decision,c(n,ntox,npend,"No","Yes","No"))
                } else if (te==Inf && td==Inf){
                  decision=rbind(decision,c(n,ntox,npend,"No","No","Yes"))
                } else if (te!=-Inf && td==-Inf){
                  if(sum(prior.p==rep(1/3,3))==3){
                    decision=rbind(decision,c(n,ntox,npend,paste("ESS>",te),paste("ESS<=",te),"No"))
                  } else {decision=rbind(decision,c(n,ntox,npend,paste("WESS>",te),paste("WESS<=",te),"No"))}
                } else if (te==Inf && td!=Inf){
                  if(sum(prior.p==rep(1/3,3))==3){
                    decision=rbind(decision,c(n,ntox,npend,"No",paste("ESS>=",td),paste("ESS<",td)))
                  } else {decision=rbind(decision,c(n,ntox,npend,"No",paste("WESS>=",td),paste("WESS<",td)))}
                } else {
                  if(sum(prior.p==rep(1/3,3))==3){
                    decision=rbind(decision,c(n,ntox,npend,paste("ESS>",te),paste0("ESS in [",td,",",te,"]"),paste("ESS<",td)))
                  } else {decision=rbind(decision,c(n,ntox,npend,paste("WESS>",te),paste0("WESS in [",td,",",te,"]"),paste("WESS<",td)))}
                }
              }
              
              if(round(te,1)==round(te.prev,1) && round(td,1)==round(td.prev,1) && te==-Inf && td==-Inf){
                decision[nrow(decision)-1,]=c(n,ntox,paste("<",npend+1),"Yes","No","No")
                decision=decision[-nrow(decision),]
              } else if(round(te,1)==round(te.prev,1) && round(td,1)==round(td.prev,1) && td==Inf && te==Inf){
                decision[nrow(decision)-1,]=c(n,ntox,paste("<",npend+1),"No","No","Yes")
                decision=decision[-nrow(decision),]
              } else if(round(te,1)==round(te.prev,1) && round(td,1)==round(td.prev,1) && td==-Inf && te==Inf){
                decision[nrow(decision)-1,]=c(n,ntox,paste("<",npend+1),"No","Yes","No")
                decision=decision[-nrow(decision),]
              } else if(round(te,1)==round(te.prev,1) && round(td,1)==round(td.prev,1) && te!=-Inf && td==-Inf){
                if(sum(prior.p==rep(1/3,3))==3){
                  decision[nrow(decision)-1,]=c(n,ntox,paste(npend.old,"~",npend),paste("ESS>",te),paste("ESS<=",te),"No")
                } else {decision[nrow(decision)-1,]=c(n,ntox,paste(npend.old,"~",npend),paste("WESS>",te),paste("WESS<=",te),"No")}
                decision=decision[-nrow(decision),]
              } else if(round(te,1)==round(te.prev,1) && round(td,1)==round(td.prev,1) && te==Inf && td!=Inf){
                if(sum(prior.p==rep(1/3,3))==3){
                  decision[nrow(decision)-1,]=c(n,ntox,paste(npend.old,"~",npend),"No",paste("ESS>=",td),paste("ESS<",td))
                } else {                decision[nrow(decision)-1,]=c(n,ntox,paste(npend.old,"~",npend),"No",paste("WESS>=",td),paste("WESS<",td))}
                decision=decision[-nrow(decision),]
              } else if(round(te,1)==round(te.prev,1) && round(td,1)==round(td.prev,1)){
                if(sum(prior.p==rep(1/3,3))==3){
                  decision[nrow(decision)-1,]=c(n,ntox,paste(npend.old,"~",npend),paste("ESS>",te),paste0("ESS in [",td,",",te,"]"),paste("ESS<",td))
                } else {decision[nrow(decision)-1,]=c(n,ntox,paste(npend.old,"~",npend),paste("WESS>",te),paste0("WESS in [",td,",",te,"]"),paste("WESS<",td))}
                decision=decision[-nrow(decision),]
              } else {
                npend.old<-npend
              }
              
              te.prev = te;
              td.prev = td;
            }
            npend=npend+1;
            #if(npend<=n & npend>(n/2)){decision=rbind(decision,c(n,ntox,paste(">=",npend),"Suspend","Suspend","Suspend"))}
          }
        
        
    }
    #decision=noquote(decision)
    decision = decision[as.numeric(decision[,1])<=n.earlystop,]
    colnames(decision) =c ("#Patients", "#DLTs observed", "#Pending patients", "Ecalation     ", "Stay          ","De-escalation          ")

    rownames(decision)=1:dim(decision)[1]
    out=list()
      out=list(lambda_e = lambda1,lambda_d = lambda2, boundary_tab=decision)
    
    if(extrasafe) {
      stopbd=NULL;
      ntrt=NULL;
      for(n in 1:npts) {
        ntrt = c(ntrt, n);
        if(n<3) { stopbd = c(stopbd, NA); }  # require treating at least 3 patients before stop a trial
        else {
          for(ntox in 1:n){ #determine stopping boundary, prior beta(1,1) is used in beta-binomial model
            if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli-offset) {stopneed=1; break;}
          }
          if(stopneed==1) { stopbd = c(stopbd, ntox); }else { stopbd = c(stopbd, NA); }
        }
      }
      
      stopboundary = rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)];
      rownames(stopboundary) = c("Number of patients treated at the lowest dose  ", "Stop the trial if # of DLT >=        ");
      colnames(stopboundary) = rep("", min(npts, n.earlystop));
      out = c(out,list(target=target, cutoff=cutoff.eli-offset, stop_boundary=stopboundary))
    }

    return(noquote(out))
}


# Reproduce Table 1, target=30%, TITE-keyboard design
get.boundary.tite(target=0.3,ncohort=4,cohortsize=3, design=1)  

# Reproduce Table S1, target=30%, TITE-mTPI design
get.boundary.tite(target=0.3,ncohort=4,cohortsize=3, design=2)  
