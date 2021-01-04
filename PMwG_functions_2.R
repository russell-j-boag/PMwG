###Simulator###
####Simulate data for the PMwG Sampler with LBA####
rm(list=ls())
setwd("~/Documents/Research/Modelling Project/Scripts/recovery")
library(rtdists)
library(mvtnorm) ## For the multivariate normal.
library(MASS) ## For matrix inverse.
library(MCMCpack)
library(lme4)


drift = FALSE
threshold = FALSE
non.decision = FALSE
coherence = TRUE

n.trials = 500              #number trials per subject per conditions
n.subj = 100                 #number of subjects
conds = c("1", "2","3","4")   #cond names
n.cond = length(conds)      #number of conditions

names=c("subject","rt","response","condition", "stimulus") #names of columns
data = data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond)))
names(data)=names
data$condition = rep(conds,each = n.trials)
data$condition <- as.factor(data$condition)
data$subject = rep(1:n.subj, each = n.trials*n.cond)
data$stimulus = as.factor(rep(c("left","right"), length.out=length(data$subject)))

#mu names
if ( drift == TRUE){
  parameter.names=c("b", "A","vC.easy", "vC.hard","vE.easy","vE.hard","t0")
}
if (threshold == TRUE){
  parameter.names=c("b.easy", "b.hard", "A","vC","vE","t0")
}
if (non.decision== TRUE){
  parameter.names=c("b", "A","vC","vE","t0.easy", "t0.hard")
}
if (coherence == TRUE){
  parameter.names=c("b", "A","vC.1","vC.2","vC.3","vC.4","vE","t0")
}



n.parameters=length(parameter.names)
ptm <- array(dim = n.parameters, dimnames = list(parameter.names))
pts2 <- array(dim=c(n.parameters,n.parameters), dimnames = list(parameter.names,parameter.names))

#mu values
if (drift == TRUE){
  ptm[1:n.parameters]=log(c(0.3,0.4,4,2.5,1,2,0.1)) # (vs) 
}
if (threshold == TRUE){
  ptm[1:n.parameters]=log(c(0.3,0.6,0.3,2.5,0.9,0.1)) # (bs)
}
if (non.decision== TRUE){
ptm[1:n.parameters]=log(c(0.3,0.4,2.5,0.5,0.2,0.1)) # (t0s)
}
if (coherence== TRUE){
  ptm[1:n.parameters]=log(c(0.3,0.4,4,3,2,1,0.5,0.1)) 
}



pts2=diag(rep(.01,n.parameters)) ## Who knows??
pts2.inv=ginv(pts2) ## Because this is calculated near the end of the main loop, needs initialising for iter=1.



# bs
fast_loglike_b=function(x,data,sample=TRUE) {
  x=exp(x)
  bs=x["A"]+x[c("b.easy", "b.hard")][data$condition]## This is faster than "paste".
  if (sample) {
    for (stim in levels(data$stimulus)){
      isSR = data$stimulus==stim
    if (stim == "left"){
      v1 = x["vC"]
      v2 = x["vE"]
    } else {
      v1 = x["vE"]
      v2 = x["vC"]
    }
    out=rLBA(n=sum(isSR),A=x["A"],b=bs,t0=x["t0"],mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=TRUE)
    data$response[isSR]=c("left","right")[out$response]
    data$rt[isSR]=out$rt
    out<-data
    }
    }else {
      v1 = v2 = numeric(length(data$RT))
      for (stim in levels(data$stimulus)){
        isSR = data$stimulus==stim
        if (stim == "left"){
          v1[isSR] = x["vC"]
          v2[isSR] = x["vE"]
        } else {
          v1[isSR] = x["vE"]
          v2[isSR] = x["vC"]
        }
    out=dLBA(rt=data$rt,response=data$response,A=x["A"],b=bs,t0=x["t0"],mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=TRUE)
      }
    }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
}

acc_loglike_b=function(x,data,sample=FALSE) {
  x=exp(x)

  if (sample) {
    data$rt <- rep(NA,nrow(data))
    data$response <- rep(NA,nrow(data))
  }   else { 
    if (any(data$rt < x["t0"])) {
      return(-1e10)
    }
    out <- numeric(nrow(data))
  }
  
  for (i in 1:nrow(data)) {
    A = x["A"]
    b = x[paste0("b.",data$condition[i])]+A
    vC = x["vC"]
    vE = x["vE"]
    t0 = x["t0"]
    s=c(1,1)
    if (data$stimulus[i]=="left") {
      vs = c(vC,vE)
    } else {
      vs = c(vE,vC)
    }
    
    if (sample) {
      tmp<-rLBA(n=1,A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
      data$rt[i]=tmp$rt
      data$response[i]=c("left","right")[tmp$response]
    }    else {
      out[i]<-dLBA(rt=data$rt[i],response=data$response[i],
                   A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
    }
    
  }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
  
}


###t0
fast_loglike_t0=function(x,data,sample=TRUE) {
  x=exp(x)
  t0s = unname(x[c("t0.easy", "t0.hard")][data$condition])
  bs=x["A"]+x[c("b")]
  if (sample) {
    for (stim in levels(data$stimulus)){
      isSR = data$stimulus==stim
      if (stim == "left"){
        v1 = x["vC"]
        v2 = x["vE"]
      } else {
        v1 = x["vE"]
        v2 = x["vC"]
      }
    out=rLBA(n=sum(isSR),A=x["A"],b=bs,t0=t0s,mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=TRUE)
    data$response[isSR]=c("left","right")[out$response]
    data$rt[isSR]=out$rt
    out<-data
    } 
  }  else {
    v1 = v2 = numeric(length(data$RT))
    for (stim in levels(data$stimulus)){
      isSR = data$stimulus==stim
      if (stim == "left"){
        v1[isSR] = x["vC"]
        v2[isSR] = x["vE"]
      } else {
        v1[isSR] = x["vE"]
        v2[isSR] = x["vC"]
      }
    }
    out=dLBA(rt=data$rt,response=data$response,A=x["A"],b=bs,t0=t0s,mean_v=x[c("v1","v2")],sd_v=c(1,1),distribution="norm",silent=TRUE)
  }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
  
}

acc_loglike_t0=function(x,data,sample=FALSE) {
  x=exp(x)
  if (sample) {
    data$rt <- rep(NA,nrow(data))
    data$response <- rep(NA,nrow(data))
  } 
  else { 
    if (any(data$rt < x["t0"])) {
      return(-1e10)
    }
    out <- numeric(nrow(data))
  }
  for (i in 1:nrow(data)) {
    A = x["A"]
    b = x["b"]+A
    vC = x["vC"]
    vE = x["vE"]
    t0 = x[paste0("t0.",data$condition[i])]
    s=c(1,1)
    if (data$stimulus[i]=="left") {
      vs = c(vC,vE)
    } else {
      vs = c(vE,vC)
    }
    
    if (sample) {
      tmp<-rLBA(n=1,A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
      data$rt[i]=tmp$rt
      data$response[i]=c("left","right")[tmp$response]
    }
    else {
      out[i]<-dLBA(rt=data$rt[i],response=data$response[i],
                   A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
    }
    
  }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
}


####vs
fast_loglike_v=function(x,data,sample=TRUE) {
  x=exp(x)
  if (sample) {  
  for (cond in levels(data$condition)){
    for (stim in levels(data$stimulus)){
      isSR = data$condition==cond & data$stimulus==stim
      if (!any(isSR)) next
      if (cond == "easy") {
        if (stim == "left"){
        v1 = x["vC.easy"]
        v2 = x["vE.easy"]
      } else {
        v1 = x["vE.easy"]
        v2 = x["vC.easy"]
      }
      }
      if (cond == "hard") {
        if (stim == "left"){
          v1 = x["vC.hard"]
          v2 = x["vE.hard"]
        } else {
          v1 = x["vE.hard"]
          v2 = x["vC.hard"]
        }
      }
      
      out=rLBA(n=sum(isSR),A=x["A"],b=x["b"]+x["A"],t0=x["t0"],mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=FALSE)  
      data$response[isSR]=c("left","right")[out$response]
      data$rt[isSR]=out$rt
      out<-data
      }
  }  
}
      else {
        v1 = v2 = numeric(length(data$RT))
        for (cond in levels(data$condition)){
          for (stim in levels(data$stimulus)){
            isSR = data$condition==cond & data$stimulus==stim
            if (!any(isSR)) next
            if (cond == "easy") {
              if (stim == "left"){
                v1[isSR] = x["vC.easy"]
                v2[isSR] = x["vE.easy"]
              } else {
                v1[isSR] = x["vE.easy"]
                v2[isSR] = x["vC.easy"]
              }
            }
            if (cond == "hard") {
              if (stim == "left"){
                v1[isSR] = x["vC.hard"]
                v2[isSR] = x["vE.hard"]
              } else {
                v1[isSR] = x["vE.hard"]
                v2[isSR] = x["vC.hard"]
              }
            }
        
         
        out=dLBA(rt=data$rt,response=data$response,A=x["A"],b=x["b"],t0=x["t0"],mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=TRUE)

          }
        }
      }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
}

acc_loglike_v=function(x,data,sample=FALSE) {
  x=exp(x)
  if (sample) {
    data$rt <- rep(NA,nrow(data))
    data$response <- rep(NA,nrow(data))
  }   else { 
    if (any(data$rt < x["t0"])) {
      return(-1e10)
    }
    out <- numeric(nrow(data))
  }
  
  for (i in 1:nrow(data)) {
    A = x["A"]
    b = x["b"]+A
    vC = x[paste0("vC.",data$condition[i])]
    vE = x[paste0("vE.",data$condition[i])]
    t0 = x["t0"]
    s=c(1,1)
    if (data$stimulus[i]=="left") {
      vs = c(vC,vE)
    } else {
      vs = c(vE,vC)
    }
    
    if (sample) {
      tmp<-rLBA(n=1,A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
      data$rt[i]=tmp$rt
      data$response[i]=c("left","right")[tmp$response]
    }    else {
      out[i]<-dLBA(rt=data$rt[i],response=data$response[i],
                   A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
    }
    
  }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
  
}


fast_loglike_coherence=function(x,data,sample=TRUE) {
  x=exp(x)
  if (sample) {  
    for (cond in levels(data$condition)){
      for (stim in levels(data$stimulus)){
        isSR = data$condition==cond & data$stimulus==stim
        if (!any(isSR)) next
        if (cond == "1") {
          if (stim == "left"){
            v1 = x["vC.1"]
            v2 = x["vE"]
          } else {
            v1 = x["vE"]
            v2 = x["vC.1"]
          }
        }
       else if (cond == "2") {
          if (stim == "left"){
            v1 = x["vC.2"]
            v2 = x["vE"]
          } else {
            v1 = x["vE"]
            v2 = x["vC.2"]
          }
       }
        else if (cond == "3") {
          if (stim == "left"){
            v1 = x["vC.3"]
            v2 = x["vE"]
          } else {
            v1 = x["vE"]
            v2 = x["vC.3"]
          }
        }
        else if (cond == "4") {
          if (stim == "left"){
            v1 = x["vC.4"]
            v2 = x["vE"]
          } else {
            v1 = x["vE"]
            v2 = x["vC.4"]
          }
        }
        
        out=rLBA(n=sum(isSR),A=x["A"],b=x["b"]+x["A"],t0=x["t0"],mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=FALSE)  
        data$response[isSR]=c("left","right")[out$response]
        data$rt[isSR]=out$rt
        out<-data
      }
    }  
  }
  else {
    v1 = v2 = numeric(length(data$rt))
    for (cond in levels(data$condition)){
      for (stim in levels(data$stimulus)){
        isSR = data$condition==cond & data$stimulus==stim
        if (cond == "1") {
          if (stim == "left"){
            v1[isSR] = x["vC.1"]
            v2[isSR] = x["vE"]
          } else {
            v1[isSR] = x["vE"]
            v2[isSR] = x["vC.1"]
          }
        }
        else if (cond == "2") {
          if (stim == "left"){
            v1[isSR] = x["vC.2"]
            v2[isSR] = x["vE"]
          } else {
            v1[isSR] = x["vE"]
            v2[isSR] = x["vC.2"]
          }
        }
        else if (cond == "3") {
          if (stim == "left"){
            v1[isSR] = x["vC.3"]
            v2[isSR] = x["vE"]
          } else {
            v1[isSR] = x["vE"]
            v2[isSR] = x["vC.3"]
          }
        }
        else if (cond == "4") {
          if (stim == "left"){
            v1[isSR] = x["vC.4"]
            v2[isSR] = x["vE"]
          } else {
            v1[isSR]= x["vE"]
            v2[isSR] = x["vC.4"]
          }
        }
        out[isSR]=dLBA(rt=data$rt[isSR],response=data$response[isSR],A=x["A"],b=x["b"],t0=x["t0"],mean_v=list(v1,v2),sd_v=c(1,1),distribution="norm",silent=TRUE)
      }
    }
  }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
}


acc_loglike_coherence=function(x,data,sample=FALSE) {
  x=exp(x)
  if (sample) {
    data$rt <- rep(NA,nrow(data))
    data$response <- rep(NA,nrow(data))
  }   else { 
    # if (any(data$rt < x["t0"])) {
    #   return(-1e10)
    # }
    out <- numeric(nrow(data))
  }
  
  for (i in 1:nrow(data)) {
    A = x["A"]
    b = x["b"]+A
    vC = x[paste0("vC.",data$condition[i])]
    vE = x["vE"]
    t0 = x["t0"]
    s=c(1,1)
    if (data$stimulus[i]=="left") {
      vs = c(vC,vE)
    } else {
      vs = c(vE,vC)
    }
    
    if (sample) {
      tmp<-rLBA(n=1,A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
      data$rt[i]=tmp$rt
      data$response[i]=c("left","right")[tmp$response]
    }    else {
      out[i]<-dLBA(rt=data$rt[i],response=data$response[i],
                   A=A,b=b,mean_v=vs,sd_v=s,t0=t0,dist="norm",silent=TRUE)
    }
    
  }
  if (sample) return(data)
  if (!sample) return(sum(log(pmax(out,1e-10))))
  
}

### Here is where i make the covariance matrix 

vars = abs(ptm)/10
sigmaC = matrix(.2,nrow=length(ptm),ncol=length(ptm)) ###correlation 
diag(sigmaC)=sqrt(vars)
sigma_matrix <- sdcor2cov(sigmaC)
subj_random_effects <- t(rmvnorm(n.subj,mean=ptm,sigma=sigma_matrix))

data_gen_mu <- ptm
data_gen_sigma <- sigma_matrix
data_gen_subj_re <- subj_random_effects


tmp=split(x=data,f=data$subject)

#generate data
if (drift ==TRUE){
  for (i in 1:n.subj){
    out <- acc_loglike_v(subj_random_effects[,i],sample=TRUE,do.call(rbind.data.frame, tmp[i]))
    data$rt[data$subject==i]=out$rt
    data$response[data$subject==i]=out$response
  }
  setwd("~/Documents/Research/Modelling Project/Scripts/recovery/xFit")
  write.csv(data, "xFit_vs.csv")
  write.csv(data_gen_subj_re, "xFit_vs_alpha.csv")
  write.csv(data_gen_mu, "xFit_vs_mu.csv")
  write.csv(data_gen_sigma, "xFit_vs_sig.csv")
}

if (threshold == TRUE){
  for (i in 1:n.subj){
    out<- acc_loglike_b(subj_random_effects[,i],sample=TRUE,do.call(rbind.data.frame, tmp[i]))
    data$rt[data$subject==i]=out$rt
    data$response[data$subject==i]=out$response
  }
  setwd("~/Documents/Research/Modelling Project/Scripts/recovery/xFit")
  # write.csv(data, "xFit_bs.csv")
  # write.csv(data_gen_subj_re, "xFit_bs_alpha.csv")
  # write.csv(data_gen_mu, "xFit_bs_mu.csv")
  # write.csv(data_gen_sigma, "xFit_bs_sig.csv")
}


if (non.decision== TRUE){
  for (i in 1:n.subj){
    out<- acc_loglike_t0(subj_random_effects[,i],sample=TRUE,do.call(rbind.data.frame, tmp[i]))
    data$rt[data$subject==i]=out$rt
    data$response[data$subject==i]=out$response
  }
  setwd("~/Documents/Research/Modelling Project/Scripts/recovery/xFit")
  write.csv(data, "xFit_t0s.csv")
  write.csv(data_gen_subj_re, "xFit_t0s_alpha.csv")
  write.csv(data_gen_mu, "xFit_t0s_mu.csv")
  write.csv(data_gen_sigma, "xFit_t0s_sig.csv")
}

if (coherence ==TRUE){
  for (i in 1:n.subj){
    out <- acc_loglike_coherence(subj_random_effects[,i],sample=TRUE,do.call(rbind.data.frame, tmp[i]))
    data$rt[data$subject==i]=out$rt
    data$response[data$subject==i]=out$response
  }
  setwd("~/Documents/Research/Modelling Project/Scripts/recovery/xFit")
  write.csv(data, "xFit_big.csv")
  write.csv(data_gen_subj_re, "xFit_big_alpha.csv")
  write.csv(data_gen_mu, "xFit_big_mu.csv")
  write.csv(data_gen_sigma, "xFit_big_sig.csv")
}

