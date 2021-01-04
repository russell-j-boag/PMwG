#Functions to use with PMwG#

### from the sampled object, uses the LL (sample =TRUE) to create data. Relies on LL being correct. Can have errors with lists
generate.posterior <- function(sampled, n){
  n.posterior=n # Number of parameter samples from posterior distribution.
  pp.data=list()
  S = sampled$n_subjects
  data=sampled$data
  sampled_stage = length(sampled$samples$stage[sampled$samples$stage=="sample"])
  for (s in 1:S) {
    cat(s," ")
    iterations=round(seq(from=(sampled$samples$idx-sampled_stage) , to=sampled$samples$idx, length.out=n.posterior))
    for (i in 1:length(iterations)) {
      x <- sampled$samples$alpha[,s,iterations[i]]
      names(x) <- sampled$par_names
      tmp=sampled$ll_func(x=x,data=data[as.integer(as.numeric(data$subject))==s,],sample=TRUE)
      if (i==1) {
        pp.data[[s]]=cbind(i,tmp)
      } else {
        pp.data[[s]]=rbind(pp.data[[s]],cbind(i,tmp))
      }
    }
  }
  return(pp.data)
}
#tmp<-generate.posterior(sampled,10)
#tmp=do.call(rbind,tmp)





### DIC for a PMwG object
pmwg.DIC=function(sampled,pD=FALSE){
  nsubj=length(unique(sampled$data$subject))
  
  # the mean likelihood of the overall (sampled-stage) model, separately for each subject
  mean.like <- apply(sampled$samples$subj_ll[,sampled$samples$stage=="sample"],1,mean)
  
  # the mean of each parameter across iterations. Keep dimensions for parameters and subjects
  mean.params <- t(apply(sampled$samples$alpha,1:2,mean))
  
  # i name mean.params here so it can be used by the log_like function
  colnames(mean.params)<-sampled$par_names
  
  # log-likelihood for each subject using their mean parameter vector
  mean.params.like <- numeric(ncol(mean.params))
  data <- transform(sampled$data, subject=match(subject, unique(subject)))
  for (j in 1:nsubj) {
    mean.params.like[j] <- sampled$ll_func(mean.params[j,], data=data[data$subject==j,], sample=FALSE)
  }
  
  # Effective number of parameters
  pD <- sum(-2*mean.like + 2*mean.params.like)
  
  # Deviance Information Criterion
  DIC <- sum(-4*mean.like + 2*mean.params.like)
  
  if (pD){
    return(c("DIC"=DIC,"effective parameters"=pD))
  }else{
    return(DIC)
  }
    
}



### Checks used in sampler Doc
tmp <- sampled
dev.off()
#### Chains should be separated (not centred around 0) with low(-ish) variance. If one thick band, concern. 
matplot(t(tmp$samples$theta_mu),type="l")
#### Chains should be separated and static (not trending up or down) with some noise. If only one line is present, concern. 
matplot(t(tmp$samples$subj_ll),type="l")

#### checks the number of new particles per subject (useful for checking when it doesn't get out of adaptation stage)
x<-apply(tmp$samples$alpha[1,,-1]!=tmp$samples$alpha[1,,-(tmp$samples$idx)],1,sum)
x[order(x)]

#### check theta parameter values (this uses exp as we exp in the LL)
tmp <- exp(sampled$samples$alpha[,,sampled$samples$idx])
round(tmp,3)


#### Covariance matrix
cov<-apply(sampled$samples$theta_sig[,,sampled$samples$idx-1000:sampled$samples$idx] ,1:2, mean)
colnames(cov)<-pars
rownames(cov)<-pars
cor<-cov2cor(cov) #correlation matrix


#### checking covariance
diagonal<-apply(tmp$samples$theta_sig,3,diag)
matplot(log(t(diagonal)), type="l")

### function to get the IACT value for each parameter and the covariance off diagonal
library(LaplacesDemon)
pmwg.iact = function(sampled){
  dim = sampled$n_pars
  n.params<- (dim*dim - dim)/2 
  n_subjects = sampled$n_subjects
  theta <- sampled$samples$theta_mu[,sampled$samples$stage=="sample"]
  alpha <- sampled$samples$alpha[,,sampled$samples$stage=="sample"]
  sig <- sampled$samples$theta_sig[,,sampled$samples$stage=="sample"]  
  tmp <- sig[,,1:100]
  
  tmp<-NULL
  for (i in 1:dim){
    tmp[i]<- IAT(theta[i,])
  }
  names(tmp)<-pars
  
  tmp2<-matrix(nrow=dim, ncol = dim)
  for (i in 1:dim){
    for (j in 1:dim){
      tmp2[i,j]<- IAT(sig[i,j,])
    }
  }
  
  # for (i in 1:n_subjects){
  #   for (j in 1:dim){
  #     tmp2[]
  #   }
  # }
  # 
  
  colnames(tmp2)<-pars
  rownames(tmp2)<-paste0(pars,"_")
  tmp3 <- tmp2[lower.tri(tmp2)]
  index <-  which(lower.tri(tmp2), arr.ind = T)
  names(tmp3) <- paste(rownames(tmp2)[index[,1]], colnames(tmp2)[index[,2]], sep=", ")
  
  iact.table <- c(tmp,tmp3)
  
  return(iact.table)
} 

