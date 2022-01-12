library("MASS")
library("numDeriv")
library("Rcpp")
library("RcppEigen")
sourceCpp("~/secondary_phenotype/POE_git/function.cpp")

# -------------- Analysis of parent-of-origin effects for secondary phenotype ----------- #
# author: Chen Shuyue, USTC                                                               #
# description: functions for analysis of parent-of-origin effects for secondary phenotype #
#              using case-control mother-child pair data                                  #
# Model:                                                                                  #
# logit P(Y=1|gmm,gcc,X) = beta0 + beta1*gmm + beta2*gcc + beta3*(gcm-gcp) + beta4*x      #       
# logit P(D=1|gmm,gcc,y,x) = delta0 + delta1*gmm + delta2*gcc + delta3*y + delta4*x       #
#                                                                                         #
# parameters: Theta=(beta0,beta1,beta2,beta3,beta4,delta0,delta1,delta2,delta3,delta4,pa) #     
#                                                                                         #
# function:                                                                               #
# -HAP_POE(y, d, gmm, gcc, x, m, hap, ppi, f, N0, N1): obtain the point estimates,        #
#   standard error, Wald test statistics, and significance test p-values for parameters.  #
#                                                                                         #
# input:                                                                                  #
# -y: a n-vector of binary secondary phenotype for n children                             #
# -d: a n-vector of disease statuses for n children (1 for case and 0 for control)        #
# -gmm: a n x k matrix of genotypes for mothers (n: number of mothers; k: number of SNPs).#
#       The possible values should be 0, 1, 2.                                            #
# -gcc: a n x k matrix of genotypes for children (n: number of children; k: number of     #
#       SNPs). The possible values should be 0, 1, 2.                                     #
# -x: a n-vector of maternal covariate for n mothers                                      #
# -m: an indicator for the test locus. The possible values could be 1,...,k.              # 
#     (k: number of SNPs).                                                                #
# -hap: a l x k matrix of possible haplotypes in the population of interest.              # 
#       (l: number of possible haplotypes; k: number of SNPs).                            #
# -ppi: a l-vector of the corresponding haplotype frequencies (l: number of possible      #
#       haplotypes).                                                                      #
# -f: specified disease prevalence                                                        #
# -N0: number of cases in the population                                                  #
# -N1: number of controls in the population                                               #
#                                                                                         #
# output:                                                                                 #
# -hap: estimation and significance test results for the our proposed method HAP          #
# -sin: estimation and significance test results for the our proposed method SIN          #
# -est: point estimates                                                                   #
# -se: standard error                                                                     #
# -z.score: Wald test statistics                                                          #
# -pval: significance test p-values                                                       #
# -cov.hap: covariance matrix of the estimated parameters by method HAP                   # 
# -cov.sin: covariance matrix of the estimated parameters by method SIN                   #                                                        
# --------------------------------------------------------------------------------------- #

HAP_POE <- function(y, d, gmm, gcc, x, m, hap, ppi, f, N0, N1) {
  n.snp <- ncol(gmm)
  hap <- as.matrix(hap)
  gmm <- as.matrix(gmm)
  gcc <- as.matrix(gcc)
  
  # subset of complete data
  data <- cbind(y, d, gmm, gcc, x)
  data_nomiss <- data[complete.cases(data), ]
  y <- data_nomiss[, 1]
  d <- data_nomiss[, 2]
  gmm <- data_nomiss[, 3:(n.snp+2)]
  gcc <- data_nomiss[, (n.snp+3):(2*n.snp+2)]
  x <- data_nomiss[, 2*n.snp+3]
  
  # compatibility with haplotypes
  comp <- distinguish_data(gmm, gcc, hap)
  if (sum(which(is.na(comp)))>0) {
    gmm <- gmm[-which(is.na(comp)), ]
    gcc <- gcc[-which(is.na(comp)), ]
    x <- x[-which(is.na(comp))]
    y <- y[-which(is.na(comp))]
    d <- d[-which(is.na(comp))]
  } else {
    gmm <- gmm
    gcc <- gcc
    x <- x
    y <- y
    d <- d
  }
  
  # inference of parental origins through haplotypes
  F <- distinguishC(gmm, gcc, m, hap)
  gm <- F$gm
  gc <- F$gc
  gcm <- F$gcm
  gcp <- F$gcp
  phi <- F$phi
  
  # get initial values
  pa.ini <- iwp_mle_pa(d, gm, gc, N0, N1)
  
  Data1 <- cbind(d, gm, gc, y, x)
  Data1 <- as.data.frame(Data1)
  result1 <- glm(d~gm+gc+y+x, family = binomial(link = "logit"), data = Data1)
  delta.ini <- c(-2, result1$coefficients[-1])
  
  Data2 <-cbind(y, gm, gc, x)
  Data2 <- as.data.frame(Data2)
  result2 <-glm(y~gm+gc+x,family = binomial(link = "logit"), data = Data2)
  beta.ini <- c(-2, result2$coefficients[c(2,3)], 0, result2$coefficients[4])
  
  Theta0 <- c(beta.ini, delta.ini, pa.ini)
  n.para <- length(Theta0)
  
  # upper and lower bound of parameters
  upper = lower = Theta0
  upper[n.para] = 0.99
  lower[n.para] = 0.01
  upper[-n.para] = Theta0[-n.para] + 10
  lower[-n.para] = Theta0[-n.para] - 10
  
  # maximum likelihood estimator of Theta by HAP
  eval_f <- function(parameter) {
    beta.h <- parameter[1:5]
    delta.h <- parameter[6:10]
    pa.h <- parameter[11]
    be.optim <- log_lmp_hapC(y = y, d = d, gmm = gmm, gcc = gcc, x = x,
                             beta = beta.h, delta = delta.h, 
                             m = m, hap = hap, ppi = ppi, pa = pa.h, f = f)
    return(-be.optim)
  }
  eval_grad_f <- function(parameter) {
    beta.h <- parameter[1:5]
    delta.h <- parameter[6:10]
    pa.h <- parameter[11]
    be.optim <- log_lmp_hap_grC(y = y, d = d, gmm = gmm, gcc = gcc, x = x,
                                beta = beta.h, delta = delta.h, 
                                m = m, hap = hap, ppi = ppi, pa = pa.h, f = f)
    return(-be.optim)
  }
  fit <- optim(par = Theta0, fn = eval_f, gr = eval_grad_f, method = "L-BFGS-B", 
               lower = lower, upper = upper, hessian = TRUE)
  
  # estimation results
  est.hap <- fit$par
  V <- fit$hessian
  U <- cov_log_lmp_hap_grC(y, d, gmm, gcc, x, beta = est[1:5], delta = est[6:10], m, hap, ppi, pa = est[11], f)
  V.inv <- ginv(V)
  cov.hap <- (V.inv%*%U)%*%V.inv
  se.hap <- sqrt(diag(cov.hap))
  z.hap <- est.hap/se.hap
  pval.hap <- 1 - pchisq(z.hap^2, 1)
  
  # maximum likelihood estimator of Theta by SIN
  eval_f <- function(parameter) {
    beta.h <- parameter[1:5]
    delta.h <- parameter[6:10]
    pa.h <- parameter[11]
    be.optim <- log_lmpC(y = y, d = d, gm = gm, gc = gc, x = x, beta = beta.h, delta = delta.h, pa = pa.h, f = f)
    return(-be.optim)
  }
  eval_grad_f <- function(parameter) {
    beta.h <- parameter[1:5]
    delta.h <- parameter[6:10]
    pa.h <- parameter[11]
    be.optim <- log_lmp_grC(y = y, d = d, gm = gm, gc = gc, x = x, beta = beta.h, delta = delta.h, pa = pa.h, f = f)
    return(-be.optim)
  }
  fit <- optim(par = Theta0, fn = eval_f, gr = eval_grad_f, method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE)
  est.sin <- fit$par
  V <- fit$hessian
  U <- cov_log_lmp_grC(y, d, gm, gc, x, beta = est[1:5], delta = est[6:10], pa = est[11], f)
  V.inv <- ginv(V)
  cov.sin <- (V.inv%*%U)%*%V.inv
  se.sin <- sqrt(diag(cov.sin))
  z.sin <- est.sin/se.sin
  pval.sin <- 1 - pchisq(z.sin^2, 1)
  
  return(list(hap = data.frame(est=est.hap,se=se.hap,z.score=z.hap,pval=pval.hap),
              sin = data.frame(est=est.sin,se=se.sin,z.score=z.sin,pval=pval.sin),
              cov.hap = cov.hap, cov.sin = cov.sin))
}

# MLE estimate of pa, using inverse-weighted-probability             
iwp_mle_pa <- function(d, gm, gc, N0, N1) {
  n.mc <- table(factor(gm, 0:2)[d == 0], factor(gc, 0:2)[d == 0]) * N0 / (N0 + N1) + table(factor(gm, 0:2)[d == 1], factor(gc, 0:2)[d == 1]) * N1 / (N0 + N1)
  pa.mle <- (n.mc[1, 2] + n.mc[2, 1] + n.mc[2, 2] + 2 * n.mc[2, 3] + 2 * n.mc[3, 2] + 3 * n.mc[3, 3]) / (3 * sum(n.mc) - n.mc[2, 2])
  return(pa.mle)
}

# distinguish whether the data fit the haplotype
distinguish_data = function(gmm,gcc,hap){
  l = max(nrow(gmm),1)
  K = nrow(hap)
  comp = rep(NA,l)
  
  for(u in 1:l){
    
    for(k in 1:K){
      #mother
      hm1 = hap[k,]
      hm2 = gmm[u,] - hm1
      r = 0
      for (kk in 1:K){
        r = r + all(hm2==hap[kk,])
      }
      if(r==1){
        #pass hm1 to child
        hcm = hm1
        hcp = gcc[u,] - hcm
        t = 0
        for (kk in 1:K){
          t = t + all(hcp==hap[kk,])
        } 
        if(t==1) comp[u] = 1
      } 
    }
  }
  comp 
}



