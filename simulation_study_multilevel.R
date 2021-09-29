#########################################################################################
#########################################################################################
##
## Simulation study for evaluating the ordered probit model with sample selection
## and random intercept
##
## Consideration of different missing-data mechanisms: 
## -> MAR, MNAR selection model, MNAR non-selection model
##
## @author Angelina Hammon
##
## 29.09.2021
##
#########################################################################################
#########################################################################################


# install package for model calculation: 
install.packages("bivprob.quad_2.0.tar.gz", repos = NULL, type="source")

library(mvtnorm)
library(mice)
library(doParallel)
library(MASS)
library(lme4)
library(PanelCount)
library(doSNOW)
library(Rcpp)
library(ordinal)
library(miceadds)
library(bivprob.quad)


### true values ### 

model_bd <- vector("numeric",length=3)

rhos <- c(0,0.3,0.6,0.9)
taus <- c(0,0.1,0.3,0.5)

# Non-Heckman MNAR & MAR
rho <- rhos[1]
tau <- taus[1]

# Heckman MNAR:
#rho <- rhos[4]  # change depending on scenario
#tau <- taus[4]

for(i in 1:100000) {
  n <- 2500
  m <- 20
  nj <- n/m
  
  x1 <- rnorm(n,0,0.3)
  x2 <- rnorm(n,0,0.8)
  x3 <- rnorm(n,0,4)
  
  vc <- diag(2)
  vc[2,1] <- vc[1,2] <- rho
  eps <- rmvnorm(n,rep(0,2),vc)
  eps_o <- eps[,1]
  eps_s <- eps[,2]
  
  vc_re <- matrix(1,2,2)
  vc_re[1,1] <- 0.5
  vc_re[2,2] <- 0.9
  vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
  
  alpha <- rmvnorm(m,rep(0,2),vc_re) 
  alpha_s <- rep(alpha[,1], each=nj)  
  alpha_o <- rep(alpha[,2], each=nj)  
  
  y_star <- 1*x1+0.5*x2+eps_o+alpha_o
  z <- c(-0.75,0.5)
  z1 <- z[1]
  z2 <- z[2]
  y <- as.factor(ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
  
  r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
  r <- ifelse(r<0,0,1)
  
  group <- rep(1:m, each=nj) 
  
  data <- data.frame(y,r,x1,x2,x3,group)
  
  model <- clmm(y~x1+x2+(1|group),data=data, link="probit",nAGQ=5)

  model_bd <- model_bd + c(model$beta,VarCorr(model)$'group'[1])
  
}

model_bd <- model_bd/100000



### Functions ###

coverage <- function(value, CI.low, CI.upper){
  ifelse(CI.low <= value && CI.upper >= value,1,0)
}

rel.bias <- function(value,est) {
  rel_bias <- (est-value)/value
  return(rel_bias)
}

RMSE <- function(value,est){
  mse <- (est-value)^2
  rmse <- sqrt(mse)
  return(mse)
}


### Data Generation ###

cl <- makeSOCKcluster(8)
registerDoSNOW(cl)
clusterExport(cl,"rel.bias")
clusterExport(cl,"coverage")
clusterExport(cl,"RMSE")
clusterExport(cl,"model_bd")
clusterExport(cl,"rhos")
clusterExport(cl,"taus")
clusterEvalQ(cl, {
  library(mice)
  library(mvtnorm)
  library(doParallel)
  library(lme4)
  library(PanelCount)
  library(ordinal)
  library(MASS)
  library(Hmisc)
  library(Rcpp)
  library(RcppArmadillo)
  library(corpcor)
  library(compiler)
  library(nloptr)
  library(bivprob.quad)
  library(miceadds)
  source("function_MI_analysis.R")
  source("mice_heckman_ordered_re1.R")
  source("mice_heckman_ordered_re2.R")
})
    
    
start_time <- Sys.time()
sim <- 
  parLapply(cl = cl, 1:500, fun = function(no, mech="non-heckman",rho=rhos[1],tau=taus[1],M=5){
    n <- 2500
    
    m <- 20
    nj <- n/m
    
    x1 <- rnorm(n,0,0.3)
    x2 <- rnorm(n,0,0.8)
    x3 <- rnorm(n,0,4)
    
    
    ## Missing mechanisms ##
    
    if(mech=="heckman") {
      
      if(rho == 0) {stop("rho is equal to 0")}
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      vc_re <- matrix(1,2,2)
      vc_re[1,1] <- 0.5
      vc_re[2,2] <- 0.9
      vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
      
      alpha <- rmvnorm(m,rep(0,2),vc_re) 
      alpha_s <- rep(alpha[,1], each=nj)  
      alpha_o <- rep(alpha[,2], each=nj)  
      
      y_star <- 1*x1+0.5*x2+eps_o+alpha_o
      z <- c(-0.75,0.5)
      z1 <- z[1]
      z2 <- z[2]
      y <- as.factor(ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
      
      r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
      r <- ifelse(r<0,0,1)
      
      
    } else if(mech=="non-heckman"){
      
      if(rho != 0) {stop("rho is unequal to 0")}
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      vc_re <- matrix(1,2,2)
      vc_re[1,1] <- 0.5
      vc_re[2,2] <- 0.9
      vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
      
      alpha <- rmvnorm(m,rep(0,2),vc_re) 
      alpha_s <- rep(alpha[,1], each=nj)  
      alpha_o <- rep(alpha[,2], each=nj)  
      
      y_star <- 1*x1+0.5*x2+eps_o+alpha_o
      z <- c(-0.75,0.5)
      z1 <- z[1]
      z2 <- z[2]
      y <- as.factor(ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
      
      # nicht so "harter" Mechanismus hinsichtlich y..:
      prob <- pnorm(1.25+1.75*y_star+1.5*x1-2.5*x2+alpha_s)
      r <- rbinom(n,1,p=prob)
      
    } else{ # MAR
      if(rho != 0) {stop("rho is unequal to 0")}
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      vc_re <- matrix(1,2,2)
      vc_re[1,1] <- 0.5
      vc_re[2,2] <- 0.9
      vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
      
      alpha <- rmvnorm(m,rep(0,2),vc_re) 
      alpha_s <- rep(alpha[,1], each=nj)  
      alpha_o <- rep(alpha[,2], each=nj)  
      
      y_star <- 1*x1+0.5*x2+eps_o+alpha_o
      z <- c(-0.75,0.5)
      z1 <- z[1]
      z2 <- z[2]
      y <- as.factor(ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
      
      r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
      r <- ifelse(r<0,0,1)
      
    }
    
    
    group <- rep(1:m, each=nj) 
    
    
    ## BD ##
    
    data_comp <- data.frame(y,r,x1,x2,x3,group)
    
    bd <- clmm(y~x1+x2+(1|group),data=data_comp, link="probit",nAGQ=5)
    
    CI_bd <- confint(bd,method="Wald")[3:(length(bd$beta)+2),]
    
    cov_bd <- NULL
    for(j in 1:(length(model_bd)-1)){
      cov_bd <- c(cov_bd,coverage(model_bd[j],CI_bd[j,1],CI_bd[j,2]))
    }
    
    rmse.est_bd <- RMSE(model_bd,c(bd$beta,bd$ST$group^2))
    
    var_bd <- c((coef(summary(bd))[3:(length(bd$beta)+2),2])^2,0)
    
    
    ## Generation of missing values ##
    
    # based on upper specified mechanism in y:
    y[r==0] <- NA
    
    data <- data.frame(y,r,x1,x2,x3,group)
    
    
    ### Imputation ###
    
    ini <- mice(data,m=1,maxit=0)
    
    ## MNAR Heckman AGHQ ##
    pred_MNAR <- ini$pred
    pred_MNAR["y","r"] <- 0
    pred_MNAR["y","group"] <- -2
    
    imp_MNAR_aghq <- mice(data,m=M,maxit=1,method=c("2l.heckman1step_ord_re_aghq","","","","",""),pred=pred_MNAR,print=F,QP=rep(10,2),draw=T,excl="x3",seed=1234)
    
    glm_MNAR_aghq <- with(imp_MNAR_aghq, clmm(y~x1+x2+(1|group),link="probit",nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MNAR_aghq$analyses[[x]]$beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MNAR_aghq$analyses[[x]]))[3:(length(bd$beta)+2)]))
    re <- t(sapply(1:M,function(x)  glm_MNAR_aghq$analyses[[x]]$ST$group^2))
    pool_glm_MNAR_aghq <- MI.analysis(est,var,m=M)
    pool_glm_MNAR_aghq <- rbind(pool_glm_MNAR_aghq,c(mean(re),rep(0,3)))
    
    cov_MNAR_aghq <- NULL
    for(j in 1:ncol(est)){
      cov_MNAR_aghq <- c(cov_MNAR_aghq,coverage(model_bd[j],pool_glm_MNAR_aghq[j,2],pool_glm_MNAR_aghq[j,3]))
    }
    
    rmse.est_MNAR_aghq <- RMSE(model_bd,pool_glm_MNAR_aghq[,1])
    
    var_MNAR_aghq <- pool_glm_MNAR_aghq[,4]
    
    
    ## MNAR Heckman GHQ ##
    imp_MNAR_ghq <- mice(data,m=M,maxit=1,method=c("2l.heckman1step_ord_re_ghq","","","","",""),pred=pred_MNAR,print=F,QP=rep(10,2),draw=T,excl="x3")
    
    glm_MNAR_ghq <- with(imp_MNAR_ghq, clmm(y~x1+x2+(1|group),link="probit",nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MNAR_ghq$analyses[[x]]$beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MNAR_ghq$analyses[[x]]))[3:(length(bd$beta)+2)]))
    re <- t(sapply(1:M,function(x)  glm_MNAR_ghq$analyses[[x]]$ST$group^2))
    pool_glm_MNAR_ghq <- MI.analysis(est,var,m=M)
    pool_glm_MNAR_ghq <- rbind(pool_glm_MNAR_ghq,c(mean(re),rep(0,3)))
    
    cov_MNAR_ghq <- NULL
    for(j in 1:ncol(est)){
      cov_MNAR_ghq <- c(cov_MNAR_ghq,coverage(model_bd[j],pool_glm_MNAR_ghq[j,2],pool_glm_MNAR_ghq[j,3]))
    }
    
    rmse.est_MNAR_ghq <- RMSE(model_bd,pool_glm_MNAR_ghq[,1])
    
    var_MNAR_ghq <- pool_glm_MNAR_ghq[,4]
    
    
    ## MAR - PMM ##     
    data2 <- data
    data2$y <- as.numeric(as.character(data2$y))
    
    ini <- mice(data2,m=1,maxit=0)
    
    pred_MAR <- ini$pred
    pred_MAR[,"r"] <- 0
    pred_MAR["y","group"] <- -2
    
    imp_MAR <- mice(data2,m=M,maxit=1,method=c("2l.pmm","","","","",""),pred=pred_MAR,print=F)
    
    glm_MAR <- with(imp_MAR, clmm(as.factor(y)~x1+x2+(1|group),link="probit",nAGQ=5)) 
    est <- t(sapply(1:M,function(x) glm_MAR$analyses[[x]]$beta))
    var <- t(sapply(1:M,function(x) diag(vcov(glm_MAR$analyses[[x]]))[3:(length(bd$beta)+2)]))
    re <- t(sapply(1:M,function(x)  glm_MAR$analyses[[x]]$ST$group^2))
    pool_glm_MAR <- MI.analysis(est,var,m=M)
    pool_glm_MAR <- rbind(pool_glm_MAR,c(mean(re),rep(0,3)))
    
    cov_MAR <- NULL
    for(j in 1:ncol(est)){
      cov_MAR <- c(cov_MAR,coverage(model_bd[j],pool_glm_MAR[j,2],pool_glm_MAR[j,3]))
    }
    
    rmse.est_MAR <- RMSE(model_bd,pool_glm_MAR[,1])
    
    var_MAR <- pool_glm_MAR[,4]
    
    
    ## CC ##
    cc <- clmm(y~x1+x2+(1|group),data=data, link="probit",nAGQ=5)
    
    CI_cc <- confint(cc,method="Wald")[3:(length(cc$beta)+2),]
    
    cov_CC <- NULL
    for(j in 1:(length(model_bd)-1)){
      cov_CC <- c(cov_CC,coverage(model_bd[j],CI_cc[j,1],CI_cc[j,2]))
    }
    
    rmse.est_cc <- RMSE(model_bd,c(cc$beta,cc$ST$group^2)) 
    
    var_cc <- c((coef(summary(cc))[3:(length(cc$beta)+2),2])^2,0)
    
    
    results <- cbind(pool_glm_MAR[,1],
                     pool_glm_MNAR_aghq[,1],
                     pool_glm_MNAR_ghq[,1],
                     c(cc$beta,cc$ST$group^2), c(bd$beta,bd$ST$group^2),
                     rel.bias(model_bd,pool_glm_MAR[,1]),
                     rel.bias(model_bd,pool_glm_MNAR_aghq[,1]), 
                     rel.bias(model_bd,pool_glm_MNAR_ghq[,1]),  
                     rel.bias(model_bd, c(cc$beta,cc$ST$group^2)), rel.bias(model_bd,c(bd$beta,bd$ST$group^2)),
                     cov_MAR, cov_MNAR_aghq, cov_MNAR_ghq,
                     cov_CC, cov_bd,
                     rmse.est_MAR, rmse.est_MNAR_aghq,rmse.est_MNAR_ghq,
                     rmse.est_cc,rmse.est_bd,
                     var_MAR, var_MNAR_aghq,var_MNAR_ghq,
                     var_cc,var_bd)
    
    results[3,c(11:15)] <- 0  
    
    rm(data,data_comp,imp_MNAR_aghq,imp_MNAR_ghq,imp_MAR)   
    
    return(results)

  })
end_time <- Sys.time() 
end_time - start_time 


## Diagnostics ##

est_MAR <- rowMeans(sapply(sim,function(x) x[,1]))
est_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,2]))
est_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,3]))
est_CC <- rowMeans(sapply(sim,function(x) x[,4]))
est_bd <- rowMeans(sapply(sim,function(x) x[,5]))

bias_MAR <- rowMeans(sapply(sim,function(x) x[,6]))
bias_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,7]))
bias_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,8]))
bias_CC <- rowMeans(sapply(sim,function(x) x[,9]))
bias_bd <- rowMeans(sapply(sim,function(x) x[,10]))

cov_MAR <- rowMeans(sapply(sim,function(x) x[,11]))
cov_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,12]))
cov_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,13]))
cov_CC <- rowMeans(sapply(sim,function(x) x[,14]))
cov_bd <- rowMeans(sapply(sim,function(x) x[,15]))

mse_MAR <- rowMeans(sapply(sim,function(x) x[,16]))
mse_MNAR_aghq <- rowMeans(sapply(sim,function(x) x[,17]))
mse_MNAR_ghq <- rowMeans(sapply(sim,function(x) x[,18]))
mse_CC <- rowMeans(sapply(sim,function(x) x[,19]))
mse_bd <- rowMeans(sapply(sim,function(x) x[,20]))  

res <- cbind(est_MAR, est_MNAR_aghq,  est_MNAR_ghq, est_CC, est_bd, 
             bias_MAR, bias_MNAR_aghq, bias_MNAR_ghq, bias_CC, bias_bd, 
             cov_MAR, cov_MNAR_aghq, cov_MNAR_ghq, cov_CC, cov_bd, 
             mse_MAR, mse_MNAR_aghq, mse_MNAR_ghq, mse_CC, mse_bd)
