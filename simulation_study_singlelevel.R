#########################################################################################
#########################################################################################
##
## Simulation study for evaluating the ordered probit model with sample selection
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

library(mvtnorm)
library(mice)
library(doParallel)
library(MASS)


### true values ### 

model_bd <- vector("numeric",length=4)

# change rho depending on scenario

# Non-Heckman MNAR & MAR
rho <- 0

# Heckman MNAR:
# weak:
#rho <- 0.3
# medium:
#rho <- 0.6
# strong: 
#rho <- 0.9

for(i in 1:100000) {

  n <- 2000
  
  x1 <- rnorm(n,0,0.3)
  x2 <- rnorm(n,0,0.8)
  x3 <- rnorm(n,0,4)
  
  vc <- diag(2)
  vc[2,1] <- vc[1,2] <- rho
  eps <- rmvnorm(n,rep(0,2),vc)
  eps_o <- eps[,1]
  eps_s <- eps[,2]
  
  r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s
  r <- as.factor(ifelse(r<0,0,1))
  
  y_star <- 1*x1+0.5*x2+eps_o
  z <- c(-0.75,0.5)
  
  z1 <- z[1]
  z2 <- z[2]
  
  y <- as.factor( ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
  
  data <- data.frame(y,r,x1,x2,x3)
  
  model <- polr(y~x1+x2, data=data,method="probit")
  
  model_bd <- model_bd + c(coef(model), model$zeta)
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

MSE <- function(value,est){
  mse <- (est-value)^2
  return(mse)
}


cl <- makeCluster(15)
clusterSetRNGStream(cl)
clusterExport(cl,"rel.bias")
clusterExport(cl,"coverage")
clusterExport(cl,"MSE")
clusterExport(cl,"model_bd")
clusterEvalQ(cl, {
  library(mice)
  library(mvtnorm)
  library(MASS)
  library(corpcor)
  library(oglmx)
  library(maxLik)
  source("mice_heckman_ordered.R")
  source("function_MI_analysis.R")
})

start_time <- Sys.time()
sim <- 
  parLapply(cl = cl, 1:1000, fun = function(no, mech="non-heckman", rho=rho, M=10){
    
    n <- 2000
    
    x1 <- rnorm(n,0,0.3)
    x2 <- rnorm(n,0,0.8)
    x3 <- rnorm(n,0,4)
    
    ## Missing-data mechanisms ##
    
    if(mech=="heckman") {
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s
      r <- as.factor(ifelse(r<0,0,1))
      
      y_star <- 1*x1+0.5*x2+eps_o
      z <- c(-0.75,0.5)
      
      z1 <- z[1]
      z2 <- z[2]
  
      y <- as.factor( ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
    
    } else if(mech=="non-heckman"){
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      y_star <- 1*x1+0.5*x2+eps_o
      z <- c(-0.75,0.5)
      
      z1 <- z[1]
      z2 <- z[2]
      
      y <- as.factor( ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))

      prob <- pnorm(1.25+1.75*y_star+1.5*x1-2.5*x2)
      r <- rbinom(n,1,p=prob)
      
    } else{ #MAR
      
      vc <- diag(2)
      vc[2,1] <- vc[1,2] <- rho
      eps <- rmvnorm(n,rep(0,2),vc)
      eps_o <- eps[,1]
      eps_s <- eps[,2]
      
      r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s
      r <- as.factor(ifelse(r<0,0,1))
      
      y_star <- 1*x1+0.5*x2+eps_o
      z <- c(-0.75,0.5)
      
      z1 <- z[1]
      z2 <- z[2]
      
      y <- as.factor( ifelse(y_star <= z1 ,0,ifelse(y_star > z1 & y_star <= z2 , 1 ,ifelse(y_star > z2 , 2 ,NA))))
      
    }
      
    ## BD ##
    
    data_comp <- data.frame(y,r,x1,x2,x3)
    
    bd <-  oglmx(y ~ x1 + x2, data=data_comp,link="probit",constantMEAN=FALSE,constantSD=FALSE,delta=0,threshparam=NA)
    
    CI_bd <- confint(bd)
    
    cov_bd <- NULL
    for(j in 1:2){
      cov_bd <- c(cov_bd,coverage(model_bd[j],CI_bd[j,1],CI_bd[j,2]))
    }
    
    
    ## Generation of missing values ##
    
    y[r==0] <- NA
    
    data <- data.frame(y,r,x1,x2,x3)
    
    
    ### Imputation ###
    
    ini <- mice(data,m=1,maxit=0)
    
    ## MNAR Heckman ##
    pred_MNAR <- ini$pred
    pred_MNAR["y","r"] <- 0
    
    imp_MNAR <- mice(data,m=M,maxit=1,method=c("heckman1step_ord","","","",""),pred=pred_MNAR,print=F,excl="x3")
    glm_MNAR <- with(data=imp_MNAR,exp= oglmx(y ~ x1 + x2,link="probit",constantMEAN=FALSE,constantSD=FALSE,delta=0,threshparam=NA))
    est <- t(sapply(1:M,function(x) coef(glm_MNAR$analyses[[x]])))
    var <- t(sapply(1:M,function(x) diag(solve(-(glm_MNAR$analyses[[x]]$hessian)))))
    pool_glm_MNAR <- MI.analysis(est,var,m=M)
    
    cov_MNAR <- NULL
    for(j in 1:4){
      cov_MNAR <- c(cov_MNAR,coverage(model_bd[j],pool_glm_MNAR[j,2],pool_glm_MNAR[j,3]))
    }
    
    
    ## MAR ##
    pred_MAR <- ini$pred
    pred_MAR[,"r"] <- 0
    
    imp_MAR <- mice(data,m=M,maxit=1,method=c("polr","","","",""),pred=pred_MAR,print=F)
    glm_MAR <- with(data=imp_MAR,exp=oglmx(y ~ x1 + x2, link="probit",constantMEAN=FALSE,constantSD=FALSE,delta=0,threshparam=NA))
    est <- t(sapply(1:M,function(x) coef(glm_MAR$analyses[[x]])))
    var <- t(sapply(1:M,function(x) diag(solve(-(glm_MAR$analyses[[x]]$hessian)))))
    pool_glm_MAR <- MI.analysis(est,var,m=M)
    
    cov_MAR <- NULL
    for(j in 1:4){
      cov_MAR <- c(cov_MAR,coverage(model_bd[j],pool_glm_MAR[j,2],pool_glm_MAR[j,3]))
    }
    
    
    ## CC ##
    cc <- oglmx(y ~ x1 + x2, data=data,link="probit",constantMEAN=FALSE, constantSD=FALSE,delta=0,threshparam=NA)
    CI_cc <- confint(cc)
    
    cov_CC <- NULL
    for(j in 1:4){
      cov_CC <- c(cov_CC,coverage(model_bd[j],CI_cc[j,1],CI_cc[j,2]))
    }
    
    return(cbind(pool_glm_MAR[1:2,1], pool_glm_MNAR[,1], coef(cc),coef(bd),
                 rel.bias(model_bd,pool_glm_MAR[,1]),rel.bias(model_bd,pool_glm_MNAR[,1]), rel.bias(model_bd,coef(cc)),
                 rel.bias(model_bd,coef(bd)),
                 cov_MAR,cov_MNAR,cov_CC,cov_bd,
                 MSE(model_bd,pool_glm_MAR[,1]), MSE(model_bd,pool_glm_MNAR[,1]), MSE(model_bd,coef(cc)),
                 MSE(model_bd,coef(bd))))
    
  })
end_time <- Sys.time()
end_time - start_time

stopCluster(cl)

est_MAR <- rowMeans(sapply(sim,function(x) x[,1]))
est_MNAR <- rowMeans(sapply(sim,function(x) x[,2]))
est_CC <- rowMeans(sapply(sim,function(x) x[,3]))
est_bd <- rowMeans(sapply(sim,function(x) x[,4]))

bias_MAR <- rowMeans(sapply(sim,function(x) x[,5]))
bias_MNAR <- rowMeans(sapply(sim,function(x) x[,6]))
bias_CC <- rowMeans(sapply(sim,function(x) x[,7]))
bias_bd <- rowMeans(sapply(sim,function(x) x[,8]))

cov_MAR <- rowMeans(sapply(sim,function(x) x[,9]))
cov_MNAR <- rowMeans(sapply(sim,function(x) x[,10]))
cov_CC <- rowMeans(sapply(sim,function(x) x[,11]))
cov_bd <- rowMeans(sapply(sim,function(x) x[,12]))

mse_MAR <- rowMeans(sapply(sim,function(x) x[,13]))
mse_MNAR <- rowMeans(sapply(sim,function(x) x[,14]))
mse_CC <- rowMeans(sapply(sim,function(x) x[,15]))
mse_bd <- rowMeans(sapply(sim,function(x) x[,16]))  

res <- cbind(est_MAR, est_MNAR, est_CC, est_bd, bias_MAR, bias_MNAR, bias_CC, bias_bd, cov_MAR, cov_MNAR, cov_CC, cov_bd, mse_MAR, mse_MNAR, mse_CC, mse_bd)
