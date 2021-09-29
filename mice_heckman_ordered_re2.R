########################################################################
########################################################################
##
## MI function for imputing MNAR hierarchical ordered data (AGHQ)
## (usable with R package mice)
##
## needs installation and loading of the bivprob.quad package
##
## @author Angelina Hammon
##
## 29.09.2021
##
########################################################################
########################################################################


# loading necessary packages:
require(MASS)
require(mvtnorm)
require(corpcor)
require(maxLik)
require(compiler)
require(Matrix)
require(bivprob.quad)
require(VGAM)

# function taken from mice package from Stef van Buuren: 
mice_multilevel_imputation_draw_random_effects <- function( mu, Sigma,ridge=1E-20 ){  
  dim_Sigma <- dim(Sigma)
  ngr <- dim_Sigma[3]
  NR <- dim_Sigma[1]
  # draw random effects
  u <- matrix(0, nrow=ngr, ncol=NR)
  for(i in 1:ngr){
    #-- compute covariance matrix with ridge
    Sigma1 <- Sigma[,,i] + diag(ridge,NR)
    # Cholesky matrix of Sigma1
    Sigma_chol <- chol(Sigma1)
    # simulation
    rn <- stats::rnorm(NR, mean=0, sd=1)
    u[i,] <- mu[i,] + as.vector( t( Sigma_chol ) %*% rn )
  }
  return(u)
}



# censored bivariate probit model with random intercept using AGHQ:
mice.impute.2l.heckman1step_ord_re_aghq <- function(y, ry, x, type, QP=rep(10,2), draw=T, excl="",...) { 
  
  nam <- colnames(x) 
  
  group <- x[,type==-2]  
  gr_unique <- unique(group)
  gr_name <- nam[type==-2]
  
  type <- type[-which(nam==gr_name)]
  nam <- nam[-which(nam==gr_name)]
  
  pos.exc <- grep(excl,nam)
  
  ind_s <- length(nam)+1   
  ind_y <- length(nam)-length(pos.exc)
  
  ry <- as.numeric(ry)
  data <- data.frame(x,y,ry)
  colnames(data) <- c(colnames(x),"y","ry")
  
  
  ## Calculate one-step heckman model considering random intercept ##
  exc <- pos.exc 
  
  sel <- paste(nam, collapse="+")
  sel <- paste(sel,"+(1|",gr_name,")",sep="")
  sel <- as.formula(paste("ry","~",sel,sep=""))
  out <- paste(nam[-exc], collapse="+")
  out <- paste(out,"+(1|",gr_name,")",sep="")
  out <- as.formula(paste("y","~",out,sep=""))
  
  XS <- data.matrix(cbind(rep(1,length(y)),x[,-which(colnames(x)==gr_name)]))
  XO <- data.matrix(x[,-c(which(colnames(x)==gr_name),grep(excl,colnames(x)))])
  
  # estimating the model:
  heckit <- loglik_ordbivprob_est_aghq(data=data,yS=ry, yO=y, XO=XO, XS=XS, kS=ind_s, kO=ind_y,sel = sel, out=out, group=group, QP=QP) 
  
  phi <- as.matrix(nearPD(ginv(heckit$hessian))$mat)
 
  q_hat <- heckit$estimates
  
  linpredS_int <- XS %*% q_hat[1:ind_s]  
  linpredO_int <- XO %*% q_hat[(ind_s+1):(ind_s+ind_y)] 
  
  group_aux <- c(0,cumsum(table(as.integer(factor(data[,gr_name])))))
  
  y <- as.numeric(as.character(y))
  levels <- as.numeric(as.character(unique(y[order(y)])))
  
  fy <- as.factor(y)
  nc <- length(levels(fy))
  
  z <-  q_hat[(ind_s+ind_y+2):(ind_s+ind_y+2+nc-2)]
  for(i in 2:length(z)){
    z[i] <- exp(z[i])+z[i-1]
  }
  
  # calculate mode of random intercepts: 
  start2 <- c(0,0)
  modes <- t(sapply(1:length(table(data[,gr_name])), function(x) {
    optim(par=start2, fn=mode_ord, gr=NULL, method="Nelder-Mead", x=x, yS=ry, yO=y, linpredS=linpredS_int, linpredO=linpredO_int, 
          group_aux=group_aux,rho=tanh(q_hat[ind_s+ind_y+1]),z=z,tau=tanh(q_hat[length(q_hat)]) ,
          sigma1=exp(q_hat[ind_s+ind_y+2+nc-2+1]),sigma2=exp(q_hat[ind_s+ind_y+nc-2+2]),levels=levels)$par
  })) 
  
  # hessian of mode:
  hess_mode <- lapply(1:nrow(modes), function(x) { 
    numericHessian(f=mode_ord, grad=NULL, t0=modes[x,], x=x, yS=ry, yO=y, linpredS=linpredS_int, linpredO=linpredO_int, group_aux=group_aux,
                   rho=tanh(q_hat[ind_s+ind_y+1]),z=z,tau=tanh(q_hat[length(q_hat)]) ,
                   sigma1=exp(q_hat[ind_s+ind_y+2+nc-2+1]),sigma2=exp(q_hat[ind_s+ind_y+nc-2+2]),levels=levels)
  })
  
  
  if(any(sapply(1:nrow(modes),function(x) isSymmetric(hess_mode[[x]])) == F)) {
    hess_mode <- lapply(1:nrow(modes),function(x) as.matrix(forceSymmetric(hess_mode[[x]])))
  }
  
  cov_mode <- lapply(1:nrow(modes),function(x) make.positive.definite(ginv(hess_mode[[x]])))
  
  NR <- ncol(modes)
  ngr <- length(gr_unique)
  re <- matrix(0, nrow=ngr, ncol=NR)     
  re <- as.matrix(modes)        
  
  pv <- array(0, dim=c(NR,NR,ngr))
  
  for(i in 1:ngr) {
    pv[,,i] <- cov_mode[[i]]  
  }                                    
  
  
  ## Draw ALL parameter values ##
  if(draw==T){
    q_star <- rmvnorm(1,q_hat,phi,method = "chol")
    beta_s_star <- q_star[1:ind_s]          
    beta_star <- q_star[(ind_s+1):(ind_s+ind_y)] 
    # retransform parameters:
    rho_star <- tanh(q_star[ind_s+ind_y+1])
    z_star <-  q_star[(ind_s+ind_y+2):(ind_s+ind_y+2+nc-2)]
    for(i in 2:length(z_star)){
      z_star[i] <- exp(z_star[i])+z_star[i-1]
    }
    sigma1_star <- exp(q_star[(ind_s+ind_y+2+nc-2+1)])
    sigma2_star <- exp(q_star[(ind_s+ind_y+2+nc-2+2)])
    tau_star <- tanh(q_star[length(q_star)]) 
  }else{ ## just draw betas ##
    q_star <- rmvnorm(1,q_hat,phi,method = "chol")
    beta_s_star <- q_star[1:ind_s]          
    beta_star <- q_star[(ind_s+1):(ind_s+ind_y)] 
    rho_star <- tanh(q_hat[ind_s+ind_y+1])
    z_star <-  q_hat[(ind_s+ind_y+2):(ind_s+ind_y+2+nc-2)]
    for(i in 2:length(z_star)){
      z_star[i] <- exp(z_star[i])+z_star[i-1]
    }
    sigma1_star <- exp(q_hat[(ind_s+ind_y+2+nc-2+1)])
    sigma2_star <- exp(q_hat[(ind_s+ind_y+2+nc-2+2)])
    tau_star <- tanh(q_hat[length(q_hat)])
  }
  
  
  ## Calculate prob based on the new parameter values ## 
  
  # draw random intercepts:
  nj <- as.vector(table(group))
  
  alpha <- mice_multilevel_imputation_draw_random_effects( mu=re, Sigma=pv, ridge=1E-6)
  alpha_s <- rep(alpha[,1], times=nj)  
  alpha_o <- rep(alpha[,2], times=nj) 
  
  # calculate linear predictors considering random intercepts: 
  lin.predS <- XS %*% beta_s_star + alpha_s
  lin.predO <- XO %*% beta_star + alpha_o   
  
  lin.predS <- as.numeric(lin.predS)[!ry]
  lin.predO <- as.numeric(lin.predO)[!ry]
  
  z_star <- c(-Inf,z_star,+Inf)
  
  p <- sapply(1:nc, function(j){ sapply(1:sum(!ry),function(i){ 
      (pbinorm(-lin.predS[i],(z_star[j+1]-lin.predO[i]), cov12 = rho_star) - 
       pbinorm(-lin.predS[i],(z_star[j]-lin.predO[i]), cov12 = rho_star))/pnorm(-lin.predS[i]) })
  })
  
  
  ## Draw Ymis based on new probability ##
  un <- rep(runif(sum(!ry)), each = nc)
  
  draws <- un > apply(p, 1, cumsum)
  idx <- 1 + apply(draws, 2, sum)
  
  vec <- levels(fy)[idx]
  
  return(vec)
}

mice.impute.2l.heckman1step_ord_re_aghq <- compiler::cmpfun(mice.impute.2l.heckman1step_ord_re_aghq)
