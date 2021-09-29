
MI.analysis <- function(Q.hat,U.hat,m){
  if (class(Q.hat)=="matrix") {
    Q.bar <- colSums(Q.hat)/m
    U.bar <- colSums(U.hat)/m
    B <- colSums((Q.hat - matrix(1, nrow=m, ncol=1) %*% Q.bar)^2)/(m-1)       
  }
  else{
    Q.bar <- sum(Q.hat)/m
    U.bar <- sum(U.hat)/m
    B <- (1/(m-1))*sum((Q.hat-Q.bar)^2) 
  }
  
  T <- U.bar+B+B/m
  df <- (m-1)*(1+(m/(m+1))*U.bar/B)^2
  CIlow <- Q.bar-qt(0.975,df)*sqrt(T)
  CIupper <- Q.bar+qt(0.975,df)*sqrt(T)
  r <- (B+B/m)/U.bar
  FMI <- (r+2/(df+3))/(1+r)
  
  return(cbind(Q.bar,CIlow,CIupper,T))
}

