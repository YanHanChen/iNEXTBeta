boot_beta_one <- function (data)
{
  poolDa <- rowSums(data)
  Sobs <- sum(poolDa > 0)
  n <- sum(poolDa)
  f1 <- sum(poolDa == 1)
  f2 <- sum(poolDa == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
  apply(data, 2, function(x){
    newp <- chaoUtility:::bootp_one_abu(Spec = x,zero = T)
    if(length(newp) != length(x)){
      p0 <- newp[length(newp)]
      out <- newp[1:length(x)]
      out <- c(out,rep(0,f0.hat))
      zrs <- which(out==0)
      chosen <- sample(x = zrs,size = min(length(newp) - length(x),length(zrs)),replace = F)
      out[chosen] <- (1-sum(out))/length(chosen)

      # out <- out/sum(out)
      out
    }else{
      out <- c(newp,rep(0,f0.hat))
      out
    }
  })
}

Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 1:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      #ifelse(A==0,NA,A^(1/(1-q)))
      A^(1/(1-q))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}
