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
      out[chosen] <- p0
      out <- out/sum(out)
      out
    }else{
      out <- c(newp,rep(0,f0.hat))
      out
    }
  })
}

