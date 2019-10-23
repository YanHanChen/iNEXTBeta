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

iNEXT_b_abu <- function(data)

function (Spec, zero = TRUE)
{
  Sobs <- sum(Spec > 0)
  n <- sum(Spec)
  f1 <- sum(Spec == 1)
  f2 <- sum(Spec == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2)
  A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
  a <- f1/n * A
  b <- sum(Spec/n * (1 - Spec/n)^n)
  if (f0.hat == 0) {
    w <- 0
    if (sum(Spec > 0) == 1) {
      warning("This site has only one species. Estimation is not robust.")
    }
  }else {
    w <- a/b
  }
  if (zero == FALSE)
    Spec <- Spec[Spec > 0]
  Prob.hat <- Spec/n * (1 - w * (1 - Spec/n)^n)
  Prob.hat.Unse <- rep(a/ceiling(f0.hat), ceiling(f0.hat))
  return(c(Prob.hat, Prob.hat.Unse))
}

sapply(1:ncol(data), function(i){
  x <- data[,i]
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
    out <- c(x/sum(x),rep(0,f0.hat))
    out
  }
  print(i)
})
