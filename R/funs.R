#' \code{iNEXT_beta}: compute the rarefaction and extrapolation of beta diveristy.
#' @param x a list consist of N data.frame/matrix describing species-by-assemblage/plot abundance. Note that
#' the species in each element must exactly match including specpes order. Use \code{data(abundata)} to see data example.
#' @param q a numeric value or vector specifying the diversity order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed.
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default endpoint and knots.
#' @param endpoint an integer specifying the sample size that is the endpoint for rarefaction/extrapolation.
#' If NULL, then endpoint = double reference sample size.
#' @param knots an integer specifying the number of equally-spaced knots (say K, default is 40) between size 1 and the endpoint;
#' each knot represents a particular sample size for which diversity estimate will be calculated. If the endpoint is smaller than the reference sample size,
#' then \code{iNEXT_beta()} computes only the rarefaction esimates for approximately K evenly spaced knots. If the endpoint is larger than the reference sample size,
#' then \code{iNEXT_beta()} computes rarefaction estimates for approximately K/2 evenly spaced knots between sample size 1 and the reference sample size,
#' and computes extrapolation estimates for approximately K/2 evenly spaced knots between the reference sample size and the endpoint.
#' @param se a logical variable to calculate the bootstrap standard error and conf confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param nboot an integer specifying the number of replications.
#' @return a list of 3 iNEXT objects, representing gamma, alpha and beta diversity respectively.
#' the profiles of all six classes of evenness indices listed in Chao and Ricotta (2019) Ecology paper.
#' @examples
#' data(abundata)
#' out <- iNEXT_beta(x = abundata, q = c(0,1,2), datatype = "abundance")
#' @import dplyr
#' @importFrom iNEXT iNEXT
#' @export
iNEXT_beta <- function(x, q = c(0,1,2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, se = FALSE, conf = 0.95,
                       nboot = 30){
  if(class(x)=="data.frame" & datatype =="abundance"){
    x <- list(region1 = x)
  }
  # if some community contains zero species?
  if(class(x)== "list"){
    if(is.null(names(x))) nms <- paste0("region",1:length(x))else nms <- names(x)
    n_sp <- sapply(x,ncol)
    mydata <- lapply(x, function(y) y[rowSums(y)>0,])
  }
  if(is.null(conf)) conf <- 0.95
  tmp <- qnorm(1 - (1 - conf)/2)
  eachR <- function(data,nm,N_site){
    gdata <- rowSums(data)
    adata <- as.matrix(data) %>% as.numeric
    adata <- adata[adata>0]
    gamma <- iNEXT(x = gdata,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                   se = FALSE,conf = conf,nboot=nboot)$iNextEst
    alpha <- iNEXT(x = adata,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                   se = FALSE,conf = conf,nboot=nboot)$iNextEst
    alpha$qD <- alpha$qD/N_site
    beta <- gamma
    beta$qD <- gamma$qD/alpha$qD

    if(se==TRUE & nboot>1){
      boot_pop <- boot_beta_one(data)
      ses <- sapply(1:nboot, function(i){
        x_bt <- sapply(1:ncol(data),function(k) rmultinom(n = 1,size = sum(data[,k]),prob = boot_pop[,k]))
        x_bt <- x_bt[rowSums(x_bt)>0,]
        ga_da_bt <- rowSums(x_bt)
        al_da_bt <- as.numeric(x_bt)
        al_da_bt <- x_bt[x_bt>0]
        gamma_bt <- iNEXT(x = ga_da_bt,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                          se = FALSE,conf = conf,nboot=nboot)$iNextEst %>% select(qD,SC)
        alpha_bt <- iNEXT(x = al_da_bt,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                          se = FALSE,conf = conf,nboot=nboot)$iNextEst %>% select(qD,SC)
        alpha_bt$qD <- alpha_bt$qD/N_site
        beta_bt <- gamma_bt$qD/alpha_bt$qD
        out_bt <- cbind(gamma = gamma_bt,alpha = alpha_bt, beta.qD = beta_bt) %>% as.matrix()
        out_bt
      },simplify = "array") %>% apply(., 1:2, sd) %>% data.frame()
    }else { ses <- data.frame(gamma.qD=rep(0,nrow(gamma)),gamma.SC=rep(0,nrow(gamma)),alpha.qD=rep(0,nrow(gamma)),
                              alpha.SC=rep(0,nrow(gamma)),beta.qD=rep(0,nrow(gamma)))}
    gamma <- gamma %>% mutate(qD.LCL = qD - tmp * ses$gamma.qD, qD.UCL = qD + tmp * ses$gamma.qD,
                     SC.LCL = SC - tmp * ses$gamma.SC, SC.UCL = SC + tmp * ses$gamma.SC, Region = nm)
    alpha <- alpha %>% mutate(qD.LCL = qD - tmp * ses$alpha.qD, qD.UCL = qD + tmp * ses$alpha.qD,
                     SC.LCL = SC - tmp * ses$alpha.SC, SC.UCL = SC + tmp * ses$alpha.SC, Region = nm)
    beta <- beta %>% mutate(qD.LCL = qD - tmp * ses$beta.qD, qD.UCL = qD + tmp * ses$beta.qD,
                     SC.LCL = SC - tmp * ses$gamma.SC, SC.UCL = SC + tmp * ses$gamma.SC, Region = nm)
    list(gamma = gamma, alpha = alpha, beta = beta)
  }
  test <- lapply(1:length(x), function(i) eachR(data = mydata[[i]],nm = nms[i],N_site = n_sp[i]))
}


gamma <- iNEXT(x = gamma_data,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
               se = FALSE,conf = conf,nboot=nboot)
alpha <- iNEXT(x = alpha_data,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
               se = FALSE,conf = conf,nboot=nboot)
beta <- gamma
if(se==TRUE & nboot>1){
  boot_pop <- lapply(1:length(x),function(i) boot_beta_one(data = x[[i]]))
  bt <- lapply(1:nboot, function(i){
    x_bt <- sapply(1:length(x), function(j){
      sapply(1:ncol(x[[j]]),function(k) rmultinom(n = 1,size = sum(x[[j]][,k]),prob = boot_pop[[j]][,k]))
    })
    ga_da_bt <- lapply(x_bt, rowSums)
    ga_da_bt <- lapply(ga_da_bt, function(i) i[i>0])
    names(ga_da_bt) <- nms
    al_da_bt <- lapply(x_bt, function(i) as.matrix(i) %>% as.numeric)
    al_da_bt <- lapply(al_da_bt, function(i) i[i>0])
    names(al_da_bt) <- nms
    gamma_bt <- iNEXT(x = ga_da_bt,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                      se = FALSE,conf = conf,nboot=nboot)$iNextEst %>% lapply(.,function(y) y$qD)
    alpha_bt <- iNEXT(x = al_da_bt,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                      se = FALSE,conf = conf,nboot=nboot)$iNextEst %>% lapply(.,function(y) y$qD)
    out_bt <- lapply(1:length(x),function(j){
      beta_bt <- n_sp[j]*gamma_bt[[j]]/alpha_bt[[j]]
      cbind(gamma = gamma_bt[[j]],alpha = alpha_bt[[j]]/n_sp[j], beta = beta_bt)
    })
    out_bt
  })
  se <- lapply(1:length(x), function(i) {
    sapply(bt, function(y) y[[i]] ,simplify = "array") %>% apply(., 1:2, sd)
  })
}else { se <- 0}
for(i in 1:length(nms)){
  tmp <- qnorm(1 - (1 - conf)/2)
  alpha$iNextEst[[i]]$qD <- alpha$iNextEst[[i]]$qD/n_sp[i]
  beta$iNextEst[[i]]$qD <- gamma$iNextEst[[i]]$qD/alpha$iNextEst[[i]]$qD

  gamma$iNextEst[[i]]$qD.LCL <- gamma$iNextEst[[i]]$qD - tmp * se[[i]]$gamma
  gamma$iNextEst[[i]]$qD.UCL <- gamma$iNextEst[[i]]$qD + tmp * se[[i]]$gamma
  alpha$iNextEst[[i]]$qD.LCL <- alpha$iNextEst[[i]]$qD - tmp * se[[i]]$alpha
  alpha$iNextEst[[i]]$qD.UCL <- alpha$iNextEst[[i]]$qD + tmp * se[[i]]$alpha
  beta$iNextEst[[i]]$qD.LCL <- beta$iNextEst[[i]]$qD - tmp * se[[i]]$beta
  beta$iNextEst[[i]]$qD.UCL <- beta$iNextEst[[i]]$qD + tmp * se[[i]]$beta

}
return(list(gamma = gamma, alpha = alpha, beta = beta))

#' \code{ggiNEXT_beta}: plot the outcome of \code{iNEXT_beta} based on the \code{ggiNEXT} function.
#' @param x the outcome of \code{iNEXT_beta}
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (type = 1); sample completeness curve (type = 2);
#' coverage-based rarefaction/extrapolation curve (type = 3).
#' @param se a logical variable to calculate the bootstrap standard error and conf confidence interval.
#' @return a list containing 3 ggplot2 object
#' @examples
#' data(abundata)
#' out <- iNEXT_beta(x = abundata, q = c(0,1,2), datatype = "abundance")
#' ggout1 <- ggiNEXT_beta(x = out, type = 1)
#' ggout2 <- ggiNEXT_beta(x = out, type = 2)
#' ggout3 <- ggiNEXT_beta(x = out, type = 3)
#' @importFrom iNEXT ggiNEXT
#' @export
ggiNEXT_beta <- function(x, type = 1, se = FALSE){
  gamma_p <- ggiNEXT(x$gamma,facet.var = "order",type = type,se = F)+ggtitle("Gamma diversity")
  alpha_p <- ggiNEXT(x$alpha,facet.var = "order",type = type,se = F)+ggtitle("Alpha diversity")
  beta_p <- ggiNEXT(x$beta,facet.var = "order",type = type,se = F)+ggtitle("Beta diversity")
  return(list(gamma_p,alpha_p,beta_p))
}


