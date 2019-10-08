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
iNEXT_beta <- function(x, q = 0, datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, se = TRUE, conf = 0.95,
                       nboot = 30){
  if(class(x)== "list"){
    if(is.null(x)) nms <- paste0(region,1:length(x))else nms <- names(x)
    n_sp <- sapply(x,ncol)
    gamma_data <- lapply(x, rowSums) %>% do.call(cbind,.)
    gamma_data <- gamma_data[rowSums(gamma_data)>0,]
    colnames(gamma_data) <- nms
    alpha_data <- lapply(x, as.numeric) %>% do.call(cbind,.)
    alpha_data <- alpha_data[rowSums(alpha_data)>0,]
    colnames(alpha_data) <- nms
  }
  gamma <- iNEXT(x = gamma_data,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = 80,
                 se = se,conf = conf,nboot=nboot)
  alpha <- iNEXT(x = alpha_data,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = 80,
                 se = se,conf = conf,nboot=nboot)
  beta <- gamma
  for(i in 1:length(nms)){
    alpha$iNextEst[[i]]$qD <- alpha$iNextEst[[i]]$qD/n_sp[i]
    beta$iNextEst[[i]]$qD <- gamma$iNextEst[[i]]$qD/alpha$iNextEst[[i]]$qD
  }
  return(list(gamma = gamma, alpha = alpha, beta = beta))
}

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


