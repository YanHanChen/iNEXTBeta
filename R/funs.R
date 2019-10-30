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
#' @return a list consists of N elements, where N is the number of region. Each element is a list containing 5 tables of gamma. alpha and beta diversities,
#' and 2 types of dissimialrities (Sorensen and Jaccard), respectively.
#' @examples
#' data(abundata)
#' out <- iNEXT_beta(x = abundata, q = c(0,1,2), datatype = "abundance",se = TRUE)
#' @import dplyr
#' @importFrom iNEXT iNEXT
#' @export
iNEXT_beta <- function(x, q = c(0,1,2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, se = TRUE, conf = 0.95,
                       nboot = 30){
  if(class(x)=="data.frame" | class(x)=="matrix" ){
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
    C_ <- beta %>% mutate(C_ = ifelse(order==1,log(qD)/log(N_site),(qD^(1-order) - 1)/(N_site^(1-order)-1))) %>%
      select(m, method, order, C_, SC)
    U_ <- beta %>% mutate(U_ = ifelse(order==1,log(qD)/log(N_site),(qD^(order-1) - 1)/(N_site^(order-1)-1))) %>%
      select(m, method, order, U_, SC)

    if(se==TRUE & nboot>1){
      boot_pop <- boot_beta_one(data)
      ses <- sapply(1:nboot, function(i){
        x_bt <- sapply(1:ncol(data),function(k) rmultinom(n = 1,size = sum(data[,k]),prob = boot_pop[,k]))
        x_bt <- x_bt[rowSums(x_bt)>0,]
        ga_da_bt <- rowSums(x_bt)
        al_da_bt <- as.numeric(x_bt)
        al_da_bt <- x_bt[x_bt>0]
        gamma_bt <- iNEXT(x = ga_da_bt,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                          se = FALSE,conf = conf,nboot=nboot)$iNextEst %>% select(order,qD,SC)
        qs <- gamma_bt$order
        gamma_bt <- gamma_bt %>% select(-order)
        alpha_bt <- iNEXT(x = al_da_bt,q = q,datatype = datatype,size = size,endpoint = endpoint,knots = knots,
                          se = FALSE,conf = conf,nboot=nboot)$iNextEst %>% select(qD,SC)
        alpha_bt$qD <- alpha_bt$qD/N_site
        beta_bt <- gamma_bt$qD/alpha_bt$qD
        beta_bt <- data.frame(q = qs, beta.qD = beta_bt) %>%
          mutate(C_ = ifelse(q==1,log(beta.qD)/log(N_site),(beta.qD^(1-q) - 1)/(N_site^(1-q)-1) ),
                 U_ = ifelse(q==1,log(beta.qD)/log(N_site),(beta.qD^(q-1) - 1)/(N_site^(q-1)-1) )) %>% select(-q)
        out_bt <- cbind(gamma = gamma_bt,alpha = alpha_bt, beta_bt) %>% as.matrix()
        out_bt
      },simplify = "array") %>% apply(., 1:2, sd) %>% data.frame()
    }else { ses <- data.frame(gamma.qD=rep(0,nrow(gamma)),gamma.SC=rep(0,nrow(gamma)),alpha.qD=rep(0,nrow(gamma)),
                              alpha.SC=rep(0,nrow(gamma)),beta.qD=rep(0,nrow(gamma)), C_=rep(0,nrow(gamma)),
                              U_=rep(0,nrow(gamma)))}
    gamma <- gamma %>% mutate(qD.LCL = qD - tmp * ses$gamma.qD, qD.UCL = qD + tmp * ses$gamma.qD,
                     SC.LCL = SC - tmp * ses$gamma.SC, SC.UCL = SC + tmp * ses$gamma.SC, Region = nm)
    alpha <- alpha %>% mutate(qD.LCL = qD - tmp * ses$alpha.qD, qD.UCL = qD + tmp * ses$alpha.qD,
                     SC.LCL = SC - tmp * ses$alpha.SC, SC.UCL = SC + tmp * ses$alpha.SC, Region = nm)
    beta <- beta %>% mutate(qD.LCL = qD - tmp * ses$beta.qD, qD.UCL = qD + tmp * ses$beta.qD,
                     SC.LCL = gamma$SC.LCL, SC.UCL = gamma$SC.UCL, Region = nm)
    C_ <- C_ %>% mutate(diss.LCL = C_ - tmp * ses$C_, diss.UCL = C_ + tmp * ses$C_,
                        SC.LCL = gamma$SC.LCL, SC.UCL = gamma$SC.UCL, Region = nm)
    U_ <- U_ %>% mutate(diss.LCL = U_ - tmp * ses$U_, diss.UCL = U_ + tmp * ses$U_,
                        SC.LCL = gamma$SC.LCL, SC.UCL = gamma$SC.UCL, Region = nm)
    list(gamma = gamma, alpha = alpha, beta = beta, Sorensen = C_, Jaccard = U_)
  }
  out <- lapply(1:length(x), function(i) eachR(data = mydata[[i]],nm = nms[i],N_site = n_sp[i]))
  names(out) <- nms
  out
}


#' \code{ggiNEXT_beta}: plot the outcome of \code{iNEXT_beta} based on the \code{ggiNEXT} function.
#' @param x the outcome of \code{iNEXT_beta}
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (type = 1); sample completeness curve (type = 2);
#' coverage-based rarefaction/extrapolation curve (type = 3).
#' @return a list containing 2 ggplot2 object, \code{$betaDiv} for beta diversity and \code{$diss} for dissimilarity
#' @examples
#' data(abundata)
#' out <- iNEXT_beta(x = abundata, q = c(0,1,2), datatype = "abundance")
#' ggout1 <- ggiNEXT_beta(output = out, type = 1)
#' ggout2 <- ggiNEXT_beta(output = out, type = 2)
#' ggout3 <- ggiNEXT_beta(output = out, type = 3)
#' @import dplyr
#' @import ggplot2
#' @export
ggiNEXT_beta <- function(output, type = 1, goal = "Diversity"){
  gamma <- lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
  alpha <- lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
  beta <- lapply(output, function(y) y[["beta"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Beta") %>% as_tibble()
  a <- rbind(gamma,alpha,beta)
  if(type == 1){
    z <- a %>% select(-c(SC,SC.LCL,SC.UCL)) %>% rename(x=m,y=qD,y.LCL=qD.LCL,y.UCL=qD.UCL)
    labx <- "Number of individual"
    laby <- "Diveristy"
  }else if (type==2){
    z <- a %>% select(-c(qD,qD.LCL,qD.UCL)) %>% rename(x=m,y=SC,y.LCL=SC.LCL,y.UCL=SC.UCL) %>%
      filter(order == unique(order)[1])
    labx <- "Number of individual"
    laby <- "Coverage"
  }else if (type==3){
    z <- a %>% select(-c(m,SC.LCL,SC.UCL)) %>% rename(x=SC,y=qD,y.LCL=qD.LCL,y.UCL=qD.UCL)
    labx <- "Coverage"
    laby <- "Diveristy"
  }
  z$div_type <- factor(z$div_type,levels = c("Gamma","Alpha","Beta"))
  z_IE <- z[z$method!="observed",]
  z_IE$method <- factor(z_IE$method,levels = c("interpolated","extrapolated"))
  p <- ggplot(data = z,aes(x = x,y = y, col = Region)) + geom_line(data = z_IE,aes(linetype = method),size = 1.2)+
    facet_grid(div_type~order,scales = "free_y")+
    geom_ribbon(aes(ymin=y.LCL, ymax=y.UCL, fill=Region, colour=NULL), alpha=0.3)+
    geom_point(data = z[z$method=="observed",],aes(shape=Region),size=3)+
    theme_bw()+theme(legend.position = "bottom",legend.title = element_blank())+xlab(labx)+ylab(laby)

  C_ <- lapply(output, function(y) y[["Sorensen"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Sorensen") %>%
    as_tibble() %>% rename(dissimilarity = C_)
  U_ <- lapply(output, function(y) y[["Jaccard"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Jaccard") %>%
    as_tibble() %>% rename(dissimilarity = U_)
  b <-rbind(C_,U_)
  if(type == 1){
    z <- b %>% select(-c(SC,SC.LCL,SC.UCL)) %>% rename(x=m,y=dissimilarity,y.LCL=diss.LCL,y.UCL=diss.UCL)
    labx <- "Number of individual"
    laby <- "Dissimilarity"
  }else if (type==3){
    z <- b %>% select(-c(m,SC.LCL,SC.UCL)) %>% rename(x=SC,y=dissimilarity,y.LCL=diss.LCL,y.UCL=diss.UCL)
    labx <- "Coverage"
    laby <- "Dissimilarity"
  }
  z$div_type <- factor(z$div_type,levels = c("Sorensen","Jaccard"))
  z_IE <- z[z$method!="observed",]
  z_IE$method <- factor(z_IE$method,levels = c("interpolated","extrapolated"))
  p2 <- ggplot(data = z,aes(x = x,y = y, col = Region)) + geom_line(data = z_IE,aes(linetype = method),size = 1.2)+
    facet_grid(div_type~order,scales = "free_y")+
    geom_ribbon(aes(ymin=y.LCL, ymax=y.UCL, fill=Region, colour=NULL), alpha=0.3)+
    geom_point(data = z[z$method=="observed",],aes(shape=Region),size=3)+
    theme_bw()+theme(legend.position = "bottom",legend.title = element_blank())+xlab(labx)+ylab(laby)
  return(list(betaDiv = p, diss = p2))
}


