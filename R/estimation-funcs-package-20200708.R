dimm <- function(x, ...) UseMethod("dimm")

print.dimm <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nVariance:\n")
  print(x$vcov)
}

summary.dimm <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.dimm"
  return(res) 
}

print.summary.dimm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}

objective.dimm <- function(object, ...){
  return(objective.dimm.compute(object))
}

objective.dimm.compute <- function(object, ...){
  p <- dim(object$covariates)[2]-1
  full.coefficients <- object$full.coefficients
  if (object$analysis=="PCL" & object$family=="gaussian"){
    func1 <- function(beta, cov, block_y, block_x, id, m, n){logCLnormal_par(beta, cov, block_y, block_x, m, n)}
    if (object$corstr=="CS"){
      func4 <- function(beta, sigma, rho, block_y, block_x, id, m, n){lapply(eenormalCSvar_par(beta, sigma, rho, block_y, block_x, m, n), function(x) x/n)}
      d <- 2
    }
    if (object$corstr=="AR1"){
      func4 <- function(beta, sigma, rho, block_y, block_x, id, m, n){lapply(eenormalAR1var_par(beta, sigma, rho, block_y, block_x, m, n), function(x) x/n)}
      d <- 2
    }
    if (object$corstr=="independence"){
      func4 <- function(beta, sigma, block_y, block_x, id, m, n){lapply(eenormalindvar_par(beta, sigma, block_y, block_x, m, n), function(x) x/n)}
      d <- 1
    }
  }
  if (object$analysis=="GEE"){
    if(object$corstr!="independence"){
      A <- matrix(1:(object$J*object$K*2), 2)
      full.coefficients <- c(object$full.coefficients[1:p],object$full.coefficients[-c(1:p)][c(A[2:1,])]) 
    }
    if(family=="gaussian") {
      family2 <- gaussian()
      if (object$corstr=="CS"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("exch", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (object$corstr=="AR1"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("ar1", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (object$corstr=="independence"){
        d <- 1
        lower_lim <- c(rep(-Inf,p), 1e-5)
        upper_lim <- rep(Inf,p+1)
        func45 <- function(beta, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("independence", CORSTRS, -1)
          alpha <- vector("numeric",0)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+1,p+1)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
    }
    if(family=="binomial") {
      family2 <- binomial()
      if (object$corstr=="CS"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("exch", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (object$corstr=="AR1"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("ar1", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (object$corstr=="independence"){
        d <- 1
        lower_lim <- c(rep(-Inf,p), 1e-5)
        upper_lim <- rep(Inf,p+1)
        func45 <- function(beta, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("independence", CORSTRS, -1)
          alpha <- vector("numeric",0)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+1,p+1)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat <- -Hessmat/n
        }
      } 
    }
  }
  if (object$analysis=="QIF"){
    if (object$corstr=="CS" | object$corstr=="AR1") d <- 2
    if (object$corstr=="independence") d <- 1
    if (object$corstr=="CS+AR1") d <- 3
  }
  obj <- list()
  obj$call <- object$call
  N <- length(object$subject_indicator)
  if(object$analysis=="QIF"){
    response <- cbind(matrix(object$response[,-which(colnames(object$response) == object$id)], 
                             length(object$response_indicator), N), object$response_indicator)
    colnames(response) <- c(object$subject_indicator,"response_indicator")
    Q_rank <- diqif.obj.eval(covariates=object$covariates, response, family=object$family, corstr=object$corstr, 
                        beta=as.vector(object$coefficients), response_indicator=object$response_indicator, 
                        subject_indicator=object$subject_indicator, J=object$J, K=object$K, d, comb_scheme=object$comb_scheme)
  } else {
    response <- cbind(matrix(object$response[,-which(colnames(object$response) == object$id)], 
                             length(object$response_indicator), N), object$response_indicator)
    colnames(response) <- c(object$subject_indicator,"response_indicator")
    Q_rank <- dimm.obj.eval(covariates=object$covariates, response, family=object$family, corstr=object$corstr, 
                            beta=as.vector(object$coefficients), full.coefficients, response_indicator=object$response_indicator, 
                            subject_indicator=object$subject_indicator, J=object$J, K=object$K, d, func=func4, comb_scheme=object$comb_scheme)
  }
  obj$Q <- Q_rank$Q
  obj$df <- Q_rank$rank-length(as.vector(object$coefficients))
  obj$p.value <- 1-as.numeric(pchisq(obj$Q, obj$df))
  class(obj) <- "objective.dimm"
  return(obj)
}

diqif.obj.eval <- function(covariates, response, family, corstr, beta, response_indicator, subject_indicator, J, K, d, comb_scheme){
  p <- dim(covariates)[2]-1
  psi_g <- matrix(0, length(subject_indicator), J*K*p*d)
  psi_g_sum <- vector("numeric", J*K*p*d)
  p_1 <- p*d
  for (k in 1:K) {
    for (j in 1:J) {
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.data.frame(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
      id <- block_x[,which(colnames(block_x)=="id")]
      block_x <- as.matrix(block_x[, -which(colnames(block_x)=="id")])
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      beta_j <- beta[(p*(comb_scheme[((k-1)*J+j)]-1)+1):(p*(comb_scheme[((k-1)*J+j)]-1)+p)]
      if (corstr == "CS") qif_block_fit <- QIF_eval(block_x, c(block_y), nobs=rep(m,n), family=family, corstr="exchangeable", beta_j)
      if (corstr == "AR1") qif_block_fit <- QIF_eval(block_x, c(block_y), nobs=rep(m,n), family=family, corstr="AR-1", beta_j)
      if (corstr == "independence") qif_block_fit <- QIF_eval(block_x, c(block_y), nobs=rep(m,n), family=family, corstr="independence", beta_j)
      if (corstr == "CS+AR1") qif_block_fit <- QIF_eval(block_x, c(block_y), nobs=rep(m,n), family=family, corstr="CS+AR1", beta_j)
      psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
            ((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- t(do.call(cbind, qif_block_fit$gi_list))
      psi_g_sum[((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- qif_block_fit$g_sub
    }
  }
  
  unordered_psi_g <- psi_g
  unordered_psi_g_sum <- psi_g_sum
  psi_g <- unordered_psi_g[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
  psi_g_sum <- unordered_psi_g_sum[order(sapply(comb_scheme, function(x) rep(x,p_1)))]

  rank <- qr(psi_g)$rank
  if(dim(psi_g)[2] != rank){
    psi_g_pca <- prcomp(psi_g, center = FALSE,scale. = FALSE)
    old_psi_g <- psi_g
    old_psi_g_sum <- psi_g_sum
    psi_g <- psi_g_pca$x[,1:rank]
    psi_g_sum <- (psi_g_sum %*% psi_g_pca$rotation)[,1:rank]
  } 
  
  V_psi <- t(psi_g)%*%psi_g
  W <- solve(V_psi)
  rank <- qr(V_psi)$rank
  
  Q <- psi_g_sum %*% W %*% psi_g_sum 
  return(list(Q=Q,rank=rank))
}

dimm.obj.eval <- function(covariates, response, family, corstr, beta, full.coefficients, response_indicator, subject_indicator, J, K, d, func, comb_scheme){
  p <- dim(covariates)[2]-1
  psi_g <- matrix(0, length(subject_indicator), J*K*(p+d))
  psi_g_sum <- vector("numeric", J*K*(p+d))
  for (k in 1:K) {
    for (j in 1:J) {
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.data.frame(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
      id <- block_x[,which(colnames(block_x)=="id")]
      block_x <- as.matrix(block_x[, -which(colnames(block_x)=="id")])
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      beta_j <- beta[(p*(comb_scheme[((k-1)*J+j)]-1)+1):(p*(comb_scheme[((k-1)*J+j)]-1)+p)]
      if(d==2){
        psi_g_combined <- func(beta_j, full.coefficients[(j-1+(k-1)*J)*d+length(beta)+1], full.coefficients[(j-1+(k-1)*J)*d+length(beta)+2], 
                               block_y=block_y, block_x=block_x, id=id, m=m, n=n) 

        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))*n
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[c(p+1, p+2)])))*n

      }
      if(d==1){
        psi_g_combined <- func(beta_j, full.coefficients[(j-1+(k-1)*J)*d+length(beta)+1], 
                               block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
        
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))*n
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[p+1])))*n

      }
    }
  }
  V_psi <- t(psi_g)%*%psi_g/length(subject_indicator)
  W <- solve(V_psi)
  rank <- K*J*p
  Q <- colMeans(psi_g) %*% W %*% colMeans(psi_g)*length(subject_indicator)
  return(list(Q=Q,rank=rank))
}

print.objective.dimm <- function(x, ...){
  res <- cbind(Q = x$Q,
               df = x$df,
               p.value = x$p.value )
  cat("Call:\n")
  print(x$call)
  cat("Testing fit of model:\n")
  print(res)
}

CS.CL.estimation <- function(par, block_y, block_x, id, m, n, div, func, funcderiv){
  ## This is the function to be minimized over the set of parameters within each block
  -func(beta=par[1:dim(block_x)[2]], cov=par[dim(block_x)[2]+1]^2*(matrix(par[dim(block_x)[2]+2],m,m)-diag(par[dim(block_x)[2]+2]-1,m,m)), 
        block_y, block_x, id, m, n)/div
}

CS.CL.estimation.deriv <- function(par, block_y, block_x, id, m, n, div, func, funcderiv){
  ## This is the derivative of CS.CL.estimation
  -Reduce("+", funcderiv(beta=par[1:dim(block_x)[2]], sigma=par[dim(block_x)[2]+1], rho=par[dim(block_x)[2]+2], block_y, block_x, id, m, n))/(div)
}

AR1.CL.estimation <- function(par, block_y, block_x, id, m, n, div, func, funcderiv){
  ## This is the function to be minimized over the set of parameters within each block
  -func(beta=par[1:(dim(block_x)[2])], cov=par[dim(block_x)[2]+1]^2*par[dim(block_x)[2]+2]^abs(outer(1:m, 1:m , "-")), 
        block_y, block_x, id, m, n)/div
}

AR1.CL.estimation.deriv <- function(par, block_y, block_x, id, m, n, div, func, funcderiv){
  ## This is the derivative of AR1.CL.estimation
  -Reduce("+", funcderiv(beta=par[1:dim(block_x)[2]], sigma=par[dim(block_x)[2]+1], rho=par[dim(block_x)[2]+2], block_y, block_x, id, m, n))/(div)
}

independence.CL.estimation <- function(par, block_y, block_x, id, m, n, div, func, funcderiv){
  ## This is the function to be minimized over the set of parameters within each block
  -func(beta=par[1:(dim(block_x)[2])], cov=diag(par[(dim(block_x)[2]+1)]^2,m), block_y, block_x, id, m, n)/div
}

independence.CL.estimation.deriv <- function(par, block_y, block_x, id, m, n, div, func, funcderiv){
  ## This is the function to be minimized over the set of parameters within each block
  -Reduce("+", funcderiv(beta=par[1:dim(block_x)[2]], sigma=par[dim(block_x)[2]+1], block_y, block_x, id, m, n))/(div)
}

psi.g.mean <- function(par, d, response, covariates, J, K, N, func){
  p <- dim(covariates)[2]-1
  psi_g <- matrix(0, N, J*K*(p+d))
  beta <- par[1:p]
  n_k <- vector("numeric", K)
  for (k in 1:K){
    for (j in 1:J){
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.data.frame(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
      id <- block_x[,which(colnames(block_x)=="id")]
      block_x <- as.matrix(block_x[, -which(colnames(block_x)=="id")])
      
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      n_k[k] <- n
      if(d==2){
        psi_g_combined <- func(beta, par[(j-1+(k-1)*J)*d+p+1], par[(j-1+(k-1)*J)*d+p+2], 
                                block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
        
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[c(p+1, p+2)])))
      }
      if(d==1){
        psi_g_combined <- func(beta, par[(j-1+(k-1)*J)*d+p+1], 
                               block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
        
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
          t(sapply(psi_g_combined, function(x) as.matrix(x[p+1])))
      }
    }
  }
  return(colMeans(psi_g))
}

psi.g.deriv.mean <- function(par, d, response, covariates, J, K, N, funcderiv){
  p <- dim(covariates)[2]-1
  S <- matrix(0, J*K*d+p, J*K*(p+d))
  beta <- par[1:p]
  for (k in 1:K){
    for (j in 1:J){
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.matrix(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
      id <- block_x[,which(colnames(block_x)=="id")]
      block_x <- as.matrix(block_x[,  -which(colnames(block_x)=="id")])
      
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      if(d==2){
        little_s <- t(funcderiv(beta, par[(j-1+(k-1)*J)*d+p+1], par[(j-1+(k-1)*J)*d+p+2], block_y, block_x, id, m, n)*n/N )
      }
      if(d==1){
        little_s <- t(funcderiv(beta, par[(j-1+(k-1)*J)*d+p+1], block_y, block_x, id, m, n)*n/N)
      }
      S[1:p,((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- little_s[1:p,1:p]
      S[1:p,(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- little_s[1:p,(p+1):(p+d)]
      S[((j-1+(k-1)*J)*d+p+1):((j-1+(k-1)*J)*d+p+d),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- little_s[(p+1):(p+d),1:p]
      S[((j-1+(k-1)*J)*d+p+1):((j-1+(k-1)*J)*d+p+d),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- little_s[(p+1):(p+d),(p+1):(p+d)]
    }
  }
  return(S)
}

combined.estimation.mean <- function(par, d, W, response, covariates, J, K, N, div, func, funcderiv){
  psi_g_mean <- psi.g.mean(par, d, response, covariates, J, K, N, func)
  return(N * psi_g_mean %*% W %*% psi_g_mean/div)
}

estimation.deriv.mean <- function(par, d, W, response, covariates, J, K, N, div, func, funcderiv){
  2 * N * psi.g.mean(par, d, response, covariates, J, K, N, func) %*% W %*% 
    t(psi.g.deriv.mean(par, d, response, covariates, J, K, N, funcderiv))/div
}

estimation.deriv2.mean <- function(par, d, W, response, covariates, J, K, N, div, func, funcderiv){
  2 * N * psi.g.deriv.mean(par, d, response, covariates, J, K, N, funcderiv) %*% W %*%
    t(psi.g.deriv.mean(par, d, response, covariates, J, K, N, funcderiv))/div
}

ridge.estimate.mean <- function(psi_list, MCLE, MCLE_mean, response, covariates, J, K, N, funcderiv, folds, lam){
  S <- list()
  W <- list()
  lam_list <- vector("numeric", K)
  scale <- matrix(0, dim(covariates)[2]-1, dim(covariates)[2]-1)
  main <- matrix(0, dim(covariates)[2]-1, 1)
  for (k in 1:K){
    S[[k]] <- list()
    for (j in 1:J){
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[colnames(response)==k],]
      block_x <- as.matrix(ids_block_x[rep(response[,which(colnames(response)=="response_indicator")], 
                                           length(unique(ids_block_x[,"id"])))==j,
                                       -which(colnames(ids_block_x)=="id")])
      
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      S[[k]][[j]] <- Reduce("+", funcderiv(matrix(MCLE[[k]][[j]],m,m), 
                                           block_y=block_y, block_x=block_x, 
                                           m=m, n=n))/n
    }
    S_k <- do.call(cbind, S[[k]])
    if (length(lam) != 1){
      CV <- vector("numeric", length(lam))
      for(f in 1:length(lam))
        CV[f] <- modified.cholesky.cv(psi_list[[k]], J, dim(covariates)[2]-1, N, lam[f], folds)
      lam_list[k] <- lam[which.min(CV)] 
    } else {
      lam_list[k] <- lam
    }
    W[[k]] <- modified.cholesky(psi_list[[k]], J, dim(covariates)[2]-1, N, lam_list[k])
    for(j in 1:J){
      scale <- scale + n*S_k %*% W[[k]][,((j-1)*(dim(covariates)[2]-1)+1):(j*(dim(covariates)[2]-1))] %*% S[[k]][[j]]
      main <- main + n*S_k %*% W[[k]][,((j-1)*(dim(covariates)[2]-1)+1):(j*(dim(covariates)[2]-1))] %*% S[[k]][[j]] %*%
        MCLE_mean[[k]][[j]]
    }
  }
  return(list(coefficients=solve(scale)%*%main, 
              W=as.matrix(Matrix::bdiag(lapply(W, function(x) matrix(unlist(x), ncol=J*(dim(covariates)[2]-1), byrow=TRUE)))), 
              lam=lam_list))
}

modified.cholesky <- function(psi_list, J, p, N, lam){
  psi_mat <- c()
  D <- matrix(0, p*J, p*J)
  Gamma <- matrix(0, p*J, p*J)
  diag(Gamma) <- 1
  for(j in 1:J){
    psi_mat <- rbind(psi_mat, matrix(unlist(psi_list[[j]]), p, N))
  }
  psi_mat <- t(psi_mat)
  for(r in 2:dim(psi_mat)[2]){
    ridge <- MASS::lm.ridge(psi_mat[,r] ~ 0 + psi_mat[,1:(r-1)], lambda=lam)
    Gamma[r, 1:(r-1)] <- -as.matrix(ridge$coef)
    D[r,r] <- var(psi_mat[,r] - psi_mat[,1:(r-1)]%*%as.matrix(ridge$coef))
  }
  D[1,1] <- var(psi_mat[,1])
  return((t(Gamma)%*%solve(D)%*%Gamma))
}

modified.cholesky.cv <- function(psi_list, J, p, N, lam, folds){
  psi_mat <- c()
  for(j in 1:J){
    psi_mat <- rbind(psi_mat, matrix(unlist(psi_list[[j]]), p, N))
  }
  psi_mat <- t(psi_mat)
  s_v <- c()
  CV <- vector("numeric", folds)
  for(f in 1:folds){
    set.seed(f)
    indices <- sample(setdiff(seq(1:dim(psi_mat)[1]), s_v), floor(dim(psi_mat)[1]/folds), replace=FALSE)
    s_v <- c(s_v, indices)
    psi_mat_sub <- psi_mat[setdiff(seq(1:dim(psi_mat)[1]), indices),]
    D <- matrix(0, p*J, p*J)
    Gamma <- matrix(0, p*J, p*J)
    diag(Gamma) <- 1
    for(r in 2:dim(psi_mat_sub)[2]){
      ridge <- MASS::lm.ridge(psi_mat_sub[,r] ~ 0 + psi_mat_sub[,1:(r-1)], lambda=lam)
      Gamma[r, 1:(r-1)] <- -as.matrix(ridge$coef)
      D[r,r] <- var(psi_mat_sub[,r] - psi_mat_sub[,1:(r-1)]%*%as.matrix(ridge$coef))
    }
    D[1,1] <- var(psi_mat_sub[,1])
    W <- (t(Gamma)%*%solve(D)%*%Gamma)
    V <- cov(psi_mat_sub)
    CV[f] <- length(indices)*log(det(V)) + sum(apply(psi_mat[indices,],1,function(x) t(x) %*% W %*% x))
  }
  return(sum(CV, na.rm=T)/folds)
}

dimm <- function(formula, data, id=id, response_indicator=NULL, subject_indicator=NULL, family, corstr, analysis=NULL, cluster=NULL, comb_scheme=NULL, ...){
  cl <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mf[,"(id)"] <- data[,id]
  mt <- attr(mf, "terms")
  response <- data.frame(id=mf[,"(id)"], model.response(mf, "numeric"))
  covariates <- data.frame(id=mf[,"(id)"], model.matrix(mt, mf))
  colnames(covariates)[match("X.Intercept.",colnames(covariates))] <- "(Intercept)"

  if (is.null(response_indicator)){
    warning(gettextf("Response indicator is null. Using one response block"), domain = NA)
    response_indicator <- as.vector(rep(1, table(response[,"id"])[1]))
  }
  if (is.null(subject_indicator)){
    warning(gettextf("Subject indicator is null. Using one subject block"), domain=NA)
    subject_indicator <- as.vector(rep(1, length(unique(response[,"id"]))))
  }
  if (is.empty.model(mt)) stop("Model is empty")
  if (family != "gaussian" & family!="binomial") {
    warning(gettextf("family = '%s' is not supported. Using gaussian", family), domain = NA)
    family <- "gaussian"
  }
  if (corstr != "AR1" & corstr != "CS" & corstr != "independence" & corstr != "CS+AR1") stop("corstr must be supported type")
  method <- "exact"
  if (is.null(analysis)) stop("Block analysis must be supported type")
  if (analysis != "PCL" & analysis != "GEE" & analysis != "QIF") stop("Block analysis must be supported type")
  if (analysis=="PCL" & 1 %in% table(response_indicator))
    stop(paste("Minimum of two repeat measurements required for pairwise composite likelihood in blocks",
               paste(which(table(response_indicator)==1), collapse=", ")))
  if(is.null(comb_scheme)) comb_scheme <- rep(1,length(unique(response_indicator))*length(unique(subject_indicator)))
  if(length(comb_scheme)!=(length(unique(response_indicator))*length(unique(subject_indicator)))) stop("Combination scheme vector is not of correct length. Please specify vector of length JK.")
  
  response <- cbind(matrix(response[,-which(colnames(response) == "id")], length(response_indicator), length(subject_indicator)), response_indicator)
  
  if (is.null(cluster)){
    output <- dimm.compute.mean(response, covariates, response_indicator, subject_indicator, family, corstr, method, analysis, comb_scheme)
  }
  if (!is.null(cluster)){
    if (cluster == 1){
      output <- dimm.compute.mean(response, covariates, response_indicator, subject_indicator, family, corstr, method, analysis, comb_scheme)
    } else {
      output <- dimm.compute.mean.parallel(response, covariates, response_indicator, subject_indicator, 
                                           family, corstr, method, analysis, cluster, comb_scheme)
      }
  }
  if(length(unique(comb_scheme))==1) {
    names(output$coefficients) <- colnames(covariates)[-1]
  } else {
    names(output$coefficients) <- c(t(sapply(colnames(covariates)[-1], function(x) paste(x, unique(comb_scheme), sep="-"))))
  }
  colnames(output$vcov) <- names(output$coefficients)
  rownames(output$vcov) <- names(output$coefficients)
  output$response <- data.frame(id=mf[,"(id)"], model.response(mf, "numeric"))
  output$covariates <- data.frame(id=mf[,"(id)"], model.matrix(mt, mf))
  colnames(output$response)[1] <- id
  output$id <- id
  output$comb_scheme <- comb_scheme
  
  output <- c(output, list(call=cl, formula=formula))
  class(output) <- "dimm"
  return(output)
}

dimm.compute.mean <- function(response, covariates, response_indicator, subject_indicator, family, corstr, method, analysis, comb_scheme, ...){
  time1 <- proc.time()
  output <- list()
  
  colnames(response) <- c(subject_indicator,"response_indicator")
  
  J <- length(unique(response_indicator))
  K <- length(unique(subject_indicator))
  N <- length(subject_indicator)
  p <- dim(covariates)[2]-1
  
  if (analysis=="PCL" & family=="gaussian"){
    func1 <- function(beta, cov, block_y, block_x, id, m, n){logCLnormal_par(beta, cov, block_y, block_x, m, n)}
    if (corstr=="CS"){
      func4 <- function(beta, sigma, rho, block_y, block_x, id, m, n){eenormalCSvar_par(beta, sigma, rho, block_y, block_x, m, n)}
      func5 <- function(beta, sigma, rho, block_y, block_x, id, m, n){
        -Reduce("+", psi_g_CSderivvar_par(beta, sigma, rho, block_y, block_x, m, n))
      }
      d <- 2
      lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
      upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
    }
    if (corstr=="AR1"){
      func4 <- function(beta, sigma, rho, block_y, block_x, id, m, n){eenormalAR1var_par(beta, sigma, rho, block_y, block_x, m, n)}
      func5 <- function(beta, sigma, rho, block_y, block_x, id, m, n){
        -Reduce("+", psi_g_AR1derivvar_par(beta, sigma, rho, block_y, block_x, m, n))
      }
      d <- 2
      lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
      upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
    }
    if (corstr=="independence"){
      func4 <- function(beta, sigma, block_y, block_x, id, m, n){eenormalindvar_par(beta, sigma, block_y, block_x, m, n)}
      func5 <- function(beta, sigma, block_y, block_x, id, m, n){
        -Reduce("+", psi_g_indderivvar_par(beta, sigma, block_y, block_x, m, n))
      }
      d <- 1
      lower_lim <- c(rep(-Inf,p), 1e-5)
      upper_lim <- rep(Inf,p+1)
    }
  }
  if (analysis=="GEE"){
    if(family=="gaussian") {
      family2 <- gaussian()
      if (corstr=="CS"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("exch", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="AR1"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("ar1", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="independence"){
        d <- 1
        lower_lim <- c(rep(-Inf,p), 1e-5)
        upper_lim <- rep(Inf,p+1)
        func45 <- function(beta, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("independence", CORSTRS, -1)
          alpha <- vector("numeric",0)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+1,p+1)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
    }
    if(family=="binomial") {
      family2 <- binomial()
      if (corstr=="CS"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("exch", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="AR1"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("ar1", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="independence"){
        d <- 1
        lower_lim <- c(rep(-Inf,p), 1e-5)
        upper_lim <- rep(Inf,p+1)
        func45 <- function(beta, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("independence", CORSTRS, -1)
          alpha <- vector("numeric",0)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+1,p+1)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat <- -Hessmat/n
        }
      } 
    }
  }
  if (analysis=="QIF"){
    if (corstr=="CS" | corstr=="AR1") d <- 2
    if (corstr=="independence") d <- 1
    if (corstr=="CS+AR1") d <- 3
  }
  
  S <- list()
  if(analysis != "QIF") {
    psi_g <- matrix(0, N, J*K*(p+d))
    MCLE <- list()
    p_1 <- p+d
  } else {
    psi_g <- matrix(0, N, J*K*(p*d))
    p_1 <- p*d
  }
  MCLE_mean <- list()
  n_k <- vector("numeric", K)
  print("Computing block coefficients.", quote=FALSE)
  time1 <- proc.time()-time1
  time_jk <- list()
  for (k in 1:K) {
    S[[k]] <- list()
    if (analysis != "QIF") MCLE[[k]] <- list()
    MCLE_mean[[k]] <- list()
    time_jk[[k]] <- list()
    for (j in 1:J) {
      time_j <- proc.time()
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.data.frame(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
      id <- block_x[,which(colnames(block_x)=="id")]
      block_x <- as.matrix(block_x[, -which(colnames(block_x)=="id")])
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      init_betas <- coef(glm(c(block_y) ~ 0 + block_x, family=family))
      
      if (corstr == "CS") {
        if (analysis=="PCL"){
          div <- choose(m,2)*n
          init_cov_parameters <- c(1,0.5)
          optimization <- optim(par=c(init_betas, init_cov_parameters), fn=CS.CL.estimation, gr=CS.CL.estimation.deriv, block_y=block_y, 
                                block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                lower=lower_lim, upper=upper_lim, control=list(maxit=500))
          init_betas <- optimization$par[1:p]
          init_cov_parameters <- optimization$par[(p+1):(p+2)]
          optimization_2 <- optim(par=c(init_betas, init_cov_parameters), fn=CS.CL.estimation, gr=CS.CL.estimation.deriv, block_y=block_y, 
                                  block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                  lower=lower_lim, upper=upper_lim)
          MCLE_jk <- optimization_2$par
          MCLE[[k]][[j]] <- MCLE_jk[p+1]^2*(matrix(MCLE_jk[p+2],m,m)-diag(MCLE_jk[p+2]-1,m,m))
          
          psi_g_combined <- func4(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], 
                                  block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(psi_g_combined, function(x) as.matrix(x[c(p+1, p+2)])))
          
          MCLE_mean[[k]][[j]] <- MCLE_jk
          S[[k]][[j]] <- func5(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], block_y, block_x, id, m, n)
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
        if (analysis=="GEE"){
          gee_block_fit <- geepack::geese.fit(x=block_x, y=c(block_y), id=id, family=family2, corstr = "exch", scale.fix=FALSE)
          MCLE_jk <- as.vector(c(gee_block_fit$beta, gee_block_fit$gamma, gee_block_fit$alpha))
          MCLE_mean[[k]][[j]] <- MCLE_jk
          MCLE[[k]][[j]] <- (matrix(gee_block_fit$alpha,m,m)-diag(gee_block_fit$alpha-1,m,m))*gee_block_fit$gamma
          
          HessGrad <- func45(gee_block_fit$beta, gee_block_fit$alpha, gee_block_fit$gamma, block_y, block_x, id, m, n)
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- do.call(rbind, HessGrad[6+seq(1,n*3,3)])*n
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- t(matrix(do.call(rbind, HessGrad[6+setdiff(1:(n*3),seq(1,n*3,3))]), 2,n))*n
          
          S[[k]][[j]] <- matrix(0,p+d,p+d)
          S[[k]][[j]][1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          S[[k]][[j]][1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          S[[k]][[j]][1:p,(p+d)] <- t(matrix(HessGrad[[4]],1,p))
          S[[k]][[j]][(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          S[[k]][[j]][(p+1),(p+d)] <- matrix(HessGrad[[5]],1,1)
          S[[k]][[j]][(p+d),(p+d)] <- matrix(HessGrad[[6]],1,1)
          S[[k]][[j]] <- -S[[k]][[j]]*n
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
        if (analysis=="QIF"){
          qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="exchangeable", init_betas, maxit=5000, tol=1e-6)
          MCLE_jk <- as.vector(qif_block_fit$beta_sub)
          MCLE_mean[[k]][[j]] <- MCLE_jk
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- t(do.call(cbind, qif_block_fit$gi_list))
          
          S[[k]][[j]] <- qif_block_fit$G_sub
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
      }
      if (corstr == "AR1") {
        if (analysis=="PCL"){
          div <- choose(m,2)*n
          init_cov_parameters <- c(1,0.5)
          optimization <- optim(par=c(init_betas, init_cov_parameters), fn=AR1.CL.estimation, gr=AR1.CL.estimation.deriv, block_y=block_y, 
                                block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                lower=lower_lim, upper=upper_lim, control=list(maxit=500))
          init_betas <- optimization$par[1:p]
          init_cov_parameters <- optimization$par[(p+1):(p+2)]
          optimization_2 <- optim(par=c(init_betas, init_cov_parameters), fn=AR1.CL.estimation, gr=AR1.CL.estimation.deriv, block_y=block_y, 
                                  block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                  lower=lower_lim, upper=upper_lim)
          MCLE_jk <- optimization_2$par
          MCLE[[k]][[j]] <- MCLE_jk[p+1]^2*MCLE_jk[p+2]^abs(outer(1:m, 1:m , "-"))
          
          psi_g_combined <- func4(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], 
                                  block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(psi_g_combined, function(x) as.matrix(x[c(p+1, p+2)])))
          
          MCLE_mean[[k]][[j]] <- MCLE_jk
          S[[k]][[j]] <- func5(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], block_y, block_x, id, m, n)
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
        if (analysis=="GEE"){
          gee_block_fit <- geepack::geese.fit(x=block_x, y=c(block_y), id=id, family=family2, corstr = "ar1", scale.fix=FALSE)
          MCLE_jk <- as.vector(c(gee_block_fit$beta, gee_block_fit$gamma, gee_block_fit$alpha))
          MCLE_mean[[k]][[j]] <- MCLE_jk
          MCLE[[k]][[j]] <- (gee_block_fit$alpha^abs(outer(1:m, 1:m , "-")))*gee_block_fit$gamma
          
          HessGrad <- func45(gee_block_fit$beta, gee_block_fit$alpha, gee_block_fit$gamma, block_y, block_x, id, m, n)
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- do.call(rbind, HessGrad[6+seq(1,n*3,3)])*n
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- t(matrix(do.call(rbind, HessGrad[6+setdiff(1:(n*3),seq(1,n*3,3))]), 2,n))*n
          
          S[[k]][[j]] <- matrix(0,p+d,p+d)
          S[[k]][[j]][1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          S[[k]][[j]][(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          S[[k]][[j]][(p+d),1:p] <- matrix(HessGrad[[4]],1,p)
          S[[k]][[j]][(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          S[[k]][[j]][(p+d),(p+1)] <- matrix(HessGrad[[5]],1,1)
          S[[k]][[j]][(p+d),(p+d)] <- matrix(HessGrad[[6]],1,1)
          S[[k]][[j]] <- -S[[k]][[j]]*n
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
        if (analysis=="QIF"){
          qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="AR-1", init_betas, maxit=5000, tol=1e-6)
          MCLE_jk <- as.vector(qif_block_fit$beta_sub)
          MCLE_mean[[k]][[j]] <- MCLE_jk
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- t(do.call(cbind, qif_block_fit$gi_list))
          
          S[[k]][[j]] <- qif_block_fit$G_sub
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
      }
      if (corstr == "independence"){
        if (analysis=="PCL"){
          div <- choose(m,2)*n
          init_cov_parameters <- 1
          optimization <- optim(par=c(init_betas, init_cov_parameters), fn=independence.CL.estimation, gr=independence.CL.estimation.deriv, block_y=block_y, 
                                block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                lower=lower_lim, upper=upper_lim, control=list(maxit=500))
          init_betas <- optimization$par[1:p]
          init_cov_parameters <- optimization$par[(p+1)]
          optimization_2 <- optim(par=c(init_betas, init_cov_parameters), fn=independence.CL.estimation, gr=independence.CL.estimation.deriv, block_y=block_y, 
                                  block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                  lower=lower_lim, upper=upper_lim)
          MCLE_jk <- optimization_2$par
          MCLE[[k]][[j]] <- diag(MCLE_jk[p+1]^2, m)
          
          psi_g_combined <- func4(MCLE_jk[1:p], MCLE_jk[p+1], 
                                  block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(psi_g_combined, function(x) as.matrix(x[1:p])))
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(psi_g_combined, function(x) as.matrix(x[p+1])))
          
          MCLE_mean[[k]][[j]] <- MCLE_jk
          S[[k]][[j]] <- func5(MCLE_jk[1:p], MCLE_jk[p+1], block_y, block_x, id, m, n)
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
        if (analysis=="GEE"){
          gee_block_fit <- geepack::geese.fit(x=block_x, y=c(block_y), id=id, family=family2, corstr = "independence", scale.fix=FALSE)
          MCLE_jk <- as.vector(c(gee_block_fit$beta, gee_block_fit$gamma))
          MCLE_mean[[k]][[j]] <- MCLE_jk
          MCLE[[k]][[j]] <- diag(gee_block_fit$gamma,m,m)
          
          HessGrad <- func45(gee_block_fit$beta, gee_block_fit$gamma, block_y, block_x, id, m, n)
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- do.call(rbind, HessGrad[6+seq(1,n*3,3)])*n
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- do.call(rbind, HessGrad[6+1+seq(1,n*3,3)])*n
          
          S[[k]][[j]] <- matrix(0,p+d,p+d)
          S[[k]][[j]][1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          S[[k]][[j]][(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          S[[k]][[j]][(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          S[[k]][[j]] <- -S[[k]][[j]]*n
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
        if (analysis=="QIF"){
          qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="independence", init_betas, maxit=5000, tol=1e-6)
          MCLE_jk <- as.vector(qif_block_fit$beta_sub)
          MCLE_mean[[k]][[j]] <- MCLE_jk
          
          psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- t(do.call(cbind, qif_block_fit$gi_list))
          
          S[[k]][[j]] <- qif_block_fit$G_sub
          time_jk[[k]][[j]] <- proc.time()-time_j
        }
      }
      if (corstr == "CS+AR1"){
        qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="CS+AR1", init_betas, maxit=5000, tol=1e-6)
        MCLE_jk <- as.vector(qif_block_fit$beta_sub)
        MCLE_mean[[k]][[j]] <- MCLE_jk
        
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
              ((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- t(do.call(cbind, qif_block_fit$gi_list))
        
        S[[k]][[j]] <- qif_block_fit$G_sub
        time_jk[[k]][[j]] <- proc.time()-time_j
      }
      
      n_k[k] <- n
    }
    time_jk[[k]] <- do.call(rbind, time_jk[[k]])
  }
  time_list <- do.call(rbind, time_jk)
  time_max <- time_list[which.max(time_list[,1]),]
  time2 <- proc.time()
  
  unordered_psi_g <- psi_g
  unordered_S <- S
  
  print("Computing combined estimate.", quote=FALSE)
  
  if (method=="exact" & analysis != "QIF"){
    psi_g <- psi_g[,c(order(sapply(comb_scheme, function(x) rep(x,p))),(J*K*p+1):(J*K*p_1))]
    S_all <- do.call(rbind, lapply(unordered_S, function(S_k) do.call(rbind, S_k)))[c(c(matrix(1:(J*K*p_1), nrow=p_1)[1:p,]),
                                                                                          c(matrix(1:(J*K*p_1), nrow=p_1)[(p+1):p_1,])),]
    S_all_reordered <- S_all[c(order(sapply(comb_scheme, function(x) rep(x,p))),(J*K*p+1):(J*K*p_1)),]
    
    S_matrix <- as.matrix(Matrix::bdiag(Matrix::bdiag(
      lapply(sort(unique(comb_scheme)), function(c) {
        S_all_reordered[c(sapply(comb_scheme, function(x) rep(x,p))[order(sapply(comb_scheme, function(x) rep(x,p)))]==c,
                           rep(FALSE, J*K*d)),1:p, drop=FALSE]
      })),
      Matrix::bdiag(
      lapply(1:(J*K), function(c){
        S_all_reordered[(J*K*p+(c-1)*d+1):(J*K*p+c*d),(p+1):p_1, drop=FALSE]
        })
      )))
    
    S_MCLE_mean_list <- mapply(function(S_k,MCLE_k) mapply(function(x,y) x %*% y, S_k, MCLE_k, SIMPLIFY=FALSE), S, MCLE_mean, SIMPLIFY = FALSE)
    S_MCLE_mean_all <- do.call(rbind, lapply(S_MCLE_mean_list, function(S_k) do.call(rbind, S_k)))[c(c(matrix(1:(J*K*p_1), nrow=p_1)[1:p,]),
                                                                                                         c(matrix(1:(J*K*p_1), nrow=p_1)[(p+1):p_1,]))]
    S_MCLE_mean_matrix <- S_MCLE_mean_all[c(order(sapply(comb_scheme, function(x) rep(x,p))),(J*K*p+1):(J*K*p_1))]
    
    V_psi <- t(psi_g)%*%psi_g/N
    output$W <- solve(V_psi)
    
    scale <- t(S_matrix) %*% output$W %*% S_matrix
    main <- t(S_matrix) %*% output$W %*% matrix(S_MCLE_mean_matrix, ncol=1)
    
    output$coefficients <- as.vector(solve(scale)[1:(p*length(unique(comb_scheme))),]%*%main)
    output$vcov <- solve(scale)[1:(p*length(unique(comb_scheme))),1:(p*length(unique(comb_scheme)))]*N
    output$full.coefficients <- as.vector(solve(scale)%*%main)
  } 
  
  if (method=="exact" & analysis == "QIF"){
    psi_g <- unordered_psi_g[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
    S_all_reordered <- t(do.call(rbind, lapply(unordered_S, function(S_k) do.call(rbind, S_k))))[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
    
    rank <- qr(psi_g)$rank
    if(dim(psi_g)[2] != rank){
      psi_g_pca <- prcomp(psi_g, center = FALSE,scale. = FALSE)
      old_psi_g <- psi_g
      old_J <- J
      old_K <- K
      old_S <- S
      ta <- cumsum(table(comb_scheme))
      psi_g <- psi_g_pca$x[,1:rank]
      S_matrix <- (as.matrix(Matrix::bdiag(lapply(sort(unique(comb_scheme)), function(c) {
        S_all_reordered[,sapply(comb_scheme, function(x) rep(x,p_1))[order(sapply(comb_scheme, function(x) rep(x,p_1)))]==c]
      }))) %*% psi_g_pca$rotation )[,1:rank]
      S_MCLE_mean_list <- mapply(function(S_k,MCLE_k) mapply(function(x,y) x %*% y, S_k, MCLE_k, SIMPLIFY=FALSE), S, MCLE_mean, SIMPLIFY = FALSE)
      S_MCLE_mean_matrix <- (t(do.call(rbind, lapply(S_MCLE_mean_list, function(S_k) 
        do.call(rbind, S_k))))[,order(sapply(comb_scheme, function(x) rep(x,p_1)))] %*% 
          psi_g_pca$rotation)[,1:rank]
      
    } else {
      S_matrix <- as.matrix(Matrix::bdiag(lapply(sort(unique(comb_scheme)), function(c) {
        S_all_reordered[,sapply(comb_scheme, function(x) rep(x,p_1))[order(sapply(comb_scheme, function(x) rep(x,p_1)))]==c]
      })))
      S_MCLE_mean_list <- mapply(function(S_k,MCLE_k) mapply(function(x,y) x %*% y, S_k, MCLE_k, SIMPLIFY=FALSE), S, MCLE_mean, SIMPLIFY = FALSE)
      S_MCLE_mean_matrix <- t(do.call(rbind, lapply(S_MCLE_mean_list, function(S_k) 
        do.call(rbind, S_k))))[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
    }
    
    V_psi <- t(psi_g)%*%psi_g/N
    output$W <- solve(V_psi)
    
    scale <- S_matrix %*% output$W %*% t(S_matrix)
    main <- S_matrix %*% output$W %*% matrix(S_MCLE_mean_matrix, ncol=1)
    output$coefficients <- as.vector(solve(scale)%*%main)
    output$vcov <- solve(scale)*N
  } 
  
  if (method=="iterative" | method=="ridge"){
    if(method=="iterative"){
      div <- choose(N,2)
      init_pars <- c(colMeans(do.call(rbind, lapply(MCLE_mean, function(x) do.call(rbind, lapply(x, function(y) y[1:p]))))), 
                     unlist(lapply(MCLE_mean, function(x) lapply(x, function(y) y[(p+1):(p+d)]))))
      optimization <- optim(par=init_pars, 
                            combined.estimation.mean, gr=estimation.deriv.mean, 
                            d=d, W=output$W, response=response, 
                            covariates=covariates, J=J, K=K, N=N, div=div, func=func4, funcderiv=func5, 
                            method="L-BFGS-B", lower=lower_lim, upper=upper_lim, control=list(maxit=500))
      init_pars <- optimization$par
      output$coefficients <- optim(par=init_pars, 
                                   combined.estimation.mean, gr=estimation.deriv.mean, 
                                   d=d, W=output$W, response=response, 
                                   covariates=covariates, J=J, K=K, N=N, div=div, func=func4, funcderiv=func5, 
                                   method="L-BFGS-B", lower=lower_lim, upper=upper_lim)$par
    }
    
    sensitivity_psi <- matrix(0, J*K*d+p, J*K*(p+d))
    phi_gh <- matrix(0, N, J*K*(p+d))
    for(k in 1:K) {
      for (j in 1:J) {
        block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                          nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
        ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
        block_x <- as.matrix(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,
                                         -which(colnames(ids_block_x)=="id")])
        m <- dim(block_y)[1]
        n <- dim(block_y)[2]
        if(d==2){
          phi_gh_combined <- func4(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], output$coefficients[p+(j-1+(k-1)*J)*d+2], 
                                  block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[1:p])))
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[c(p+1, p+2)])))
          
          little_s <- t(func5(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], 
                              output$coefficients[p+(j-1+(k-1)*J)*d+2], block_y, block_x, id, m, n)*n/N)
        }
        if(d==1){
          phi_gh_combined <- func4(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], 
                                   block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                 ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[1:p])))
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                 (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[p+1])))
          
          little_s <- t(func5(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], 
                                           block_y, block_x, id, m, n)*n/N)
        }
        sensitivity_psi[1:p,((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- little_s[1:p,1:p]
        sensitivity_psi[1:p,(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- little_s[1:p,(p+1):(p+d)]
        sensitivity_psi[((j-1+(k-1)*J)*d+p+1):((j-1+(k-1)*J)*d+p+d),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- little_s[(p+1):(p+d),1:p]
        sensitivity_psi[((j-1+(k-1)*J)*d+p+1):((j-1+(k-1)*J)*d+p+d),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- little_s[(p+1):(p+d),(p+1):(p+d)]
      }
    }
    variability_psi <- matrix(0, J*K*(p+d), J*K*(p+d))
    for(i in 1:dim(phi_gh)[1]){
      variability_psi <- variability_psi + phi_gh[i,]%o%phi_gh[i,]
    }
    variability_psi <- variability_psi/N
    
    output$vcov <- solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))%*%
      sensitivity_psi%*%output$W%*%variability_psi%*%output$W%*%t(sensitivity_psi)%*%
      solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))/N
    
    full_coefficients <- output$coefficients
    output$coefficients <- output$coefficients[1:p]
    output$vcov <- output$vcov[1:p,1:p]
    output$full.coefficients <- full_coefficients
  }
  time2 <- proc.time()-time2
  output$family <- family
  output$analysis <- analysis
  output$corstr <- corstr
  output$response_indicator <- response_indicator
  output$subject_indicator <- subject_indicator
  output$J <- J
  output$K <- K
  output$MCLE <- list()
  output$MCLE$beta <- MCLE_mean
  if(analysis != "QIF") output$MCLE$vcov <- MCLE
  output$time <- time1+time2+time_max
  
  return(output)
}

dimm.compute.mean.parallel <- function(response, covariates, response_indicator, subject_indicator, 
                                        family, corstr, method, analysis, cluster, comb_scheme, ...){
  time1 <- proc.time()
  output <- list()
  
  colnames(response) <- c(subject_indicator,"response_indicator")
  
  J <- length(unique(response_indicator))
  K <- length(unique(subject_indicator))
  N <- length(subject_indicator)
  p <- dim(covariates)[2]-1
  
  if (analysis=="PCL" & family=="gaussian"){
    func1 <- function(beta, cov, block_y, block_x, id, m, n){logCLnormal_par(beta, cov, block_y, block_x, m, n)}
    if (corstr=="CS"){
      func4 <- function(beta, sigma, rho, block_y, block_x, id, m, n){eenormalCSvar_par(beta, sigma, rho, block_y, block_x, m, n)}
      func5 <- function(beta, sigma, rho, block_y, block_x, id, m, n){
        -Reduce("+", psi_g_CSderivvar_par(beta, sigma, rho, block_y, block_x, m, n))
      }
      d <- 2
      lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
      upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
    }
    if (corstr=="AR1"){
      func4 <- function(beta, sigma, rho, block_y, block_x, id, m, n){eenormalAR1var_par(beta, sigma, rho, block_y, block_x, m, n)}
      func5 <- function(beta, sigma, rho, block_y, block_x, id, m, n){
        -Reduce("+", psi_g_AR1derivvar_par(beta, sigma, rho, block_y, block_x, m, n))
      }
      d <- 2
      lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
      upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
    }
    if (corstr=="independence"){
      func4 <- function(beta, sigma, block_y, block_x, id, m, n){eenormalindvar_par(beta, sigma, block_y, block_x, m, n)}
      func5 <- function(beta, sigma, block_y, block_x, id, m, n){
        -Reduce("+", psi_g_indderivvar_par(beta, sigma, block_y, block_x, m, n))
      }
      d <- 1
      lower_lim <- c(rep(-Inf,p), 1e-5)
      upper_lim <- rep(Inf,p+1)
    }
  }
  if (analysis=="GEE"){
    if(family=="gaussian") {
      family2 <- gaussian()
      if (corstr=="CS"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("exch", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="AR1"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("ar1", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="independence"){
        d <- 1
        lower_lim <- c(rep(-Inf,p), 1e-5)
        upper_lim <- rep(Inf,p+1)
        func45 <- function(beta, gamma, block_y, block_x, id, m, n) {
          family2 <- gaussian()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("independence", CORSTRS, -1)
          alpha <- vector("numeric",0)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+1,p+1)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
    }
    if(family=="binomial") {
      family2 <- binomial()
      if (corstr=="CS"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("exch", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="AR1"){
        d <- 2
        lower_lim <- c(rep(-Inf,p), 1e-5, -1+1e-5)
        upper_lim <- c(rep(Inf,p), Inf, 1-(1e-5))
        func45 <- function(beta, alpha, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("ar1", CORSTRS, -1)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, alpha, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, alpha, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, alpha, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+2,p+2)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+2),1:p] <- matrix(HessGrad[[4]],1,p)
          Hessmat[1:p,(p+2)] <- t(matrix(HessGrad[[4]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat[(p+1),(p+2)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+1)] <- matrix(HessGrad[[5]],1,1)
          Hessmat[(p+2),(p+2)] <- matrix(HessGrad[[6]],1,1)
          Hessmat <- -Hessmat/n
        }
      }
      if (corstr=="independence"){
        d <- 1
        lower_lim <- c(rep(-Inf,p), 1e-5)
        upper_lim <- rep(Inf,p+1)
        func45 <- function(beta, gamma, block_y, block_x, id, m, n) {
          family2 <- binomial()
          LINKS <- c("identity", "logit", "probit", "cloglog", "log", 
                     "inverse", "fisherz", "lwybc2", "lwylog")
          VARIANCES <- c("gaussian", "binomial", "poisson", "Gamma")
          CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", 
                       "userdefined", "fixed")
          clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
          clusz <- c(clusnew[1], diff(clusnew))
          maxclsz <- max(clusz)
          waves <- as.integer(unlist(sapply(clusz, function(x) 1:x)))
          mean.link.v <- pmatch(family2$link, LINKS, -1, TRUE)
          cor.link.v <- sca.link.v <- 1
          variance.v <- pmatch(family2$family, VARIANCES, -1, TRUE)
          corstrv <- pmatch("independence", CORSTRS, -1)
          alpha <- vector("numeric",0)
          
          getHessGrads(y=c(block_y), x=block_x, offset=rep(0,dim(block_x)[1]), doffset=rep(0,dim(block_x)[1]), 
                       w=rep(1,dim(block_x)[1]), linkwave=as.integer(rep(1, dim(block_x)[1])),  
                       zsca=matrix(1, dim(block_x)[1], 1), zcor=geepack::genZcor(clusz, waves, corstrv), corp=as.double(waves), 
                       clusz, geestr = list(length(mean.link.v), as.integer(mean.link.v), 
                                            as.integer(variance.v), 1, 1, as.integer(FALSE)), 
                       cor=list(as.integer(corstrv), maxclsz), par=list(beta, alpha, gamma)) 
        }
        func4 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          Gradmat <- matrix(unlist(HessGrad[-c(1:6)]), length(c(beta, gamma)), n)
          return(lapply(seq_len(n), function(i) Gradmat[,i]))
        }
        func5 <- function(beta, gamma, block_y, block_x, id, m, n){
          HessGrad <- func45(beta, gamma, block_y, block_x, id, m, n)
          p <- length(beta)
          Hessmat <- matrix(0,p+1,p+1)
          Hessmat[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
          Hessmat[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
          Hessmat[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
          Hessmat[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
          Hessmat <- -Hessmat/n
        }
      } 
    }
  }
  if (analysis=="QIF"){
    if (corstr=="CS" | corstr=="AR1") d <- 2
    if (corstr=="independence") d <- 1
    if (corstr=="CS+AR1") d <- 3
  }
  
  S <- list()
  if(analysis != "QIF") {
    psi_g <- matrix(0, N, J*K*(p+d))
    MCLE <- list()
    p_1 <- p+d
  } else {
    psi_g <- matrix(0, N, J*K*(p*d))
    p_1 <- p*d
  }
  MCLE_mean <- list()
  n_k <- vector("numeric", K)
  time_list <- list()
  time1 <- proc.time()-time1
  
  sock <- parallel::makeCluster(rep("localhost", 4), type = "SOCK")
  doParallel::registerDoParallel(sock)
  
  print("Computing block coefficients.", quote=FALSE)
  MCLE_psilist_n <- foreach::foreach(k=1:K) %:% foreach::foreach(j=1:J) %dopar% {
    time_j <- proc.time()
    block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                      nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
    ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
    block_x <- as.data.frame(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
    id <- block_x[,which(colnames(block_x)=="id")]
    block_x <- as.matrix(block_x[,  -which(colnames(block_x)=="id")])
    m <- dim(block_y)[1]
    n <- dim(block_y)[2]
    init_betas <- coef(glm(c(block_y) ~ 0 + block_x, family=family))
    
    if (corstr == "CS") {
      if (analysis=="PCL"){
        div <- choose(m,2)*n
        init_cov_parameters <- c(1,0.5)
        optimization <- optim(par=c(init_betas, init_cov_parameters), fn=CS.CL.estimation, gr=CS.CL.estimation.deriv, block_y=block_y, 
                              block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                              lower=lower_lim, upper=upper_lim, control=list(maxit=500))
        init_betas <- optimization$par[1:p]
        init_cov_parameters <- optimization$par[(p+1):(p+2)]
        optimization_2 <- optim(par=c(init_betas, init_cov_parameters), fn=CS.CL.estimation, gr=CS.CL.estimation.deriv, block_y=block_y, 
                                block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                lower=lower_lim, upper=upper_lim)
        MCLE_jk <- optimization_2$par
        MCLE_matrix <- MCLE_jk[dim(block_x)[2]+1]^2*(matrix(MCLE_jk[dim(block_x)[2]+2],m,m)-diag(MCLE_jk[dim(block_x)[2]+2]-1,m,m))
        
        psi_g_combined <- func4(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], 
                                block_y=block_y, block_x=block_x, id=id, m=m, n=n)
        
        return(list(MCLE_matrix,
                    MCLE_jk,
                    t(sapply(psi_g_combined, function(x) as.matrix(x[1:p]))),
                    t(sapply(psi_g_combined, function(x) as.matrix(x[c(p+1, p+2)]))),
                    func5(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], block_y, block_x, id, m, n),
                    n,
                    proc.time()-time_j)) 
      }
      if (analysis=="GEE"){
        gee_block_fit <- geepack::geese.fit(x=block_x, y=c(block_y), id=id, family=family2, corstr = "exch", scale.fix=FALSE)
        MCLE_jk <- as.vector(c(gee_block_fit$beta, gee_block_fit$gamma, gee_block_fit$alpha))
        MCLE_matrix <- (matrix(gee_block_fit$alpha,m,m)-diag(gee_block_fit$alpha-1,m,m))*gee_block_fit$gamma
        
        HessGrad <- func45(gee_block_fit$beta, gee_block_fit$alpha, gee_block_fit$gamma, block_y, block_x, id, m, n)
        
        Hess <- matrix(0,p+d,p+d)
        Hess[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
        Hess[1:p,(p+1)] <- t(matrix(HessGrad[[2]],1,p))
        Hess[1:p,(p+d)] <- t(matrix(HessGrad[[4]],1,p))
        Hess[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
        Hess[(p+1),(p+d)] <- matrix(HessGrad[[5]],1,1)
        Hess[(p+d),(p+d)] <- matrix(HessGrad[[6]],1,1)
        Hess <- -Hess*n
        return(list(MCLE_matrix,
                    MCLE_jk,
                    do.call(rbind, HessGrad[6+seq(1,n*3,3)])*n,
                    t(matrix(do.call(rbind, HessGrad[6+setdiff(1:(n*3),seq(1,n*3,3))]), 2,n))*n,
                    Hess,
                    n,
                    proc.time()-time_j))
      }
      if (analysis=="QIF"){
        qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="exchangeable", init_betas, maxit=5000, tol=1e-6)
        return(list(NA, as.vector(qif_block_fit$beta_sub), t(do.call(cbind, qif_block_fit$gi_list)), NA, qif_block_fit$G_sub, n, proc.time()-time_j))
      }
    }
    if (corstr == "AR1") {
      if (analysis=="PCL"){
        div <- choose(m,2)*n
        init_cov_parameters <- c(1,0.5)
        optimization <- optim(par=c(init_betas, init_cov_parameters), fn=AR1.CL.estimation, gr=AR1.CL.estimation.deriv, block_y=block_y, 
                              block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                              lower=lower_lim, upper=upper_lim, control=list(maxit=500))
        init_betas <- optimization$par[1:p]
        init_cov_parameters <- optimization$par[(p+1):(p+2)]
        optimization_2 <- optim(par=c(init_betas, init_cov_parameters), fn=AR1.CL.estimation, gr=AR1.CL.estimation.deriv, block_y=block_y, 
                                block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                lower=lower_lim, upper=upper_lim)
        MCLE_jk <- optimization_2$par
        MCLE_matrix <- MCLE_jk[p+1]^2*MCLE_jk[p+2]^abs(outer(1:m, 1:m , "-"))
        psi_g_combined <- func4(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], 
                                block_y=block_y, block_x=block_x, id=id, m=m, n=n)
        
        return(list(MCLE_matrix,
                    MCLE_jk,
                    t(sapply(psi_g_combined, function(x) as.matrix(x[1:p]))),
                    t(sapply(psi_g_combined, function(x) as.matrix(x[c(p+1, p+2)]))),
                    func5(MCLE_jk[1:p], MCLE_jk[p+1], MCLE_jk[p+2], block_y, block_x, id, m, n),
                    n,
                    proc.time()-time_j)) 
      }
      if (analysis=="GEE"){
        gee_block_fit <- geepack::geese.fit(x=block_x, y=c(block_y), id=id, family=family2, corstr = "ar1", scale.fix=FALSE)
        MCLE_jk <- as.vector(c(gee_block_fit$beta, gee_block_fit$gamma, gee_block_fit$alpha))
        MCLE_matrix <- (gee_block_fit$alpha^abs(outer(1:m, 1:m , "-")))*gee_block_fit$gamma
        
        HessGrad <- func45(gee_block_fit$beta, gee_block_fit$alpha, gee_block_fit$gamma, block_y, block_x, id, m, n)
        
        Hess <- matrix(0,p+d,p+d)
        Hess[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
        Hess[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
        Hess[(p+d),1:p] <- matrix(HessGrad[[4]],1,p)
        Hess[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
        Hess[(p+d),(p+1)] <- matrix(HessGrad[[5]],1,1)
        Hess[(p+d),(p+d)] <- matrix(HessGrad[[6]],1,1)
        Hess <- -Hess*n
        return(list(MCLE_matrix,
                    MCLE_jk,
                    do.call(rbind, HessGrad[6+seq(1,n*3,3)])*n,
                    t(matrix(do.call(rbind, HessGrad[6+setdiff(1:(n*3),seq(1,n*3,3))]), 2,n))*n,
                    Hess,
                    n,
                    proc.time()-time_j))
      }
      if (analysis=="QIF"){
        qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="AR-1", init_betas, maxit=5000, tol=1e-6)
        return(list(NA, as.vector(qif_block_fit$beta_sub), t(do.call(cbind, qif_block_fit$gi_list)), NA, qif_block_fit$G_sub, n, proc.time()-time_j))
      }
    }
    if (corstr == "independence") {
      if (analysis=="PCL"){
        div <- choose(m,2)*n
        init_cov_parameters <- 1
        optimization <- optim(par=c(init_betas, init_cov_parameters), fn=independence.CL.estimation, gr=independence.CL.estimation.deriv, block_y=block_y, 
                              block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                              lower=lower_lim, upper=upper_lim, control=list(maxit=500))
        init_betas <- optimization$par[1:p]
        init_cov_parameters <- optimization$par[(p+1)]
        optimization_2 <- optim(par=c(init_betas, init_cov_parameters), fn=independence.CL.estimation, gr=independence.CL.estimation.deriv, block_y=block_y, 
                                block_x=block_x, id=id, m=m, n=n, div=div, func=func1, funcderiv=func4, method="L-BFGS-B",
                                lower=lower_lim, upper=upper_lim)
        MCLE_jk <- optimization_2$par
        MCLE_matrix <- diag(MCLE_jk[p+1]^2, m)
        psi_g_combined <- func4(MCLE_jk[1:p], MCLE_jk[p+1], block_y=block_y, block_x=block_x, m=m, n=n)
        
        return(list(MCLE_matrix,
                    MCLE_jk,
                    t(sapply(psi_g_combined, function(x) as.matrix(x[1:p]))),
                    t(sapply(psi_g_combined, function(x) as.matrix(x[p+1]))),
                    func5(MCLE_jk[1:p], MCLE_jk[p+1], block_y, block_x, id, m, n),
                    n,
                    proc.time()-time_j)) 
      }
      if (analysis=="GEE"){
        gee_block_fit <- geepack::geese.fit(x=block_x, y=c(block_y), id=id, family=family2, corstr = "independence", scale.fix=FALSE)
        MCLE_jk <- as.vector(c(gee_block_fit$beta, gee_block_fit$gamma))
        MCLE_matrix <- diag(gee_block_fit$gamma,m,m)
        
        HessGrad <- func45(gee_block_fit$beta, gee_block_fit$gamma, block_y, block_x, id, m, n)
        
        Hess <- matrix(0,p+d,p+d)
        Hess[1:p,1:p] <- matrix(HessGrad[[1]],p,p)
        Hess[(p+1),1:p] <- matrix(HessGrad[[2]],1,p)
        Hess[(p+1),(p+1)] <- matrix(HessGrad[[3]],1,1)
        Hess <- -Hess*n
        return(list(MCLE_matrix,
                    MCLE_jk,
                    do.call(rbind, HessGrad[6+seq(1,n*3,3)])*n,
                    do.call(rbind, HessGrad[6+1+seq(1,n*3,3)])*n,
                    Hess,
                    n,
                    proc.time()-time_j))
      }
      if (analysis=="QIF"){
        qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="independence", init_betas, maxit=5000, tol=1e-6)
        return(list(NA, as.vector(qif_block_fit$beta_sub), t(do.call(cbind, qif_block_fit$gi_list)), NA, qif_block_fit$G_sub, n, proc.time()-time_j))
      }
    }
    if (corstr == "CS+AR1"){
      qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr="CS+AR1", init_betas, maxit=5000, tol=1e-6)
      return(list(NA, as.vector(qif_block_fit$beta_sub), t(do.call(cbind, qif_block_fit$gi_list)), NA, qif_block_fit$G_sub, n, proc.time()-time_j))
    }
  }
  time2 <- proc.time()
  for(k in 1:K){
    if(analysis != "QIF") MCLE[[k]] <- list()
    MCLE_mean[[k]] <- list()
    S[[k]] <- list()
    n_k[k] <- MCLE_psilist_n[[k]][[1]][[6]]
    time_list[[k]] <- list()
    for(j in 1:J){
      if(analysis != "QIF") MCLE[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[1]]
      MCLE_mean[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[2]]
      if(analysis=="QIF") {
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
              ((j-1+(k-1)*J)*(p_1)+1):((j-1+(k-1)*J)*(p_1)+p_1)] <- MCLE_psilist_n[[k]][[j]][[3]]
      } else {
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- MCLE_psilist_n[[k]][[j]][[3]]
        psi_g[min(which(subject_indicator==k)):max(which(subject_indicator==k)),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
          MCLE_psilist_n[[k]][[j]][[4]] 
      }
      S[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[5]]
      time_list[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[7]]
    }
    time_list[[k]] <- do.call(rbind, time_list[[k]])
  }
  time_list <- do.call(rbind, time_list)
  time_max <- time_list[which.max(time_list[,1]),]
  
  unordered_psi_g <- psi_g
  unordered_S <- S
  
  print("Computing combined estimate.", quote=FALSE)
  if (method=="exact" & analysis != "QIF"){
    psi_g <- psi_g[,c(order(sapply(comb_scheme, function(x) rep(x,p))),(J*K*p+1):(J*K*p_1))]
    S_all <- do.call(rbind, lapply(unordered_S, function(S_k) do.call(rbind, S_k)))[c(c(matrix(1:(J*K*p_1), nrow=p_1)[1:p,]),
                                                                                      c(matrix(1:(J*K*p_1), nrow=p_1)[(p+1):p_1,])),]
    S_all_reordered <- S_all[c(order(sapply(comb_scheme, function(x) rep(x,p))),(J*K*p+1):(J*K*p_1)),]
    
    S_matrix <- as.matrix(Matrix::bdiag(Matrix::bdiag(
      lapply(sort(unique(comb_scheme)), function(c) {
        S_all_reordered[c(sapply(comb_scheme, function(x) rep(x,p))[order(sapply(comb_scheme, function(x) rep(x,p)))]==c,
                          rep(FALSE, J*K*d)),1:p, drop=FALSE]
      })),
      Matrix::bdiag(
        lapply(1:(J*K), function(c){
          S_all_reordered[(J*K*p+(c-1)*d+1):(J*K*p+c*d),(p+1):p_1, drop=FALSE]
        })
      )))
    
    S_MCLE_mean_list <- mapply(function(S_k,MCLE_k) mapply(function(x,y) x %*% y, S_k, MCLE_k, SIMPLIFY=FALSE), S, MCLE_mean, SIMPLIFY = FALSE)
    S_MCLE_mean_all <- do.call(rbind, lapply(S_MCLE_mean_list, function(S_k) do.call(rbind, S_k)))[c(c(matrix(1:(J*K*p_1), nrow=p_1)[1:p,]),
                                                                                                     c(matrix(1:(J*K*p_1), nrow=p_1)[(p+1):p_1,]))]
    S_MCLE_mean_matrix <- S_MCLE_mean_all[c(order(sapply(comb_scheme, function(x) rep(x,p))),(J*K*p+1):(J*K*p_1))]
    
    V_psi <- t(psi_g)%*%psi_g/N
    output$W <- solve(V_psi)
    
    scale <- t(S_matrix) %*% output$W %*% S_matrix
    main <- t(S_matrix) %*% output$W %*% matrix(S_MCLE_mean_matrix, ncol=1)
    
    output$coefficients <- as.vector(solve(scale)[1:(p*length(unique(comb_scheme))),]%*%main)
    output$vcov <- solve(scale)[1:(p*length(unique(comb_scheme))),1:(p*length(unique(comb_scheme)))]*N
    output$full.coefficients <- as.vector(solve(scale)%*%main)
  } 
  
  if (method=="exact" & analysis == "QIF"){
    psi_g <- unordered_psi_g[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
    S_all_reordered <- t(do.call(rbind, lapply(unordered_S, function(S_k) do.call(rbind, S_k))))[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
    
    rank <- qr(psi_g)$rank
    if(dim(psi_g)[2] != rank){
      psi_g_pca <- prcomp(psi_g, center = FALSE,scale. = FALSE)
      old_psi_g <- psi_g
      old_J <- J
      old_K <- K
      old_S <- S
      ta <- cumsum(table(comb_scheme))
      psi_g <- psi_g_pca$x[,1:rank]
      S_matrix <- (as.matrix(Matrix::bdiag(lapply(sort(unique(comb_scheme)), function(c) {
        S_all_reordered[,sapply(comb_scheme, function(x) rep(x,p_1))[order(sapply(comb_scheme, function(x) rep(x,p_1)))]==c]
      }))) %*% psi_g_pca$rotation )[,1:rank]
      S_MCLE_mean_list <- mapply(function(S_k,MCLE_k) mapply(function(x,y) x %*% y, S_k, MCLE_k, SIMPLIFY=FALSE), S, MCLE_mean, SIMPLIFY = FALSE)
      S_MCLE_mean_matrix <- (t(do.call(rbind, lapply(S_MCLE_mean_list, function(S_k) 
        do.call(rbind, S_k))))[,order(sapply(comb_scheme, function(x) rep(x,p_1)))] %*% 
          psi_g_pca$rotation)[,1:rank]
      
    } else {
      S_matrix <- as.matrix(Matrix::bdiag(lapply(sort(unique(comb_scheme)), function(c) {
        S_all_reordered[,sapply(comb_scheme, function(x) rep(x,p_1))[order(sapply(comb_scheme, function(x) rep(x,p_1)))]==c]
      })))
      S_MCLE_mean_list <- mapply(function(S_k,MCLE_k) mapply(function(x,y) x %*% y, S_k, MCLE_k, SIMPLIFY=FALSE), S, MCLE_mean, SIMPLIFY = FALSE)
      S_MCLE_mean_matrix <- t(do.call(rbind, lapply(S_MCLE_mean_list, function(S_k) 
        do.call(rbind, S_k))))[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
    }
    
    V_psi <- t(psi_g)%*%psi_g/N
    output$W <- solve(V_psi)
    
    scale <- S_matrix %*% output$W %*% t(S_matrix)
    main <- S_matrix %*% output$W %*% matrix(S_MCLE_mean_matrix, ncol=1)
    output$coefficients <- as.vector(solve(scale)%*%main)
    output$vcov <- solve(scale)*N
  } 
  
  if (method=="iterative" | method=="ridge"){
    if (method=="iterative"){
      div <- choose(N,2)
      init_pars <- c(colMeans(do.call(rbind, lapply(MCLE_mean, function(x) do.call(rbind, lapply(x, function(y) y[1:p]))))), 
                     unlist(lapply(MCLE_mean, function(x) lapply(x, function(y) y[(p+1):(p+d)]))))
      optimization <- optim(par=init_pars, 
                            combined.estimation.mean, gr=estimation.deriv.mean, 
                            d=d, W=output$W, response=response, 
                            covariates=covariates, J=J, K=K, N=N, div=div, func=func4, funcderiv=func5, 
                            method="L-BFGS-B", lower=lower_lim, upper=upper_lim, control=list(maxit=500))
      init_pars <- optimization$par
      output$coefficients <- optim(par=init_pars, 
                                   combined.estimation.mean, gr=estimation.deriv.mean, 
                                   d=d, W=output$W, response=response, 
                                   covariates=covariates, J=J, K=K, N=N, div=div, func=func4, funcderiv=func5, 
                                   method="L-BFGS-B", lower=lower_lim, upper=upper_lim)$par
    } 
    
    sensitivity_psi <- matrix(0, J*K*d+p, J*K*(p+d))
    phi_gh <- matrix(0, N, J*K*(p+d))
    for(k in 1:K) {
      for (j in 1:J) {
        block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                          nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
        ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
        block_x <- as.matrix(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,
                                         -which(colnames(ids_block_x)=="id")])
        m <- dim(block_y)[1]
        n <- dim(block_y)[2]
        if(d==2){
          phi_gh_combined <- func4(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], output$coefficients[p+(j-1+(k-1)*J)*d+2], 
                                   block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                 ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[1:p])))
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                 (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[c(p+1, p+2)])))
          
          little_s <- t(func5(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], 
                              output$coefficients[p+(j-1+(k-1)*J)*d+2], block_y, block_x, id, m, n)*n/N)
        }
        if(d==1){
          phi_gh_combined <- func4(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], 
                                   block_y=block_y, block_x=block_x, id=id, m=m, n=n) 
          
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                 ((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[1:p])))
          phi_gh[min(which(subject_indicator==k)):max(which(subject_indicator==k)),
                 (J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- 
            t(sapply(phi_gh_combined, function(x) as.matrix(x[p+1])))
          
          little_s <- t(func5(output$coefficients[1:p], output$coefficients[p+(j-1+(k-1)*J)*d+1], 
                              block_y, block_x, id, m, n)*n/N)
        }
        sensitivity_psi[1:p,((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- little_s[1:p,1:p]
        sensitivity_psi[1:p,(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- little_s[1:p,(p+1):(p+d)]
        sensitivity_psi[((j-1+(k-1)*J)*d+p+1):((j-1+(k-1)*J)*d+p+d),((j-1+(k-1)*J)*(p)+1):((j-1+(k-1)*J)*(p)+p)] <- little_s[(p+1):(p+d),1:p]
        sensitivity_psi[((j-1+(k-1)*J)*d+p+1):((j-1+(k-1)*J)*d+p+d),(J*K*p+(j-1+(k-1)*J)*d+1):(J*K*p+(j-1+(k-1)*J)*d+d)] <- little_s[(p+1):(p+d),(p+1):(p+d)]
      }
    }
    variability_psi <- matrix(0, J*K*(p+d), J*K*(p+d))
    for(i in 1:dim(phi_gh)[1]){
      variability_psi <- variability_psi + phi_gh[i,]%o%phi_gh[i,]
    }
    variability_psi <- variability_psi/N
    
    output$vcov <- solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))%*%
      sensitivity_psi%*%output$W%*%variability_psi%*%output$W%*%t(sensitivity_psi)%*%
      solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))/N
    full_coefficients <- output$coefficients
    output$coefficients <- output$coefficients[1:p]
    output$vcov <- output$vcov[1:p,1:p]
    output$full.coefficients <- full_coefficients
  }
  time2 <- proc.time()-time2
  parallel::stopCluster(sock)
  print("Cluster stopped.", quote=FALSE)
  foreach::registerDoSEQ()
  
  output$family <- family
  output$analysis <- analysis
  output$corstr <- corstr
  output$response_indicator <- response_indicator
  output$subject_indicator <- subject_indicator
  output$J <- J
  output$K <- K
  output$MCLE <- list()
  output$MCLE$beta <- MCLE_mean
  if(analysis != "QIF") output$MCLE$vcov <- MCLE
  output$time <- time1+time2+time_max

  return(output)
}
