IPM_data <- read.csv("newdata.csv")

centauree_data <- IPM_data[!is.na(IPM_data$Size0Mars) & !is.na(IPM_data$Age),]
centauree_data$Age[centauree_data$Age > 8] <- 8
growthdata <- centauree_data[!is.na(centauree_data$Size1Mars), ]
growthdata <- growthdata[growthdata$Size1Mars != 0, ]



bxcx <- function(object, interval=c(-2,2), verbose=FALSE) {
  y <- object$y
  objfn <- function(lambda, return_fit=FALSE) {
    if (lambda==0) {
      newy <- log(y)
    } else newy <- (y^lambda-1)/lambda
    refit <- update_resp(object, newresp=newy)
    if (return(fit)) {
      return(refit)
    } else {
      if (verbose) print(logLik(refit))
      return( - logLik(refit))
    }
  }
  
  optr <- optim(0,fn=objfn,interval=interval)
  objfn(lambda = optr$par, return_fit = TRUE)
}

boxcx <- function(object, interval=c(-2,2), verbose=FALSE) {
  data <- object$data
  Size0Mars <- data$Size0Mars
  Size1Mars <- data$Size1Mars
  objfn <- function(lambda, return_fit=FALSE) {
    if (lambda==0) {
      data$Size0Mars <- log(Size0Mars)
      data$Size1Mars <- log(Size1Mars)
    } else {
      data$Size0Mars <- (Size0Mars^lambda-1)/lambda
      data$Size1Mars <- (Size1Mars^lambda-1)/lambda
    }
    refit <- update(object, data=data)
    if (return_fit) {
      return(refit)
    } else {
      logjac <- (lambda-1)*sum(log(Size1Mars))
      resu <- logLik(refit)+logjac
      if (verbose) print(resu)
      return( - resu)
    }
  }
  
  optr <- optimize(f=objfn,interval=interval)
  return(list(optr=optr,
              fit=objfn(lambda = optr$minimum, return_fit = TRUE),
              objfn=objfn))
}

powtr <- function(object, interval=c(-0.1,0.9), verbose=FALSE) {
  data <- object$data
  Size0Mars <- data$Size0Mars
  Size1Mars <- data$Size1Mars
  objfn <- function(lambda, return_fit=FALSE) {
    data$Size0Mars <- Size0Mars^lambda
    data$Size1Mars <- Size1Mars^lambda
    refit <- update(object, data=data)
    if (return_fit) {
      return(refit)
    } else {
      logjac <- sum((lambda-1)*log(Size1Mars)+log(lambda))
      resu <- logLik(refit)+logjac
      if (verbose) print(resu)
      return( - resu)
    }
  }
  
  optr <- optimize(f=objfn,interval=interval)
  return(list(optr=optr,
              fit=objfn(lambda = optr$minimum, return_fit = TRUE),
              objfn=objfn))
}



growthdata$size0resid <- growthdata$Size0Mars

(Growthglm2b <- fitme(Size1Mars**0.4242424 ~ 1 +
                        poly(Size0Mars**0.4242424,3) + bs(Age,degree=2,knots=6.5) +
                        (Age|year) + (1|Pop),
                      resid.model = ~ log(Size0Mars)+log(Age),
                      data=growthdata))

(essaiLMMfull <- fitme(Size1Mars ~ 1 +
                     poly(Size0Mars,3) + bs(Age,degree=2,knots=6.5) +
                     (Age|year) + (1|Pop),
                   resid.model = ~ log(size0resid)+log(Age)+
                     (1|year) + (1|Pop),
                   data=growthdata))
(essaiLMM <- fitme(Size1Mars ~ 1 +
                     poly(Size0Mars,3) + bs(Age,degree=2,knots=6.5) +
                     (1|year) + (1|Pop),
                   resid.model = ~ log(size0resid)+log(Age),
                   data=growthdata))
gof(essaiLMM)
ESSAI <- ESSAIxx <- boxcx(essaiLMMfull,verbose=TRUE)
ESSAI <- ESSAIpo <- powtr(essaiLMMfull,verbose=TRUE)
ESSAI2 <- powtr(essaiLMM)
gof(ESSAI$fit)
gof(Growthglm2b)
ESSAI$objfn(0.4242424, return_fit = TRUE)
ESSAI$objfn(ESSAI$optr$minimum, return_fit = TRUE)
