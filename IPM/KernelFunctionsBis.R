#Fecundity
numbofcap <- function(x,a) {
  fake_data$Size0Mars <- unique(x)
  sortie <- fake_data %>% mutate(CapituleNbr = predict(Cptlglm1, newdata = fake_data)[,1])
  sortie2 <- mean(sortie$CapituleNbr[fake_data$age0==a])
  return(matrix(sortie2,MatrixDim,MatrixDim,byrow = T ))
}


numbofcap=function(x,a) {
  fake_data$age0 <- a  
  fake_data$Size0Mars <- unique(x)
  sortie=predict(Cptlglm1, fake_data)
  sortie2=aggregate(sortie,list(fake_data$Size0Mars),mean)
  return(matrix(sortie2$V1,MatrixDim,MatrixDim,byrow = T ))
}

#Survival Probability
sx0 <- function(x,a) {
  fake_data$Size0Mars <- unique(x)
  sortie <- fake_data %>% mutate(Survival = predict(Survglm1, newdata = fake_data)[,1])
  sortie2 <- mean(sortie$Survival[fake_data$age0==a])
  return(sortie2)}


#Flowering Probability 
#Beta = will be used to modify the intercept of the flowering function
#obs_beta = observed value
#extract Beta0
load("obs_beta")
load("se_obs_beta")


#Flowering function for The survival-growth kernel
flr0 <- function(x,a,Beta=0) {
  fake_data$Size0Mars <- unique(x)
  sortie <- fake_data %>% mutate(Flowering = predict(Flowglm1, newdata = fake_data)[,1] + Beta)
  sortie2 <- mean(sortie$Flowering[fake_data$age0==a])
  return(sortie2)}


#Flowering function for fecundity fyx0 - same function the difference is in the format of the output
flr1 <- function(x,a,Beta=0) {
  fake_data$Size0Mars <- unique(x)
  sortie <- fake_data %>% mutate(Flowering = predict(Flowglm1, newdata = fake_data)[,1] + Beta)
  sortie2 <- mean(sortie$Flowering[fake_data$age0==a])
  return(matrix(sortie2,MatrixDim,MatrixDim,byrow = T ))
}


#Growth function
Gyx0 <- function (y,x,a) {
  fake_data$Size0Mars <- unique(x)
  sortie <- fake_data %>% mutate(Size1 = exp(predict(Growthglm1, newdata = fake_data)[,1]))
  sortie2 <- mean(sortie$Size1[fake_data$age0==a])  
  M <- matrix(sortie2,MatrixDim,MatrixDim, byrow=T)
  SD <- sd(sortie$Size1[fake_data$age0==a])
  return(dnorm(y,mean=M,sd=SD))
}


#Combine fecundity and flowering probability 
fyx0 <- function(y, x ,a, Estbl, Beta=0, intervalle) {
  #flowering probability:
  p.flow <-flr1(x,a, Beta)
  #number of capitula per flowering plant
  n.captl <-numbofcap(x,a)
  #Seedlings Size distribution
  ProporSeedlSize=densSeedl(y, intervalle)
  
  Fyx0=ifelse((p.flow*n.captl*Estbl*ProporSeedlSize) >=0,p.flow*n.captl*Estbl*ProporSeedlSize,0)
  return(Fyx0)
}
