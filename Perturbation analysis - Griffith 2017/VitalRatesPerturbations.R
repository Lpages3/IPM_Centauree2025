setwd("/media/loic/Commun/0Travail/Stage 2025 ISEM/Code/IPM")

#Load packages
library(spaMM)
library(popbio)
library(splines)
library(tidyverse)
#Load functions

source("/media/loic/Commun/0Travail/Stage 2025 ISEM/Code/Perturbation analysis - Griffith 2017/Perturbation Functions (Copy).R")

#Load Data
IPM_data <- read.csv("newdata.csv")
AgeMax <- 8
IPM_data$Age[IPM_data$Age > AgeMax] <- AgeMax
centauree_data <- IPM_data[!is.na(IPM_data$Size0Mars) & !is.na(IPM_data$Age),]

#Kernel options

#Kernel options
MatrixDim <- 50
RecrRate <- 0.338

n <- MatrixDim #Kernel dimensions
Iterations <- 1000
p <- 1e-5 #numerical perturbation amount
minsize <-  0.5
maxsize <-  25
maxsize1 <- 15
n.size <- n
h <-  (maxsize - minsize) / n.size
h1 <-  (maxsize1 - minsize) / n.size
b <- minsize + c(0:n.size) * h

#sizes for an individual of age a
#midpoint
ymid.a = 0.5 * (b[1:n.size] + b[2:(n.size + 1)])
#intervals
yI.a <- matrix(nrow = n.size, ncol = 2)
for (i in 1:n.size){
  yI.a[i,1] <- b[i]
  yI.a[i,2] <- b[i+1]
}

b1 = minsize + c(0:n.size) * h1
#sizes for an individual of age 1
#midpoint
ymid.1 = 0.5 * (b1[1:n.size] + b1[2:(n.size + 1)])
#intervals
yI.1 <- matrix(nrow = n.size, ncol = 2)
for (i in 1:n.size){
  yI.1[i,1] <- b1[i]
  yI.1[i,2] <- b1[i+1]
}

#Load fitted model objects
load("ModelsAIC")

## Creation of fake data
AgeMax <- 8
year <- 1995:2022
Pop <- c("Po","Au","Pe","E1","E2","Cr")
taille_range <- seq(0.5, 25, by = 0.5) 
age_range <- 1:AgeMax
# 
# fake_data <- expand.grid(
#   year = annees,
#   Pop = populations,
#   Size0Mars = NA,
#   Age = NA)


Pop=gl(6,length(year),labels = Pop)
rer=data.frame(Pop,year)
rer$Index=1:(length(year)*6)
fake_data=do.call("rbind", replicate(MatrixDim, rer, simplify = FALSE))
fake_data=fake_data[order(fake_data$Index),]
fake_data$Size0Mars=NA
fake_data$Age=NA

rm(rer, Pop, year)

annees <- 1995:2022
populations <- c("Po","Au","Pe","E1","E2","Cr")

#Fecundity
numbofcap <- function(x, a) {
  fake_data <- expand.grid(
    year = annees,
    Pop = populations,
    Size0Mars = x,
    Age = a)
  mu <- exp(predict(Cptlglm1, fake_data))
  mu2 <- aggregate(mu, list(fake_data$Size0Mars), mean)
  nbcap <- ifelse(mu2$V1>=0,mu2$V1,0)
  return(matrix(nbcap, length(x), length(x), byrow = T))
}

#Survival Probability for age 1
sx01 <- function(x,a) {
  fake_data <- expand.grid(
    year = annees,
    Pop = populations,
    Size0Mars = x,
    Age = a)
  mu <- predict(Survglm11, fake_data, allow.new.levels = T)
  mu2 <- aggregate(mu, list(fake_data$Size0Mars), mean)
  return(mu2$V1)
}

#Survival Probability for age > 1
sx02 <- function(x, a) {
  fake_data <- expand.grid(
    year = annees,
    Pop = populations,
    Size0Mars = x,
    Age = a)
  mu <- predict(Survglm12, fake_data, allow.new.levels = T)
  mu2 <- aggregate(mu, list(fake_data$Size0Mars), mean)
  return(mu2$V1)
}

#Flowering function for The survival-growth kernel
flr0 <- function(x, a) {
  fake_data <- expand.grid(
    year = annees,
    Pop = populations,
    Size0Mars = unique(x),
    Age = a)
  mu <- plogis(predict(Flowglm1, newdata = fake_data, allow.new.levels = T, type = "link") )
  mu2 <- aggregate(mu, list(fake_data$Size0Mars), mean)
  return(mu2$V1)
}

#Flowering function for fecundity fyx0 - same function the difference is in the format of the output
flr1 <- function(x, a) {
  fake_data <- expand.grid(
    year = annees,
    Pop = populations,
    Size0Mars = unique(x),
    Age = a)
  mu <- plogis(predict(Flowglm1, fake_data, allow.new.levels = T, type = "link") )
  mu2 <- aggregate(mu, list(fake_data$Size0Mars), mean)
  return(matrix(mu2$V1, MatrixDim, MatrixDim, byrow = T))
}

#Growth function
Gyx1 <- function(ysup,yinf,x,a){ # Croissance CDF sans log 
  yinf[1] <- -Inf
  ysup[length(ysup)] <- Inf
  
  fake_data <- expand.grid(
    year = annees,
    Pop = populations,
    Size0Mars = unique(x),
    Age = a)
  mu <- predict(Growthglm2, fake_data, allow.new.levels = T) #Prédit les taille à t+1 possibles 
  mu2 <- aggregate(mu, list(fake_data$Size0Mars), mean) #Calcule la moyenne des taillet+1
  usd <- sqrt(exp(1.0462608+0.5367287*log(unique(x))))
  
  cdf.sup <- mapply(pnorm, mean = mu2$V1, sd = usd, MoreArgs = list(q=ysup))
  cdf.inf <- mapply(pnorm, mean = mu2$V1, sd = usd, MoreArgs = list(q=yinf))
  # cdf.sup <- mapply(pnorm, q=ysup, MoreArgs = list(mean = mu2$V1, sd = usd))
  # cdf.inf <- mapply(pnorm, q=yinf, MoreArgs = list(mean = mu2$V1, sd = usd))
  return(cdf.sup - cdf.inf)}

#Seedling size density
densSeedl <- function(ysup, yinf) {
  yinf[1] <- -Inf
  ysup[length(ysup)] <- Inf
  mu <- predict(Pltglm1, newdata = fake_data)[, 1]
  shape <- 1 / Pltglm1$phi
  mu <- mean(mu)
  scale <- mu * Pltglm1$phi
  cdf.sup <- pgamma(ysup, shape = shape, scale = scale)
  cdf.inf <- pgamma(yinf, shape = shape, scale = scale)
  return(cdf.sup - cdf.inf)
}

#Make IPM
load("IPMKernal_AIC")
Kernel.Mid.3 <- IPMKernal_AIC

allEV <- expand.grid(Age=1:8,
                     Growth=NA,
                     Retrogression=NA,
                     Survival=NA,
                     Flowering=NA,
                     Fecundity=NA,
                     SeedlingSize_small=NA,
                     SeedlingSize_large=NA,
                     Estb=NA)
a <- 2
age <- a
for(a in 1:8){
ymid <- ymid.1
if(a>1){ymid <- ymid.a}

#Get mean values for discretized vital rate functions
S <- sx01(ymid.a,a)
Fl <- flr0(ymid.a, a)
G <- Gyx1(ysup=yI.a[,2], yinf=yI.a[,1], x=ymid, a=age)
R <- numbofcap(ymid.a, a)
OG <- densSeedl(ysup=yI.1[,2],yinf=yI.1[,1])
Er <- RecrRate

# # - Perturbation Apporach 3 - # # ----

#Sensitivity and elasticity
Pert3.S <- do.Pert3.S(Kernel.Mid.3)
Pert3.Fl <- do.Pert3.Fl(Kernel.Mid.3)
Pert3.G <- do.Pert3.G(Kernel.Mid.3)
Pert3.R <- do.Pert3.R(Kernel.Mid.3)
Pert3.OG <- do.Pert3.OG(Kernel.Mid.3)
Pert3.Er <- do.Pert3.Er(Kernel.Mid.3)

#Small and large offspring size indexing
mean.OG <- mean(centauree_data$Size0Mars[centauree_data$Age==1])
mean.OG.z <- min(which(ymid.1>mean.OG))
OG.sm <- 1:(mean.OG.z-1)
OG.lg <- mean.OG.z:n

#Splitting/summing growth elasticities
Elas3.G.grow.pos <- Pert3.G$EV * lower.tri(Pert3.G$EV, diag=TRUE) * (Pert3.G$EV>=0) 
Elas3.G.grow.neg <- Pert3.G$EV * lower.tri(Pert3.G$EV, diag=TRUE) * (Pert3.G$EV<0) 
Elas3.G.retro.pos <- Pert3.G$EV * upper.tri(Pert3.G$EV) * (Pert3.G$EV>=0) 
Elas3.G.retro.neg <- Pert3.G$EV * upper.tri(Pert3.G$EV) * (Pert3.G$EV<0)

#Plot

allEV[allEV$Age==a,][-1] <- c(sum(Elas3.G.grow.pos), sum(Elas3.G.retro.neg),
                    sum(Pert3.S$EV), sum(Pert3.Fl$EV), sum(Pert3.R$EV),
                    sum(Pert3.OG$EV[OG.lg,]), sum(Pert3.OG$EV[OG.sm,]), sum(Pert3.Er$EV))
}
allEV_long <- pivot_longer(allEV, cols = -Age, names_to = "Variable", values_to = "Value")

# Création du bar plot
library(RColorBrewer)

ggplot(allEV_long, aes(x = Variable, y = Value, fill = Variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Age, ncol = 4) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Summed elasticity value by age",
       x = "Vital rates",
       y = "Ealsticity",
       fill = "Vital rates") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(allEV_long, aes(x = Age, y = Value, fill = as.factor(Age))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "Summed elasticity value by vital rate",
       x = "Age",
       y = "Elasticity",
       fill = "Age") +  
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  theme_minimal()