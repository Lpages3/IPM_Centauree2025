####### Save selected models for all vital rates ######

#### Initialization
library(spaMM)
library(tidyverse)
library(splines)

setwd("/media/loic/Commun/0Travail/2025 IPM Centauree ISEM/2025 Stage ISEM/Code")

IPM_data <- read.csv("Data/newdata.csv")

IPM_data$Age[IPM_data$Age > 8] <- 8
centauree_data <- IPM_data[!is.na(IPM_data$Size0Mars) & !is.na(IPM_data$Age),]

spaMM.options(separation_max=70)


#### Fecundity
cptldata <- IPM_data[IPM_data$Flowering!=0,]

Cptlglm1 <- fitme(log(Capitule) ~ 1 + Size0Mars + (Age|year), 
                  data=cptldata)

#Modèle avec intercept=0
Cptlglm2 <- fitme(log(Capitule) ~ 0 + Size0Mars + (0+Size0Mars|Pop), 
                  data=cptldata)

#### Survival Probability
survdata <- centauree_data[centauree_data$Flowering!=1,]
survdata1 <- survdata[survdata$Age==1,]
survdata2 <- survdata[survdata$Age>1,]

Survglm11 <- fitme(Survie ~ 1 + bs(Size0Mars,degree=2,df=4) + 
                     (Size0Mars|year),
                   family=binomial,
                   data=survdata1,
                   method="PQL/L")

Survglm12 <- fitme(Survie ~ 1 + bs(Size0Mars,df=3,degree=2) + bs(Age,degree=3,knots=6.5) + 
                     (1|Pop) + (Age|year),
                   family=binomial,
                   data=survdata2,
                   method="PQL/L")


#### Growth
growthdata <- centauree_data[!is.na(centauree_data$Size1Mars), ]
growthdata <- growthdata[growthdata$Size1Mars != 0, ]

Growthglm1 <- fitme(Size1Mars ~ 1 +
                      poly(Size0Mars,3) + bs(Age,degree=2,knots=6.5) +
                      (Size0Mars+Age|year) + (1|Pop),
                    resid.model = ~ log(Size0Mars)+log(Age),
                    data=growthdata)

Growthglm1gamma <- fitme(Size1Mars ~ 1 +
                            poly(Size0Mars,3) + poly(Age,2) +
                            (1|year) + (1|Pop),
                          resid.model = ~ log(Size0Mars)+log(Age),
                          data=growthdata, family=Gamma(log))

Growthglm2 <- fitme(I(Size1Mars**0.4343354) ~ 1 +
                      poly(I(Size0Mars**0.4343354),3) + bs(Age,degree=2,knots=6.5) +
                      (Age|year) + (1|Pop),
                    resid.model = ~ log(Size0Mars)+log(Age)+
                      (1|year) + (1|Pop),
                    data=growthdata)

#### Flowering Probability
reduitdata <- centauree_data[-sample(nrow(centauree_data[centauree_data$Age==1,]),320),]
Flowglm1 <- fitme(Flowering ~  1 + poly(Size0Mars,3) + poly(Age,2) + (Age|Pop),
                  family=binomial,
                  data=centauree_data, method="PQL/L")

Flowglm2 <- fitme(Flowering ~  1 + poly(Size0Mars,3) + poly(Age,2) + (Age|Pop),
                  family=binomial,
                  data=reduitdata, method="PQL/L")


#### Seedling size distribution
plantule_data <- centauree_data[centauree_data$Age==1,]

Pltglm1 <- fitme(Size0Mars ~ 1 + (1|year) + (1|Pop) + (1|Pop:year), 
                 data=plantule_data,
                 family = Gamma(log))


### Establishment rate

Cptlglm1 <- fitme(log(Capitule) ~ 1 + Size0Mars+ (Age|year),
                  data=cptldata)
# NbrCptl = exp(2.31070+0.06846*Size0Mars)
cptl_data_predi <- cptldata %>%
  mutate(Capitule = ifelse(is.na(Capitule), exp(2.31070+0.06846*Size0Mars), Capitule))

#Jeu de données 2 : tous les quadrats, y compris ceux où il n'y a pas de plantules
cptl_data_predi2 <- IPM_data %>%
  mutate(Capitule = ifelse(is.na(Capitule), exp(2.31070+0.06846*Size0Mars), Capitule))

plt <- IPM_data %>% 
  filter(Age==1) %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(Seedlings = sum(Age))

cptl <- cptl_data_predi %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(Capitule = sum(Capitule))

cptl2 <- cptl_data_predi2 %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(Capitule = sum(Capitule))

#Calcul direct de l'establishment rate, puis fit sur ces données
Estb <- inner_join(plt,cptl, by=join_by(Quadrat,year,Pop))

Estb <- Estb %>% mutate(EstbRate=rep(NA)) %>%
  arrange(Quadrat)

for (i in 2:length(Estb$Quadrat)){
  if (Estb$Quadrat[i]!=Estb$Quadrat[i-1]){next}
  if (Estb$year[i]!=Estb$year[i-1]+1){next}
  Estb$EstbRate[i] <- Estb$Seedlings[i]/Estb$Capitule[i-1]
}
Estb <- Estb %>%
  group_by(Pop,year) %>%
  mutate(EstbPop = mean(EstbRate,na.rm=TRUE))

Estbglm1 <- fitme(EstbRate ~ 1 , data=Estb)

#Fit du nombre de plantules à t+1
cptl$year <- cptl$year+1
cptl2$year <- cptl2$year+1
Estb <- left_join(plt, cptl, by = join_by(Quadrat, year, Pop))
Estb2 <- left_join(plt, cptl2, by = join_by(Quadrat, year, Pop))


Estbglm2 <- fitme(Seedlings ~ 1 + offset(log(Capitule)) + (1|Pop:year),
                  data = Estb,
                  family = Poisson(log))

# Estbglm3 <- fitmv(submodels=list(list(as.integer(Capitule) ~ 1+(0+mv(1,2)|Pop:year),family=poisson()),
#                                  list(Seedlings ~ 1+(0+mv(1,2)|Pop:year),family=poisson())), 
#                   data=Estb)
# 
# 
# Clapefn <- function(v, data, 
#                     formula, # the formula will contain offset(logfixef)
#                     rand.family, 
#                     return.fit=FALSE) {
#   fixef <- v[1] + v[2] * data$Capitule
#   data$logfixef <- log(fixef)
#   fit <- fitme(formula=formula, data=data,
#                rand.family=rand.family, family=poisson(log) )
#   if (return.fit) {fit} else {-logLik(fit)}
# }
# lower <- c(1e-6,1e-6)
# upper <- c(200,200)
# optr <- (spaMM::.safe_opt(
#   init=c(1,0.1),objfn = Clapefn,
#   data=Estb2, lower = lower, upper=upper,
#   verbose = FALSE,LowUp = list(lower=lower, upper=upper),
#   formula=Seedlings ~ offset(logfixef) + 0 + (1|Pop:year), # Note the +0+ ...
#   rand.family=Gamma(log),
#   return.fit=FALSE
# ))
# Seedglmm <- Clapefn(optr$solution, data=Estb2, 
#                      formula=Seedlings ~ offset(logfixef) + 0 + (1|Pop:year),
#                      rand.family=Gamma(log), return.fit=TRUE) 

#### Save all models
save(Survglm11,
     Survglm12,
     Cptlglm1,
     Growthglm1,
     Growthglm1gamma,
     Growthglm2,
     Flowglm1,
     Flowglm2,
     Pltglm1,
     Estbglm1,
     Estbglm2,
     Seedglmm,
     file="IPM/ModelsAIC.RData")

obs_beta=as.numeric(Flowglm1$fixef[1])
se_obs_beta=as.numeric(sqrt(diag(vcov(Flowglm1)))[1])

save(obs_beta, file = "IPM/obs_beta.RData")
save(se_obs_beta, file = "IPM/se_obs_beta.RData")
