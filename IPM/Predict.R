####### Save selected models for all vital rates ######

#### Initialization
library(spaMM)
library(tidyverse)
library(splines)

setwd("/media/loic/Commun/0Travail/Stage 2025 ISEM/Code/IPM")

IPM_data <- read.csv("newdata.csv")

IPM_data$Age[IPM_data$Age > 8] <- 8
centauree_data <- IPM_data[!is.na(IPM_data$Size0Mars) & !is.na(IPM_data$Age),]

spaMM.options(separation_max=70)


#### Fecundity
cptl_data <- IPM_data[IPM_data$Flowering!=0,]

Cptlglm1 <- fitme(log(Capitule) ~ 1 + Size0Mars + (Age|year), 
                  data=cptl_data)

#Modèle avec intercept=0
Cptlglm2 <- fitme(log(Capitule) ~ 0 + Size0Mars + (0+Size0Mars|Pop), 
                  data=cptl_data)

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

plt <- IPM_data %>% 
  filter(Age==1) %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(NombrePlantules = sum(Age))

cptl <- cptl_data_predi %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(Capitule = sum(Capitule))

#Calcul avec la formule d'origine utilisant les données brutes
Estb <- inner_join(plt,cptl, by=join_by(Quadrat,year,Pop))

Estb <- Estb %>% mutate(EstbRate=rep(NA)) %>%
  arrange(Quadrat)

for (i in 2:length(Estb$Quadrat)){
  if (Estb$Quadrat[i]!=Estb$Quadrat[i-1]){next}
  if (Estb$year[i]!=Estb$year[i-1]+1){next}
  Estb$EstbRate[i] <- Estb$NombrePlantules[i]/Estb$Capitule[i-1]
}
Estb <- Estb %>%
  group_by(Pop,year) %>%
  mutate(EstbPop = mean(EstbRate,na.rm=TRUE))

Estbglm1 <- fitme(EstbRate ~ 1 , data=Estb)

Estb <- left_join(plt, cptl, by = join_by(Quadrat, year, Pop))

Estbglm2 <- fitme(NombrePlantules ~ 1 + offset(log(Capitule)) + (1|Pop:year),
                  data = Estb,
                  family = Poisson(log))

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
     file="ModelsAIC")

obs_beta=as.numeric(Flowglm1$fixef[1])
se_obs_beta=as.numeric(sqrt(diag(vcov(Flowglm1)))[1])

save(obs_beta, file = "obs_beta")
save(se_obs_beta, file = "se_obs_beta")
