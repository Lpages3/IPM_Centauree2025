####### Save selected models for all vital rates ######

#### Initialization
library(spaMM)
library(tidyverse)
library(splines)

setwd("/media/loic/Commun/0Travail/Stage 2025 ISEM/Models/IPM")

IPM_data <- read.csv("newdata.csv")

IPM_data$Age[IPM_data$Age > 8] <- 8
centauree_data <- IPM_data[!is.na(IPM_data$Size0Mars) & !is.na(IPM_data$Age),]

spaMM.options(separation_max=70)


#### Fecundity
cptldata <- IPM_data[IPM_data$Flowering!=0,]

Cptlglm1 <- fitme(log(Capitule) ~ 1 + Size0Mars + (Age|year), 
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
                     (1|Pop),
                   family=binomial,
                   data=survdata2,
                   method="PQL/L")


#### Growth
growthdata <- centauree_data[!is.na(centauree_data$Size1Mars), ]
growthdata <- growthdata[growthdata$Size1Mars != 0, ]

Growthglm2 <- fitme(Size1Mars ~ 1 +
                      poly(Size0Mars,3) + bs(Age,degree=2,knots=6.5) +
                      (Size0Mars+Age|year) + (1|Pop),
                    resid.model = ~ log(Size0Mars)+log(Age),
                    data=growthdata)

Growthglm1 <- fitme(log(Size1Mars) ~ 1 + poly(log(Size0Mars),4) + poly(Age,3) + (log(Size0Mars)|year) + (log(Size0Mars)|Pop),
                    data=growthdata)

Growthglm12 <- fitme(log(Size1Mars) ~ 1 +
                       poly(Size0Mars,3) + poly(Age,3) +
                       (Size0Mars+Age|year) + (Size0Mars|Pop),
                     data=growthdata, resid.model = ~log(Size0Mars))


#### Flowering Probability

Flowglm1 <- fitme(Flowering ~  1 + poly(Size0Mars,3) + poly(Age,2) + (Age|Pop),
                  family=binomial,
                  data=centauree_data, method="PQL/L")


#### Seedling size distribution
plantule_data <- centauree_data[centauree_data$Age==1,]

Pltglm1 <- fitme(Size0Mars ~ 1 + (1|year) + (1|Pop) + (1|Pop:year), 
                 data=plantule_data,
                 family = Gamma(log))


#### Establishment rate

# Cptlglm1 <- fitme(log(Capitule) ~ 1 + Size0Mars+ (Age|year), 
#                   data=cptldata)
# NbrCptl = exp(2.31070+0.06846*Size0Mars)
# cptl_data_predi <- cptldata %>% 
#   mutate(Capitule = ifelse(is.na(Capitule), exp(2.31070+0.06846*Size0Mars), Capitule))
# 
# plt <- centauree_data_complet %>% 
#   filter(Age==1) %>% 
#   group_by(Quadrat,year,Pop) %>% 
#   summarize(NombrePlantules = sum(Age))
# 
# cptl <- cptl_data_predi %>% 
#   group_by(Quadrat,year,Pop) %>% 
#   summarize(NombreCapitules = sum(Capitule))
# 
# #Calcul avec le fs0 du modèle matriciel
# Estb <- inner_join(plt,cptl, by=join_by(Quadrat,year,Pop)) %>% 
#   group_by(year,Pop) %>% 
#   mutate(Cptl = mean(NombreCapitules)) %>% 
#   arrange(year)
# EstYear <- NULL
# for (i in 1:27){
#   EstYear[i] <- fs0M[i]/Estb$Cptl[1+6*(i-1)]
# }
# EstYear
# mean(EstYear)
# 
# #Calcul avec la formule d'origine utilisant les données brutes
# Estb <- inner_join(plt,cptl, by=join_by(Quadrat,year,Pop))
# 
# Estb <- Estb %>% mutate(EstbRate=rep(NA)) %>% 
#   arrange(Quadrat)
# 
# for (i in 2:length(Estb$Quadrat)){
#   if (Estb$Quadrat[i]!=Estb$Quadrat[i-1]){next}
#   if (Estb$year[i]!=Estb$year[i-1]+1){next}
#   Estb$EstbRate[i] <- Estb$NombrePlantules[i]/Estb$NombresCapitules[i-1]
# }
# Estb <- Estb %>%
#   group_by(Pop,year) %>% 
#   mutate(EstbPop = mean(EstbRate,na.rm=TRUE))
# 
# Estbglm1 <- fitme(EstbRate ~ 1 + (1|Pop:year), data=Estb)

#### Save all models
save(Survglm11,
     Survglm12,
     Cptlglm1,
     Growthglm1,
     Growthglm12,
     Growthglm2,
     Flowglm1,
     Pltglm1, file="ModelsAIC")

obs_beta=as.numeric(Flowglm1$fixef[1])
se_obs_beta=as.numeric(sqrt(diag(vcov(Flowglm1)))[1])

save(obs_beta, file = "obs_beta")
save(se_obs_beta, file = "se_obs_beta")