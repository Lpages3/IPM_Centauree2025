####### Creation of a fake data set with responses predicted by the models ######

#### Initialization
library(spaMM)
library(tidyverse)
library(splines)

setwd("/media/loic/Commun/0Travail/Stage 2025 ISEM/Models/IPM")

centauree_data <- read.csv("donnesIPM_short.csv")
centauree_data_complet <- read.csv("donnesIPM.csv")

#Supprimer plantes dont l'age est inconnu
centauree_data <- centauree_data[!is.na(centauree_data$age0), ]
centauree_data$age1 <- ifelse(centauree_data$Stage1=="V",centauree_data$age0+1,NA)

#Forcer l'age maximal à 8
centauree_data$age0[centauree_data$age0 > 8] <- 8

spaMM.options(separation_max=70)

#### AIC
# Fecondity
cptldata <- centauree_data_complet[centauree_data_complet$Flowering0!=0,]

Cptlglm1 <- fitme(Cptl0 ~ 1 + Size0Mars, 
                  data=cptldata)

# Survival Probability
survdata <- centauree_data[centauree_data$Flowering0!=1,]
survdata1 <- survdata[survdata$age0==1,]
survdata2 <- survdata[survdata$age0>1,]

Survglm11 <- fitme(SurvieMars ~ 1 + poly(Size0Mars,3) + 
                     (Size0Mars|year) + (1|Pop),
                   family=binomial,
                   data=survdata1,
                   method="PQL/L")

Survglm12 <- fitme(SurvieMars ~ 1 + bs(Size0Mars,df=4,degree=2) +(bs(age0,degree=3,knots=6.5)) + 
                     (age0|year) + (1|Pop),
                   family=binomial,
                   data=survdata2,
                   method="PQL/L")

# Survglm1 <- fitme(SurvieMars ~ 1+ bs(Size0Mars,df=5,degree=3) + bs(age0,degree=3,knots = 6.5)+ (Size0Mars|year) + (age0|Pop) ,
#                   family=binomial, 
#                   data=survdata,
#                   method="PQL/L")

# Growth
growthdata <- centauree_data[!is.na(centauree_data$Size1Mars), ]
growthdata <- growthdata[growthdata$Size1Mars != 0, ]

Growthglm1 <- fitme(log(Size1Mars) ~ 1 + poly(log(Size0Mars),4) + poly(age0,3) + (log(Size0Mars)|year) + (log(Size0Mars)|Pop),
                    data=growthdata)


# Flowering Probability
Flowglm1 <- fitme(Flowering0 ~  1 + poly(Size0Mars,3) + poly(age0,2) + (age0|Pop),
                  family=binomial,
                  data=centauree_data, method="PQL/L")

# Save all models
save(Survglm11,Survglm12,Cptlglm1,Growthglm1,Flowglm1, file="ModelsAIC")

obs_beta=as.numeric(Flowglm1$fixef[1])
se_obs_beta=as.numeric(sqrt(diag(vcov(Flowglm1)))[1])

save(obs_beta, file = "obs_betaAIC")
save(se_obs_beta, file = "se_obs_betaAIC")


#### BIC
# Fecondity
cptldata <- centauree_data_complet[centauree_data_complet$Flowering0!=0,]

Cptlglm1 <- fitme(Cptl0 ~ 1 + Size0Mars, 
                  data=cptldata)

# Survival Probability
survdata <- centauree_data[centauree_data$Flowering0!=1,]
survdata1 <- survdata[survdata$age0==1,]
survdata2 <- survdata[survdata$age0>1,]

Survglm11 <- fitme(SurvieMars ~ 1 + poly(Size0Mars,3) + 
                     (Size0Mars|year) + (1|Pop),
                   family=binomial,
                   data=survdata1,
                   method="PQL/L")

Survglm12 <- fitme(SurvieMars ~ 1 + bs(Size0Mars,df=4,degree=2) +(bs(age0,degree=3,knots=6.5)) + 
                     (age0|year) + (1|Pop),
                   family=binomial,
                   data=survdata2,
                   method="PQL/L")

# Survglm1 <- fitme(SurvieMars ~ 1+ bs(Size0Mars,df=5,degree=3) + bs(age0,degree=3,knots = 6.5)+ (Size0Mars|year) + (age0|Pop) ,
#                   family=binomial, 
#                   data=survdata,
#                   method="PQL/L")

# Growth
growthdata <- centauree_data[!is.na(centauree_data$Size1Mars), ]
growthdata <- growthdata[growthdata$Size1Mars != 0, ]

Growthglm1 <- fitme(log(Size1Mars) ~ 1 + poly(log(Size0Mars),4) + poly(age0,3) + (log(Size0Mars)|year) + (log(Size0Mars)|Pop),
                    data=growthdata)


# Flowering Probability
Flowglm1 <- fitme(Flowering0 ~  1 + poly(Size0Mars,3) + poly(age0,2) + (age0|Pop),
                  family=binomial,
                  data=centauree_data, method="PQL/L")

# Save all models
save(Survglm11,Survglm12,Cptlglm1,Growthglm1,Flowglm1, file="ModelsBIC")

obs_beta=as.numeric(Flowglm1$fixef[1])
se_obs_beta=as.numeric(sqrt(diag(vcov(Flowglm1)))[1])

save(obs_beta, file = "obs_betaBIC")
save(se_obs_beta, file = "se_obs_betaBIC")