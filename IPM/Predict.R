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
length(centauree_data$age0[centauree_data$age0 >= 8])
centauree_data$age0[centauree_data$age0 > 8] <- 8

spaMM.options(separation_max=70)

## Creation of fake data
AgeMax <- 8
annees <- 1995:2022
populations <- c("Po","Au","Pe","E1","E2","Cr")
taille_range <- seq(0.5, 25, by = 0.5) 
age_range <- 1:AgeMax

fake_data <- expand.grid(
  year = annees,
  Pop = populations,
  Size0Mars = taille_range,
  age0 = age_range)

fake_data1 <- fake_data[fake_data$age0==1,]
fake_data2 <- fake_data[fake_data$age0>1,]

## Add predicted responses

# Fecondity
cptldata <- centauree_data_complet[centauree_data_complet$Flowering0!=0,]

Cptlglm1 <- fitme(Cptl0 ~ 1 + Size0Mars, 
                  data=cptldata)
Cptlpredict1 <- predict(Cptlglm1, newdata = fake_data)[,1]

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
# Survpredict1 <- predict(Survglm1, newdata = fake_data)[,1]

# Growth
growthdata <- centauree_data[!is.na(centauree_data$Size1Mars), ]
growthdata <- growthdata[growthdata$Size1Mars != 0, ]
growthdata <- growthdata[!is.na(growthdata$age0),]

Growthglm1 <- fitme(log(Size1Mars) ~ 1 + poly(log(Size0Mars),4) + poly(age0,3) + (log(Size0Mars)|year) + (log(Size0Mars)|Pop),
                    data=growthdata, resid.model = ~log(Size0Mars))
Growthglm2 <- fitme(Size1Mars ~ 1 +
                      bs(Size0Mars,degree=3,df=5) + bs(age0,degree=2,knots=6.5) +
                      (Size0Mars+age0|year) + (Size0Mars|Pop),
                    resid.model = ~ log(Size0Mars),
                    data=growthdata)

Growthpredict1 <- predict(Growthglm1, newdata = fake_data)[,1]

# Flowering Probability
Flowglm1 <- fitme(Flowering0 ~  1 + poly(Size0Mars,3) + poly(age0,2) + (age0|Pop),
                  family=binomial,
                  data=centauree_data, method="PQL/L")
Flowpredict1 <- predict(Flowglm1, newdata = fake_data)[,1]


Pltglm1 <- fitme(Size0Mars ~ 1 + (1|year) + (1|Pop) + (1|Pop:year), 
                 data=plantule_data,
                 family = Gamma(log))
# Save all models
save(Survglm11,Survglm12,Cptlglm1,Growthglm1,Growthglm2,Flowglm1,Pltglm1, file="ModelsAIC")

obs_beta=as.numeric(Flowglm1$fixef[1])
se_obs_beta=as.numeric(sqrt(diag(vcov(Flowglm1)))[1])

save(obs_beta, file = "obs_beta")
save(se_obs_beta, file = "se_obs_beta")

# 
# # Add to fake_data
# predi_data <- fake_data %>% mutate(CapituleNbr = Cptlpredict1,
#                                   SurvivalProba = Survpredict1,
#                                   Size1 = exp(Growthpredict1),
#                                   FloweringProba = Flowpredict1)
# 
# save(predi_data, file = "Prediction")




