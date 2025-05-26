library(spaMM)
library(tidyverse)

### DATA
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

cptldata <- centauree_data_complet[centauree_data_complet$Flowering0!=0,]

#Number of capitula model
Cptlglm1 <- fitme(log(Cptl0) ~ 1 + Size0Mars + (age0|year), 
                  data=cptldata)


###CALCUL DE ESTB

#Calcul avec le fs0 du modèle matriciel

load("matricesperyear.RData")
Myear <- matrices

load("/media/loic/Commun/0Travail/Stage 2025 ISEM/Bibliographie/dynamique_des_pops/E1matrices.Rdata")
ME1 <- matrices
load("/media/loic/Commun/0Travail/Stage 2025 ISEM/Bibliographie/dynamique_des_pops/E2matrices.Rdata")
ME2 <- matrices
load("/media/loic/Commun/0Travail/Stage 2025 ISEM/Bibliographie/dynamique_des_pops/Crmatrices.Rdata")
MCr <- matrices
load("/media/loic/Commun/0Travail/Stage 2025 ISEM/Bibliographie/dynamique_des_pops/Aumatrices.Rdata")
MAu <- matrices
load("/media/loic/Commun/0Travail/Stage 2025 ISEM/Bibliographie/dynamique_des_pops/Pomatrices.Rdata")
MPo <- matrices
load("/media/loic/Commun/0Travail/Stage 2025 ISEM/Bibliographie/dynamique_des_pops/Pematrices.Rdata")
MPe <- matrices

annees <- 1995:2021
populations <- c("Po","Au","Pe","E1","E2","Cr")

fs0M <- expand.grid(Pop = populations,
                        year = 1995:2021,
                        fs0 = NA)

matrices <- list(MPo,MAu,MPe,ME1,ME2,MCr)

for(j in 1:6) {
  for (i in 1:27) {
    pop <- populations[j]
    year <- annees[i]
    fs0M[fs0M$year==year & fs0M$Pop==pop,]$fs0 <- matrices[[j]][[i]][1, 3]
  }
}

cptl_data_predi <- cptldata %>% 
  mutate(Cptl0 = ifelse(is.na(Cptl0), exp(2.31070+0.06846*Size0Mars), Cptl0))

cptl <- cptl_data_predi %>% 
  group_by(year,Pop) %>% 
  summarize(NombreCapitules = mean(Cptl0))

Estb <- cptl %>% 
  right_join(fs0M, by=join_by(year,Pop)) %>% 
  arrange(year) %>% 
  mutate(EstbRate = fs0/NombreCapitules)

Estbglm1 <- fitme(EstbRate ~ 1+(1|year), data=Estb)

mean(Estb$EstbRate,na.rm=TRUE)
#mean : 0.26


#Calcul avec la formule d'origine utilisant les données brutes

plt <- centauree_data_complet %>% 
  filter(age0==1) %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(NombrePlantules = sum(age0))

cptl <- cptl_data_predi %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(NombreCapitules = sum(Cptl0))

Estb <- inner_join(plt,cptl, by=join_by(Quadrat,year,Pop))

Estb <- Estb %>% mutate(EstbRate=rep(NA)) %>% 
  arrange(Quadrat)

for (i in 2:length(Estb$Quadrat)){
  if (Estb$Quadrat[i]!=Estb$Quadrat[i-1]){next}
  if (Estb$year[i]!=Estb$year[i-1]+1){next}
  Estb$EstbRate[i] <- Estb$NombrePlantules[i]/Estb$NombreCapitules[i-1]
}

Estb <- Estb %>%
  group_by(Pop,year) %>% 
  mutate(EstbPop = mean(EstbRate,na.rm=TRUE))

Estbglm1 <- fitme(EstbRate ~ 1 + (1|Pop:year), data=Estb)
