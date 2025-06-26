library(spaMM)
library(tidyverse)

setwd("/media/loic/Commun/0Travail/Stage 2025 ISEM/Code/IPM")

IPM_data <- read.csv("newdata.csv")

IPM_data$Age[IPM_data$Age > 8] <- 8
centauree_data <- IPM_data[!is.na(IPM_data$Size0Mars) & !is.na(IPM_data$Age),]

spaMM.options(separation_max=70)

#Supprimer plantes dont l'age est inconnu
# centauree_data <- centauree_data[!is.na(centauree_data$Age), ]

#Forcer l'age maximal à 8
length(centauree_data$Age[centauree_data$Age >= 8])
centauree_data$Age[centauree_data$Age > 8] <- 8


#Number of capitula model
cptldata <- centauree_data[centauree_data$Flowering!=0,]

Cptlglm1 <- fitme(log(Capitule) ~ 1 + Size0Mars + (Age|year), 
                  data=cptldata)
# 2.32734 + 0.06752x

cptl_data_predi <- cptldata %>% 
  mutate(Capitule = ifelse(is.na(Capitule), exp(2.32734+0.06752*Size0Mars), Capitule))

###CALCUL DE ESTB

#Calcul avec la formule d'origine utilisant les données brutes

plt <- centauree_data_complet %>% 
  filter(Age==1) %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(NombrePlantules = sum(Age))

cptl <- cptl_data_predi %>% 
  group_by(Quadrat,year,Pop) %>% 
  summarize(NombreCapitules = sum(Capitule))

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


#Calcul en utilisant le nombre de plantules comme variable
annees <- 1995:2022
populations <- c("E2", "E1", "Au", "Po", "Pe", "Cr")
taille_range <- seq(0.5, 25, by = 0.5)
cap_range <- seq(1, 200)

fake_data <- expand.grid(
  year = annees,
  Pop = populations,
  Size0Mars = taille_range,
  Capitule = cap_range
)

fake_data <- fake_data %>%
  mutate(Nrw = row_number())

plt <- IPM_data %>%
  filter(Age == 1) %>%
  group_by(Quadrat, year, Pop) %>%
  summarize(NombrePlantules = sum(Age)) %>%
  mutate(year = year - 1)

cptl <- cptl_data_predi %>%
  group_by(Quadrat, year, Pop) %>%
  summarize(Capitule = sum(Capitule))

Estb <- left_join(plt, cptl, by = join_by(Quadrat, year, Pop))

Estbglm1 <- fitme(NombrePlantules ~ 1 + offset(log(Capitule))+(1|Pop:year), 
                  data=Estb,
                  family = Poisson(log))


predi <- predict(Estbglm1,newdata = fake_data)

Estbpredi <- fake_data %>%
  mutate(Plantule = predi) %>%
  group_by(year, Pop) %>%
  mutate(Plantule = mean(Plantule)) %>%
  select(!c(Size0Mars, Nrw)) %>%
  summarize(Estb = mean(Plantule / Capitule))

summary(Estbpredi)

Estbpredi %>% 
ggplot(aes(x = year, y = Estb)) +
  geom_line(aes(color = as.factor(Pop)), size = 0.75) +
  theme_bw() +
  scale_color_viridis_d(option = "plasma")

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
  mutate(Capitule = ifelse(is.na(Capitule), exp(2.31070+0.06846*Size0Mars), Capitule))

cptl <- cptl_data_predi %>% 
  group_by(year,Pop) %>% 
  summarize(NombreCapitules = mean(Capitule))

Estb <- cptl %>% 
  right_join(fs0M, by=join_by(year,Pop)) %>% 
  arrange(year) %>% 
  mutate(EstbRate = fs0/NombreCapitules)

Estbglm1 <- fitme(EstbRate ~ 1+(1|year), data=Estb)

mean(Estb$EstbRate,na.rm=TRUE)
#mean : 0.26

