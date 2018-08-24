#clean routine
rm (list=ls())
setwd("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/PNAS_RESUB")
####Package Installer####
#delete "#" if needed#
#install.packages("ggplot2")
#install.packages("ggthemes")
#install.packages("plyr")
#install.packages("gridExtra")
#install.packages("betapart")
#install.packages("latticeExtra")
#install.packages("sjstats")

####load installed libraries####
library(ggplot2)
library(ggthemes)
library(plyr)
library(gridExtra)
library(latticeExtra)
library(fitdistrplus)
library(sjstats)
library(reshape2)

theme_set(theme_classic()+  theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                  strip.background = element_rect(colour ="black", fill="grey90"), 
                                  strip.text.x = element_text(face = "bold"), 
                                  panel.grid.minor = element_blank(),  
                                  panel.grid.major = element_blank()))

####load the Data####
Data <-  read.csv("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/PNAS_RESUB/newdiversities.csv", 
                  row.names = NULL, sep=";")

#add factor Period which is important for later grouping
attach(Data)
Data$Period[Age < 550 & Age > 485.4] <- "Cambrian"
Data$Period[Age < 485.4 & Age > 443.8] <- "Ordovician"
Data$Period[Age < 443.8 & Age > 419.2] <- "Silurian"
Data$Period[Age < 419.2 & Age > 358.9] <- "Devonian"
Data$Period[Age < 358.9 & Age > 298] <- "Carboniferous"
Data$Period[Age < 298.9 & Age > 252.17] <- "Permian"
Data$Period[Age < 252.17 & Age > 201.3] <- "Triassic"
Data$Period[Age < 201.3 & Age > 145.0] <- "Jurassic"
Data$Period[Age < 145.0 & Age > 66.0] <- "Cretaceous"
Data$Period[Age < 66.0 & Age > 23.03] <- "Paleogene"
Data$Period[Age < 23.03 & Age > 0] <- "Neogene"
detach(Data)
#make it factor
Data$Period.new = factor(Data$Period, levels=c("Cambrian","Ordovician","Silurian","Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene"))

####References per Formation####
Ref.Rho <- c()
Ref.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$RefForm, method="spearman", exact=FALSE)
Ref.Rho <- c(Ref.Rho, tmp$estimate)
Ref.pval <- c(Ref.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$RefForm, method="spearman", exact=FALSE)
Ref.Rho <- c(Ref.Rho, tmp$estimate)
Ref.pval <- c(Ref.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$RefForm, method="spearman", exact=FALSE)
Ref.Rho <- c(Ref.Rho, tmp$estimate)
Ref.pval <- c(Ref.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$RefForm, method="spearman", exact=FALSE)
Ref.Rho <- c(Ref.Rho, tmp$estimate)
Ref.pval <- c(Ref.pval, tmp$p.value)

Ref.infl <- cbind(Ref.Rho, Ref.pval)
rownames(Ref.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")



####Maximum Great Circle Distance####
MaxGCD.Rho <- c()
MaxGCD.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$maxGCD, method="spearman", exact=FALSE)
MaxGCD.Rho <- c(MaxGCD.Rho, tmp$estimate)
MaxGCD.pval <- c(MaxGCD.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$maxGCD, method="spearman", exact=FALSE)
MaxGCD.Rho <- c(MaxGCD.Rho, tmp$estimate)
MaxGCD.pval <- c(MaxGCD.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$maxGCD, method="spearman", exact=FALSE)
MaxGCD.Rho <- c(MaxGCD.Rho, tmp$estimate)
MaxGCD.pval <- c(MaxGCD.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$maxGCD, method="spearman", exact=FALSE)
MaxGCD.Rho <- c(MaxGCD.Rho, tmp$estimate)
MaxGCD.pval <- c(MaxGCD.pval, tmp$p.value)

MaxGCD.infl <- cbind(MaxGCD.Rho, MaxGCD.pval)
rownames(MaxGCD.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")
#write.csv(MaxGCD.infl, file = "Rho_MaxGCD")


####Median Great Circle Distance####
MedGCD.Rho <- c()
MedGCD.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$medianGCD, method="spearman", exact=FALSE)
MedGCD.Rho <- c(MedGCD.Rho, tmp$estimate)
MedGCD.pval <- c(MedGCD.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$medianGCD, method="spearman", exact=FALSE)
MedGCD.Rho <- c(MedGCD.Rho, tmp$estimate)
MedGCD.pval <- c(MedGCD.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$medianGCD, method="spearman", exact=FALSE)
MedGCD.Rho <- c(MedGCD.Rho, tmp$estimate)
MedGCD.pval <- c(MedGCD.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$medianGCD, method="spearman", exact=FALSE)
MedGCD.Rho <- c(MedGCD.Rho, tmp$estimate)
MedGCD.pval <- c(MedGCD.pval, tmp$p.value)

MedGCD.infl <- cbind(MedGCD.Rho, MedGCD.pval)
rownames(MedGCD.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")


####Duration####
Dur.Rho <- c()
Dur.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$Duration, method="spearman", exact=FALSE)
Dur.Rho <- c(Dur.Rho, tmp$estimate)
Dur.pval <- c(Dur.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$Duration, method="spearman", exact=FALSE)
Dur.Rho <- c(Dur.Rho, tmp$estimate)
Dur.pval <- c(Dur.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$Duration, method="spearman", exact=FALSE)
Dur.Rho <- c(Dur.Rho, tmp$estimate)
Dur.pval <- c(Dur.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$Duration, method="spearman", exact=FALSE)
Dur.Rho <- c(Dur.Rho, tmp$estimate)
Dur.pval <- c(Dur.pval, tmp$p.value)

Dur.infl <- cbind(Dur.Rho, Dur.pval)
rownames(Dur.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")


####NumberofCollections####
Coll.Rho <- c()
Coll.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$CollpF, method="spearman", exact=FALSE)
Coll.Rho <- c(Coll.Rho, tmp$estimate)
Coll.pval <- c(Coll.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$CollpF, method="spearman", exact=FALSE)
Coll.Rho <- c(Coll.Rho, tmp$estimate)
Coll.pval <- c(Coll.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$CollpF, method="spearman", exact=FALSE)
Coll.Rho <- c(Coll.Rho, tmp$estimate)
Coll.pval <- c(Coll.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$CollpF, method="spearman", exact=FALSE)
Coll.Rho <- c(Coll.Rho, tmp$estimate)
Coll.pval <- c(Coll.pval, tmp$p.value)

Coll.infl <- cbind(Coll.Rho, Coll.pval)
rownames(Coll.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")

####NumberofEnvironments####
Env.Rho <- c()
Env.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$Environments, method="spearman", exact=FALSE)
Env.Rho <- c(Env.Rho, tmp$estimate)
Env.pval <- c(Env.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$Environments, method="spearman", exact=FALSE)
Env.Rho <- c(Env.Rho, tmp$estimate)
Env.pval <- c(Env.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$Environments, method="spearman", exact=FALSE)
Env.Rho <- c(Env.Rho, tmp$estimate)
Env.pval <- c(Env.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$Environments, method="spearman", exact=FALSE)
Env.Rho <- c(Env.Rho, tmp$estimate)
Env.pval <- c(Env.pval, tmp$p.value)

Env.infl <- cbind(Env.Rho, Env.pval)
rownames(Env.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")

####madGCD####
madGCD.Rho <- c()
madGCD.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$madGCD, method="spearman", exact=FALSE)
madGCD.Rho <- c(madGCD.Rho, tmp$estimate)
madGCD.pval <- c(madGCD.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$madGCD, method="spearman", exact=FALSE)
madGCD.Rho <- c(madGCD.Rho, tmp$estimate)
madGCD.pval <- c(madGCD.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$madGCD, method="spearman", exact=FALSE)
madGCD.Rho <- c(madGCD.Rho, tmp$estimate)
madGCD.pval <- c(madGCD.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$madGCD, method="spearman", exact=FALSE)
madGCD.Rho <- c(madGCD.Rho, tmp$estimate)
madGCD.pval <- c(madGCD.pval, tmp$p.value)

madGCD.infl <- cbind(madGCD.Rho, madGCD.pval)
rownames(madGCD.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")

####Htaxa####
Htaxa.Rho <- c()
Htaxa.pval <- c()

tmp <- cor.test(Data$GammaForm,Data$Htaxa, method="spearman", exact=FALSE)
Htaxa.Rho <- c(Htaxa.Rho, tmp$estimate)
Htaxa.pval <- c(Htaxa.pval, tmp$p.value)
tmp <- cor.test(Data$AlphaForm,Data$Htaxa, method="spearman", exact=FALSE)
Htaxa.Rho <- c(Htaxa.Rho, tmp$estimate)
Htaxa.pval <- c(Htaxa.pval, tmp$p.value)
tmp <- cor.test(Data$BetaWForm,Data$Htaxa, method="spearman", exact=FALSE)
Htaxa.Rho <- c(Htaxa.Rho, tmp$estimate)
Htaxa.pval <- c(Htaxa.pval, tmp$p.value)
tmp <- cor.test(Data$BetaSimForm,Data$Htaxa, method="spearman", exact=FALSE)
Htaxa.Rho <- c(Htaxa.Rho, tmp$estimate)
Htaxa.pval <- c(Htaxa.pval, tmp$p.value)

Htaxa.infl <- cbind(Htaxa.Rho, Htaxa.pval)
rownames(Htaxa.infl) <- c("Gamma", "Alpha", "BetaW", "BetaSim")



####putting shit together####
all.rho <- cbind(Coll.infl, Ref.infl, Dur.infl, Env.infl, Htaxa.infl, MaxGCD.infl, MedGCD.infl, madGCD.infl) 
write.csv(all.rho, file = "Rhos")



DataRef10 <- subset(Data, RefForm < 10)

Data2 <- subset(Data, madGCD < 450)
ggplot(Data2, aes(y=BetaSimForm , x=maxGCD))+
  geom_point()+
  geom_smooth(method='lm')
+
 
  facet_grid(Period.new ~ .)