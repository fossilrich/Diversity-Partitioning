#clean routine
rm (list=ls())
#setwd("//...")


####Package Installer####
#delete "#" if needed#
#install.packages("ggplot2")
#install.packages("ggthemes")
#install.packages("plyr")
#install.packages("gridExtra")
#install.packages("betapart")
#install.packages("latticeExtra")
#install.packages("sjstats")
#install.packages("lme4")


####load installed libraries####
library(ggplot2)
library(ggthemes)
library(plyr)
library(gridExtra)
library(latticeExtra)
library(fitdistrplus)
library(sjstats)
library(reshape2)
library(lme4)

theme_set(theme_classic()+  theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                  strip.background = element_rect(colour ="black", fill="grey90"), 
                                  strip.text.x = element_text(face = "bold"), 
                                  panel.grid.minor = element_blank(),  
                                  panel.grid.major = element_blank()))

####Read the Data####
#local if downloaded to your computer
#Data <- read.csv("//diversities.csv", row.names = NULL, sep=";")
#web
Data <- read.csv(url("https://raw.githubusercontent.com/fossilrich/Diversity-Partitioning/master/diversities.csv"), row.names = NULL, sep=";")


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

####FIGURES MAIN TEXT####


####Figure 1: A-B-G-Plot#### 
ggplot(Data, aes(x=GammaForm, y=BetaSimForm)) +
  stat_smooth(method="loess", span = 1, size = 0.8, colour="darkblue", se = TRUE, alpha = 0.3)+
  geom_point(shape=19, size = 2, alpha = 0.2, colour = "darkblue")+
  stat_smooth(method="loess", span = 1, size = 0.8, colour="red3", se = TRUE, alpha = 0.3, aes(x=GammaForm, y=AlphaForm/40)) +
  geom_point(shape=19, size = 2, alpha = 0.2, colour = "red3", data = Data, aes(x=GammaForm, y=AlphaForm/40))+
  scale_x_continuous(name = "Mean gamma diversity per formation [species]", expand = c(0,1)) +
  scale_y_continuous(name = "Mean beta diversity per formation [Simpson's Metric]",  
                     sec.axis = sec_axis(~.*40, name = "Mean alpha diversity per formation [species]"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))
#ggsave(("Fig_1_Phan-abg.png"), plot = last_plot(), width=8.7, height=9, units = "cm", dpi = 300)
#ggsave(("Fig_1_Phan-abg.pdf"), plot = last_plot(), width=8.7, height=9, units = "cm")

####Figure 2: Periods A-B-G-Plot####
ggplot(Data, aes(x=GammaForm, y=AlphaForm)) +
  stat_smooth(method="lm", colour="red3", se = TRUE, alpha = 0.3, size = 0.8)+
  geom_point(shape=19, size = 2, alpha = 0.2, colour = "red3")+
  stat_smooth(method="lm", formula = y ~ log(x), colour="blue", size = 0.8, se = TRUE, alpha = 0.3, aes(x=GammaForm, y=BetaSimForm*50)) +
  geom_point(shape=19, size = 2, alpha = 0.2, colour = "darkblue", data = Data, aes(x=GammaForm, y=BetaSimForm*50))+
  scale_x_continuous(name = "Mean gamma diversity per formation [species]", expand = c(0,1)) +
  #scale_y_continuous(name = "Mean alpha (red) and beta (blue) diversity per formation", expand = c(0,1))+
  scale_y_continuous(name = "Mean alpha diversity (red) per formation [species]",  
                     sec.axis = sec_axis(~./50, name = "Mean beta diversity (blue) per formation [Simpson's Metric]" ))+
  #facet_wrap( ~ Period.new, scales = "free_x")+ #for adjusted scales
  facet_wrap( ~ Period.new)+ #with fixed scale
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8), 
        strip.text = element_text(size=9))
#ggsave(("Fig_2_Period-abg-fixed.jpg"), plot = last_plot(), width=18, height=17.5, units = "cm", dpi = 300)
#ggsave(("Fig_2_Period-abg-fixed.pdf"), plot = last_plot(),width=18, height=17.5, units = "cm")


####SUPPLEMENTARY MATERIAL####


####Figure S2 Time Series####
boundaries <- c(485.4, 443.8, 419.2, 358.9, 298.9, 252.17, 201.3, 145, 66, 23.03)
extinctions <- c(444, 372, 252, 201.5, 66)
alpha <- .4
Alphaplot <- ggplot(Data, aes(x=Age, y=AlphaForm)) + 
  geom_vline(xintercept=boundaries, colour="black", linetype="dashed")+
  geom_vline(xintercept=extinctions, colour="darkgrey", linetype="solid", alpha = 0.4, size = 3)+
  geom_point(colour = "red3", alpha = 0.3) +
  geom_smooth(span = 0.15, se=TRUE, colour ="red3",  alpha = alpha)+
  labs(title = "A")+
  scale_x_continuous(trans = "reverse", name = "Age (Ma)", position = "top", expand = c(0,3), 
                     breaks = c(50,100,150,200,250,300,350,400,450,500))+
  scale_y_continuous(name = "Alpha diversity [species]", expand = c(0,0.5))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank())

Gammaplot <- ggplot(Data, aes(x=Age, y=GammaForm)) + 
  geom_vline(xintercept=boundaries, colour="black", linetype="dashed")+
  geom_vline(xintercept=extinctions, colour="darkgrey", linetype="solid", alpha = 0.4, size = 3)+
  geom_point(colour = "orange", alpha = alpha) +
  geom_smooth(span = 0.15, se=TRUE, colour ="orange", alpha = alpha)+
  labs(title = "B")+
  scale_x_continuous(trans = "reverse", expand = c(0,3))+
  #scale_x_continuous(name = "Period", trans = "reverse",breaks = c(513, 465, 431.5, 389, 329, 275, 227, 173.3, 100, 44.5, 2),#these are no exact ages, they just make the labels appear in the middle of each period
  #                   labels = c("Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic",
  #                              "Jurassic", "Cretaceous", "Palaeogene", "Neogene"))+
  scale_y_continuous(name = "Gamma diversity [species]", expand = c(0,0.5))+
  theme(axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank())

Betaplot <-  ggplot(Data, aes(x=Age, y=BetaSimForm)) + 
  geom_vline(xintercept=boundaries, colour="black", linetype="dashed")+
  geom_vline(xintercept=extinctions, colour="darkgrey", linetype="solid", alpha = 0.4, size = 3)+
  geom_smooth(span = 0.15, se=TRUE, colour ="blue", alpha = alpha)+
  geom_point(colour = "darkblue", alpha = alpha) +
  labs(title = "C")+
  scale_x_continuous(trans = "reverse", expand = c(0,3))+
  #scale_x_continuous(name = "Period", trans = "reverse", expand = c(0,3),
  #                   breaks = c(502, 465, 431.5, 389, 329, 275, 227, 173.3, 105, 44.5, 8),#these are no exact ages, they just make the labels appear in the middle of each period
  #                   labels = c("Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous",
  #                            "Permian", "Triassic", "Jurassic", "Cretaceous", "Palaeogene", "Neogene"))+
  scale_y_continuous(name = "Beta diversity [Simpson's Metric]", expand = c(0,0.05), limits = c(NA,NA))+
  theme(axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank())

Betaplot2 <-  ggplot(Data, aes(x=Age, y=BetaWForm)) + 
  geom_vline(xintercept=boundaries, colour="black", linetype="dashed")+
  geom_vline(xintercept=extinctions, colour="darkgrey", linetype="solid", alpha = 0.4, size = 3)+
  geom_smooth(span = 0.15, se=TRUE, colour ="blue", alpha = alpha)+
  geom_point(colour = "darkblue", alpha = alpha) +
  labs(title = "D")+
  scale_x_continuous(name = "Period", trans = "reverse", expand = c(0,3),
                     breaks = c(502, 465, 431.5, 389, 329, 275, 227, 173.3, 105, 44.5, 8),#these are no exact ages, they just make the labels appear in the middle of each period
                     labels = c("Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous",
                                "Permian", "Triassic", "Jurassic", "Cretaceous", "Palaeogene", "Neogene"))+
  scale_y_continuous(name = "Beta Diversity [Whittaker's Beta]")+#, expand = c(0,0.05))+
  theme(axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(), 
        panel.ontop = FALSE,
        axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(Alphaplot, Gammaplot, Betaplot, Betaplot2, ncol=1, left="Mean diversity per formation")
#ggsave("Fig_S2_time-series.jpg", width=12, height=28, dpi = 300, units = "cm", grid.arrange(Alphaplot, Gammaplot, Betaplot,Betaplot2, ncol=1, left="Mean diversity per formation"))
#ggsave("Fig_S2_time-series.pdf", width=12, height=28, units = "cm", grid.arrange(Alphaplot, Gammaplot, Betaplot,Betaplot2, ncol=1, left="Mean diversity per formation"))

####Figure S3 Environments per Formation####
ggplot(Data, aes(x=Environments))+
  geom_bar()+
  scale_x_continuous(name = "Numbers of Environments recorded") +
  scale_y_continuous(name = "Number of formations")+
  facet_wrap(~Period.new)
#ggsave(("Fig_S3_Environments.jpg"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S3_Environments.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

####Figure S4: all Periods in one B-G-Plot####
ggplot(Data, aes(x=GammaForm, y=BetaSimForm), colour = factor(Period.new))+
  #geom_point(aes(colour = factor(Period.new)), alpha = 0.1, size = 3, shape = 19)+
  scale_x_continuous(name = "Gamma Diversity [species]")+
  scale_y_continuous(name = "Beta Diversity [Simpson's Metric]")+
  stat_smooth(method="lm", formula = y ~ log(x), se = FALSE, alpha = 0.2, aes(x=GammaForm, y=BetaSimForm, colour=factor(Period.new)))+
  scale_colour_hue("Period", h = c(20, 340))
#ggsave(("Fig_S4_beta_allinone.jpg"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S4_beta_allinone.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

####Figure S5: Average Beta####
means.BS <- ddply(Data, "Period.new", summarise, rating.mean=mean(BetaSimForm))
ggplot(Data, aes(x=BetaSimForm)) +
  scale_x_continuous(name = "Beta diversity [Simpson's Metric]")+
  scale_y_continuous(name = "Density (frequency distribution)")+
  geom_density(alpha= 0.4, fill = "darkblue", colour = NA)+
  geom_vline(data=means.BS, aes(xintercept=rating.mean),
             linetype="dashed", size=0.3, colour = "darkblue")  +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  panel.grid.major = element_blank(), 
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.y = element_text(size = 7))+
  #element_text(face = "bold"))+
  facet_grid(Period.new ~ .)
  #ggsave(("Fig_S5_average-beta.jpg"), plot = last_plot(), width=12, height=22, units = "cm", dpi = 300)
  #ggsave(("Fig_S5_average-beta.pdf"), plot = last_plot(), width=12, height=22, units = "cm")


####NULL MODEL TESTING####
####Building the Null model for a Gamma < 100 specices####
Data2 <- subset(Data, Data$Gamma < 100) ####Because this is Threshold from which beta levels off anyway

set.seed(666) # sampling is never random. make it reproducible!
n <- 200 # number of repetitions
sampsize <- 20 # sample size
#output <- data.frame(matrix(nrow=sampsize, ncol=ncol(Data2)))
for(i in 1:n){
  subData2 <- Data2[sample(1:nrow(Data2), sampsize),]
  if(i==1){
    output.2 <- subData2
  }else{output.2 <- rbind(output.2, subData2)}
}
repetition.2 <- rep(c(1:n), each=sampsize)
output.2$repetition <- as.factor(repetition)
sub.output.2 <- output.2[1:(sampsize*20),]

####Building Beta slope functions for BetaSim for Gamma < 100 species ####
slope_beta.2 <- function(xx){
  lm1.2 <- lm(xx$BetaSimForm~xx$GammaForm)
  pval.2 <- summary(lm(xx$BetaSimForm~xx$GammaForm))$coefficients[2,4]
  # add CI for the slope
  SE <- summary(lm1.2)[["coefficients"]][2,2] 
  upper = lm1.2$coefficients[[2]] + 1.96*SE
  lower = lm1.2$coefficients[[2]] - 1.96*SE
  return(data.frame(slope = lm1.2$coefficients[[2]], pvalue = pval.2, lowerCI=lower, upperCI=upper))
}
df.slope.b.2 <- ddply(Data2, .(Period.new), slope_beta.2)
output.slope_beta.2 <- ddply(output.2, .(repetition), slope_beta.2) # null model beta slope
df.slope.b.2$signif <- df.slope.b.2$pvalue<0.05

####Figure S6: Global Null model####
ggplot(output.2, aes(x=GammaForm, y=AlphaForm, group=repetition)) +
  geom_line(stat="smooth",method = "lm", alpha = 0.1, colour="red3")+
  #geom_point(shape=19, size = 2, alpha = 0.4, colour = "darkblue")+
  geom_line(stat="smooth",method = "lm", formula = y ~ log(x), colour="blue", se = FALSE, alpha = 0.2, aes(x=GammaForm, y=BetaSimForm*45))+
  #geom_point(shape=19, size = 2, alpha = 0.4, colour = "red3", data = output, aes(x=GammaForm, y=AlphaForm/40))+
  scale_x_continuous(name = "Mean gamma diversity [species] per random sample")+#,  expand = c(0,1)) +
  scale_y_continuous(name = "Mean alpha (red) diversity [species] per random sample",
                     sec.axis = sec_axis(~./45, name = "Mean beta (blue) diversity [Simpson's Metric] per random sample"))
# annotate("text", x=300, y=0.1, label=paste("repetitions of 20 random samplings: ", n))
#ggsave(("Fig_S6_global-nullmodel.jpg"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S6_global-nullmodel.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

####Figure S7: Beta Slopes####
ggplot(df.slope.b.2, aes(x=Period.new, y=slope))+
  geom_point(size=3)+
  geom_line(aes(x=1:11,y=df.slope.b.2$slope), colour = "black", size=0.7)+
  geom_errorbar(aes(x=Period.new, ymin=lowerCI, ymax=upperCI), width=0.1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(aes(yintercept = mean(output.slope_beta.2$slope)), colour = "red")+
  scale_x_discrete("Period")+
  scale_y_continuous("Slope")+
  geom_hline(aes(yintercept = mean(output.slope_beta.2$slope)+
                   2*sd(output.slope_beta.2$slope)), color="red", lty=2)+
  geom_hline(aes(yintercept = mean(output.slope_beta.2$slope)-
                   2*sd(output.slope_beta.2$slope)), color="red", lty=2)
#ggsave(("Fig_S7_Beta-slopes.jpg"), plot = last_plot(), width=12, height=12, dpi = 300, units = "cm")
#ggsave(("Fig_S7_Beta-slopes.pdf"), plot = last_plot(), width=12, height=12, units = "cm")


####Figure S8: Higher Taxa present####
occtax <- read.csv("//occtax.csv", 
                   row.names = NULL, sep=",")

# Get number of higher taxa in each collection (occtax):
table(occtax$collection_no %in% df$collectionID) # select needed collections
occtax <- occtax[occtax$collection_no%in%df$collectionID,]
## check class names:
unique(occtax$phylum)
occtax$my_taxnomy <- as.character(occtax$phylum)
#occtax$my_taxnomy[which(occtax$class=="Bivalvia")] <- "Bivalvia"
#occtax$my_taxnomy[which(occtax$class=="Gastropoda")] <- "Gastropoda"
## get the numbers
temp <- tapply(occtax$my_taxnomy, occtax$collection_no, function(x)length(unique(x)) )
ht_df <- as.data.frame(temp)
library(plyr)
ht_df <- adply(temp, c(1))
names(ht_df) <- c("collection_no", "higher_taxa_no")
occtax$mid_ma <- (occtax$max_ma + occtax$min_ma)/2
## assign periods to the collections:
occtax$Period[occtax$mid_ma < 550 & occtax$mid_ma > 485.4] <- "Cambrian"
occtax$Period[occtax$mid_ma < 485.4 & occtax$mid_ma > 443.4] <- "Ordovician"
occtax$Period[occtax$mid_ma < 443.4 & occtax$mid_ma > 419.2] <- "Silurian"
occtax$Period[occtax$mid_ma < 419.2 & occtax$mid_ma > 358.9] <- "Devonian"
occtax$Period[occtax$mid_ma < 358.9 & occtax$mid_ma > 298.9] <- "Carboniferous"
occtax$Period[occtax$mid_ma < 298.9 & occtax$mid_ma > 252.2] <- "Permian"
occtax$Period[occtax$mid_ma < 252.2 & occtax$mid_ma > 201.3] <- "Triassic"
occtax$Period[occtax$mid_ma < 201.3 & occtax$mid_ma > 145.0] <- "Jurassic"
occtax$Period[occtax$mid_ma < 145.0 & occtax$mid_ma > 66.0] <- "Cretaceous"
occtax$Period[occtax$mid_ma < 66.0 & occtax$mid_ma > 23.03] <- "Paleogene"
occtax$Period[occtax$mid_ma < 23.03 & occtax$mid_ma > 0] <- "Neogene"

table(occtax$Period)
time <- subset(occtax, select = c("collection_no", "Period"))
time <- unique(time)
ht_df <- merge(ht_df, time, all=TRUE)
ht_df_new <- unique(ht_df)
table(ht_df_new$Period)
periods <- unique(ht_df_new$Period)
results <- matrix(ncol = length(periods), nrow = 600)
for(i in 1:length(periods)) {
  temp <- ht_df_new[ht_df_new$Period == periods[i],]
  #print(nrow(temp))
  temptwo <- temp[sample(nrow(temp), 600, replace = FALSE), ] 
  results[,i] <- temptwo$higher_taxa_no
}
results <- as.data.frame(results)
colnames(results) <- periods
results.final<- melt(results) 
results.final$variable = factor(results.final$variable, levels=c("Cambrian","Ordovician","Silurian","Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "not.assigned"))
results.final <- results.final
colnames(results.final) <- c("Period", "Taxa") 

results.final <- subset(results.final, results.final$Period != "NA")
ggplot(results.final, aes(x=Taxa))+
  geom_bar()+
  scale_x_continuous(name = "numbers of higher taxa present") +
  scale_y_continuous(name = "number of collections")+
  facet_wrap(~Period)
#ggsave(("Fig_S8_taxa-present.jpg"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S8_taxa-present.pdf"), plot = last_plot(), width=20, height=19, units = "cm")




####Figure S9: Full data vs. Subset####
#make the subset

DataRef10 <- subset(Data, RefForm < 10)

ggplot(Data, aes(x=GammaForm, y=AlphaForm)) +
  stat_smooth(method="loess", span = 1, colour="red3", se = TRUE, alpha = 0.3)+
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "red3")+
  stat_smooth(method="loess", span = 1, colour="blue", se = TRUE, alpha = 0.3, aes(x=GammaForm, y=BetaSimForm*45)) +
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "darkblue", data = Data, aes(x=GammaForm, y=BetaSimForm*45))+
  #subsampled points
  stat_smooth(method="loess", span = 1, lty="longdash" , colour="darkgrey", se = TRUE, alpha = 0.4, data = DataRef10, aes(x=GammaForm, y=BetaSimForm*45))+
  geom_point(shape=21, size = 3, colour = "black", data = DataRef10, aes(x=GammaForm, y=BetaSimForm*45))+
  stat_smooth(method="loess", span = 1,lty="longdash" , colour="darkgrey", se = TRUE, alpha = 0.4, data = DataRef10, aes(x=GammaForm, y=AlphaForm)) +
  geom_point(shape=21, size = 3, colour = "black", data = DataRef10, aes(x=GammaForm, y=AlphaForm))+
  scale_x_continuous(name = "Mean gamma diversity per formation [species]", expand = c(0,1)) +
  #scale_y_continuous(name = "Mean alpha (red) and beta (blue) diversity per formation", expand = c(0,1))+
  scale_y_continuous(name = "Mean alpha diversity (red) per formation [species]",  
                     sec.axis = sec_axis(~./45, name = "Mean beta diversity (blue) per formation [Simpson's Metric]" ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank())
#ggsave(("Fig_S09_subset_vs_full.jpg"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S09_subset_vs_full.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

####Figure S10: Full data vs. Subset Periods####
ggplot(Data, aes(x=GammaForm, y=AlphaForm)) +
  stat_smooth(method="lm", colour="red3", se = TRUE, alpha = 0.3)+
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "red3")+
  stat_smooth(method="lm",formula = y ~ log(x), colour="blue", se = TRUE, alpha = 0.3, aes(x=GammaForm, y=BetaSimForm*45)) +
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "darkblue", data = Data, aes(x=GammaForm, y=BetaSimForm*45))+
  #subsampled points
  stat_smooth(method="lm", lty="longdash" , colour="darkgrey", se = TRUE, alpha = 0.4, data = DataRef10, aes(x=GammaForm, y=BetaSimForm*45))+
  geom_point(shape=21, size = 3, colour = "black", data = DataRef10, aes(x=GammaForm, y=BetaSimForm*45))+
  stat_smooth(method="lm", lty="longdash" , colour="darkgrey", se = TRUE, alpha = 0.4, data = DataRef10, aes(x=GammaForm, y=AlphaForm)) +
  geom_point(shape=21, size = 3, colour = "black", data = DataRef10, aes(x=GammaForm, y=AlphaForm))+
  scale_x_continuous(name = "Mean gamma diversity per formation [species]", expand = c(0,1)) +
  #scale_y_continuous(name = "Mean alpha (red) and beta (blue) diversity per formation", expand = c(0,1))+
  scale_y_continuous(name = "Mean alpha diversity (red) per formation [species]",  
                     sec.axis = sec_axis(~./45, name = "Mean beta diversity (blue) per formation [Simpson's Metric]" ))+
  facet_wrap( ~ Period.new)+ #with fixed scale
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank())
#ggsave(("Fig_S10_subset_vs_full-periods.jpg"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S10_subset_vs_full-periods.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

##### Figure S11: Goodness of fit for different functions####
r.squared <- c()
r.squared.log <- c()
r.squared.poly3 <- c()
r.squared.beta <- c()
r.squared.log.beta <- c()
r.squared.poly3.beta <- c()
for(i in 1:nlevels(Data$Period.new)){
  temp <- Data[Data$Period.new==levels(Data$Period.new)[i],]
  
  # linear model Alpha
  lm.temp <- lm(temp$AlphaForm ~ temp$GammaForm) 
  r.temp <- summary(lm.temp)$adj.r.squared
  r.squared <- c(r.squared, r.temp)
    # log model Alpha 
  log.temp <- lm(temp$AlphaForm ~ log(temp$GammaForm)) 
  r.temp <- summary(log.temp)$adj.r.squared
  r.squared.log <- c(r.squared.log, r.temp)
    # poly3 model Alpha 
  log.temp <- lm(temp$AlphaForm ~ poly(temp$GammaForm, 3)) 
  r.temp <- summary(log.temp)$adj.r.squared
  r.squared.poly3 <- c(r.squared.poly3, r.temp)
    # linear model Beta
  lm.temp <- lm(temp$BetaSimForm ~ temp$GammaForm) 
  r.temp <- summary(lm.temp)$adj.r.squared
  r.squared.beta <- c(r.squared.beta, r.temp)
   # log model Beta
  log.temp <- lm(temp$BetaSimForm ~ log(temp$GammaForm)) 
  r.temp <- summary(log.temp)$adj.r.squared
  r.squared.log.beta <- c(r.squared.log.beta, r.temp)
   # poly3 model Beta
  log.temp <- lm(temp$BetaSimForm ~ poly(temp$GammaForm, 3)) 
  r.temp <- summary(log.temp)$adj.r.squared
  r.squared.poly3.beta <- c(r.squared.poly3.beta, r.temp)
}
rdf <- data.frame(r.squared, r.squared.log, r.squared.poly3, r.squared.beta,
                  r.squared.log.beta, r.squared.poly3.beta)
rdf2 <- melt(rdf)
rdf2$Period <- factor(rep(levels(Data$Period.new), ncol(rdf)))
rdf2$ab <- c(rep("Alpha", nrow(rdf2)/2), rep("BetaSim", nrow(rdf2)/2))
rdf2$fit <- c(c(rep("linear", nrow(rdf2)/6), rep("log", nrow(rdf2)/6), rep("poly 3", nrow(rdf2)/6)),
              c(rep("linear", nrow(rdf2)/6), rep("log", nrow(rdf2)/6), rep("poly 3", nrow(rdf2)/6)))
rdf2$Period <- factor(rdf2$Period, levels=c("Cambrian","Ordovician","Silurian","Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene"))

ggplot(rdf2, aes(x=Period, y=value, group=variable, col=ab))+
  geom_line(aes(linetype=fit))+
  scale_linetype_discrete("Fit function")+
  scale_color_manual("Diversity", values=c("red3", "blue"))+
  scale_y_continuous("Adjusted R²")+
  scale_x_discrete("Period")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white", colour = "black"), 
        strip.background = element_rect(colour ="black", fill="grey90"), 
        strip.text.x = element_text(face = "bold"), 
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank())
#ggsave(("Fig_S11_goodness-of-fit.jpg"), plot = last_plot(), width=14, height=13, units = "cm", dpi = 300)
#ggsave(("Fig_S11_goodness-of-fit.pdf"), plot = last_plot(), width=14, height=13, units = "cm")

