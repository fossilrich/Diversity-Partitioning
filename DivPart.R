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

Data <-  read.csv("diversities.csv", row.names = NULL, sep=";")

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
  stat_smooth(method="loess", span = 1, colour="darkblue", se = TRUE, alpha = 0.3)+
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "darkblue")+
  stat_smooth(method="loess", span = 1, colour="red3", se = TRUE, alpha = 0.3, aes(x=GammaForm, y=AlphaForm/40)) +
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "red3", data = Data, aes(x=GammaForm, y=AlphaForm/40))+
  scale_x_continuous(name = "Mean gamma diversity per formation [species]", expand = c(0,1)) +
  scale_y_continuous(name = "Mean beta diversity per formation [Simpson's Metric]",  
                     sec.axis = sec_axis(~.*40, name = "Mean alpha diversity per formation [species]"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank())
#ggsave(("Fig_1_Phan-abg.png"), plot = last_plot(), width=12, height=12, units = "cm", dpi = 300)
#ggsave(("Fig_1_Phan-abg.jpg"), plot = last_plot(), width=12, height=12, units = "cm", dpi = 300)
#ggsave(("Fig_1_Phan-abg.pdf"), plot = last_plot(), width=12, height=12, units = "cm")

####Figure 2: Periods A-B-G-Plot####
ggplot(Data, aes(x=GammaForm, y=AlphaForm)) +
  stat_smooth(method="lm", colour="red3", se = TRUE, alpha = 0.3)+
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "red3")+
  stat_smooth(method="lm", formula = y ~ log(x), colour="blue", se = TRUE, alpha = 0.3, aes(x=GammaForm, y=BetaSimForm*50)) +
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "darkblue", data = Data, aes(x=GammaForm, y=BetaSimForm*50))+
  scale_x_continuous(name = "Mean gamma diversity per formation [species]", expand = c(0,1)) +
  #scale_y_continuous(name = "Mean alpha (red) and beta (blue) diversity per formation", expand = c(0,1))+
  scale_y_continuous(name = "Mean alpha diversity (red) per formation [species]",  
                     sec.axis = sec_axis(~./50, name = "Mean beta diversity (blue) per formation [Simpson's Metric]" ))+
  #facet_wrap( ~ Period.new, scales = "free_x")+ #for adjusted scales
  facet_wrap( ~ Period.new)+ #with fixed scale
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank())
#ggsave(("Fig_2_Period-abg-fixed.png"), plot = last_plot(), width=18, height=17.5, units = "cm", dpi = 300)
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

#ggsave("Fig_S2_time-series.png", width=12, height=28, dpi = 300, units = "cm", grid.arrange(Alphaplot, Gammaplot, Betaplot,Betaplot2, ncol=1, left="Mean diversity per formation"))
#ggsave("Fig_S2_time-series.pdf", width=12, height=28, units = "cm", grid.arrange(Alphaplot, Gammaplot, Betaplot,Betaplot2, ncol=1, left="Mean diversity per formation"))

####Figure S3: all Periods in one A-B-G-Plot####
ggplot(Data, aes(x=GammaForm, y=BetaSimForm), colour = factor(Period.new))+
  #geom_point(aes(colour = factor(Period.new)), alpha = 0.1, size = 3, shape = 19)+
  scale_x_continuous(name = "Gamma Diversity [species]")+
  scale_y_continuous(name = "Beta Diversity [Simpson's Metric]")+
  stat_smooth(method="lm", formula = y ~ log(x), se = FALSE, alpha = 0.2, aes(x=GammaForm, y=BetaSimForm, colour = factor(Period.new)))+
  scale_colour_hue("Period", h = c(160, 340))
#ggsave(("Fig_S3_beta_allinone.png"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S3_beta_allinone.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

####Figure S4: Average Beta####
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
  #ggsave(("Fig_S4_average-beta.png"), plot = last_plot(), width=12, height=22, units = "cm", dpi = 300)
  #ggsave(("Fig_S4_average-beta.pdf"), plot = last_plot(), width=12, height=22, units = "cm")


####NULL MODEL TESTING####
####Building the Null model####
set.seed(666) # sampling is never random. make it reproducible!
n <- 200 # number of repetitions
sampsize <- 20 # sample size
#output <- data.frame(matrix(nrow=sampsize, ncol=ncol(Data)))
for(i in 1:n){
  subData <- Data[sample(1:nrow(Data), sampsize),]
  if(i==1){
    output <- subData
  }else{output <- rbind(output, subData)}
}
repetition <- rep(c(1:n), each=sampsize)
output$repetition <- as.factor(repetition)
sub.output <- output[1:(sampsize*20),]

####Figure S5: Random Sample Plots####
ggplot(sub.output, aes(x=GammaForm, y=AlphaForm)) +
  stat_smooth(method="lm", colour="red3", se = TRUE, alpha = 0.3)+
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "red3")+
  stat_smooth(method="lm", formula = y ~ log(x), colour="darkblue", se = TRUE, alpha = 0.3, aes(x=GammaForm, y=BetaSimForm*50)) +
  geom_point(shape=19, size = 3, alpha = 0.2, colour = "darkblue", data = sub.output, aes(x=GammaForm, y=BetaSimForm*50))+
  scale_x_continuous(name = "Mean gamma diversity [species] per formation", expand = c(0,1)) +
  #scale_y_continuous(name = "mean alpha (red) and beta (blue) diversity per formation", expand = c(0,1))+
  scale_y_continuous(name = "Mean alpha diversity (red) per formation [species]",  
                     sec.axis = sec_axis(~./50, name = "Mean beta diversity (blue) per formation [Simpson's Metric]" ))+
  facet_wrap( ~ repetition)+
  #facet_wrap( ~ repetition, scales = "free")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
      panel.grid.minor = element_blank(),  
      panel.grid.major = element_blank())
#ggsave(("Fig_S5_Random-sample-plots.png"), plot = last_plot(), width=20, height=23, units = "cm", dpi = 300)
#ggsave(("Fig_S5_Random-sample-plots.pdf"), plot = last_plot(), width=20, height=23, units = "cm")

####Figure S6: Global Nullmodell####
ggplot(output, aes(x=GammaForm, y=AlphaForm, group=repetition)) +
  geom_line(stat="smooth",method = "lm", alpha = 0.1, colour="red3")+
  #geom_point(shape=19, size = 2, alpha = 0.4, colour = "darkblue")+
  geom_line(stat="smooth",method = "lm", formula = y ~ log(x), colour="blue", se = FALSE, alpha = 0.2, aes(x=GammaForm, y=BetaSimForm*45))+
  #geom_point(shape=19, size = 2, alpha = 0.4, colour = "red3", data = output, aes(x=GammaForm, y=AlphaForm/40))+
  scale_x_continuous(name = "Mean gamma diversity [species] per random sample")+#,  expand = c(0,1)) +
  scale_y_continuous(name = "Mean alpha (red) diversity [species] per random sample",
                     sec.axis = sec_axis(~./45, name = "Mean beta diversity [Simpson's Metric] per random sample"))
  #                   sec.axis = sec_axis(name = "mean alpha diversity per random sample [sp]"))+
  #ggtitle("Random sample plots")+
 # annotate("text", x=300, y=0.1, label=paste("repetitions of 20 random samplings: ", n))
#ggsave(("Fig_S6_global-nullmodel.png"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S6_global-nullmodel.pdf"), plot = last_plot(), width=20, height=19, units = "cm")


####Building Beta slope functions for BetaSim ####

slope_beta <- function(xx){
  lm1 <- lm(xx$BetaSimForm~xx$GammaForm)
  pval <- summary(lm(xx$BetaSimForm~xx$GammaForm))$coefficients[2,4]
  # add CI for the slope
  SE <- summary(lm1)[["coefficients"]][2,2] 
  upper = lm1$coefficients[[2]] + 1.96*SE
  lower = lm1$coefficients[[2]] - 1.96*SE
  return(data.frame(slope = lm1$coefficients[[2]], pvalue = pval, lowerCI=lower, upperCI=upper))
}
df.slope.b <- ddply(Data, .(Period.new), slope_beta)
output.slope_beta <- ddply(output, .(repetition), slope_beta) # null model beta slope
df.slope.b$signif <- df.slope.b$pvalue<0.05

####Figure S7 Beta Slopes####
ggplot(df.slope.b, aes(x=Period.new, y=slope))+
    geom_point(size=3)+
    geom_line(aes(x=1:11,y=df.slope.b$slope), colour = "black", size=0.7)+
    geom_errorbar(aes(x=Period.new, ymin=lowerCI, ymax=upperCI), width=0.1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    geom_hline(aes(yintercept = mean(output.slope_beta$slope)), colour = "red")+
    scale_x_discrete("Period")+
    scale_y_continuous("Slope")+
    geom_hline(aes(yintercept = mean(output.slope_beta$slope)+
                     2*sd(output.slope_beta$slope)), color="red", lty=2)+
    geom_hline(aes(yintercept = mean(output.slope_beta$slope)-
                     2*sd(output.slope_beta$slope)), color="red", lty=2)
#ggsave(("Fig_S7_beta-slopes.png"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S7_beta-slopes.pdf"), plot = last_plot(), width=20, height=19, units = "cm")

####FigureS8: Full data vs. Subset####
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
#ggsave(("Fig_S8_subset_vs_full.png"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S8_subset_vs_full.pdf"), plot = last_plot(), width=20, height=19, units = "cm")


####FigureS9: Full data vs. Subset Periods####
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
#ggsave(("Fig_S9_subset_vs_full-periods.png"), plot = last_plot(), width=20, height=19, units = "cm", dpi = 300)
#ggsave(("Fig_S9_subset_vs_full-periods.pdf"), plot = last_plot(), width=20, height=19, units = "cm")



##### Figure S10: Goodness of fit for different functions####
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
#ggsave(("Fig_S10_goodness-of-fit.png"), plot = last_plot(), width=14, height=13, units = "cm", dpi = 300)
#ggsave(("Fig_S10_goodness-of-fit.pdf"), plot = last_plot(), width=14, height=13, units = "cm")
