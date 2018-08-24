#clean routine
rm (list=ls())

setwd("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/PNAS_RESUB")
#setwd("C:/Users/fossi/Documents/pRojects/MatrixNew")

library(vegan)
library(betapart)####NEWBETA_0####

mydata <- read.csv("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/PNAS_RESUB/Vetted/Jurassic1-vetted.csv", 
                   row.names = NULL, sep=";")
#mydata <- read.csv("C:/Users/fossi/Documents/pRojects/MatrixNew/Vetted/Jurassic1-vetted.csv", 
#                   row.names = NULL, sep=";")

##standardize environments: Environmental setting of samples is erratic. Here we assign each possible category to a more robust classification scheme  
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%c("intrashelf/intraplatform reef","slope/ramp reef", "reef", "reef, buildup or bioherm","platform/shelf-margin reef", "buildup or bioherm", "perireef or subreef")] <- "reef"
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%c("deep-water indet.","deep subtidal ramp","offshore indet.","offshore shelf","deep subtidal indet.","deep subtidal shelf", "offshore", "offshore ramp")] <- "deep subtidal"
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%c( "basin", "basinal (siliceous)","basinal (siliciclastic)","basinal (carbonate)")] <- "basin"
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%c("estuary/bay", "deltaic indet.", "interdistributary bay", "estuary/bay", "delta plain", "delta front", "lagoonal/restricted shallow subtidal","paralic indet.","lagoonal", "peritidal", "prodelta","marginal marine indet.")] <- "marginal marine"
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%c("foreshore", "sand shoal", "shoreface", "transition zone/lower shoreface", "coastal indet.", "open shallow subtidal", "shallow subtidal indet.")] <- "shallow marine"
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%c("0", "","unknown", "marine indet.", "carbonate indet.")] <- "unknown"
levels(mydata$collections.environment)[levels(mydata$collections.environment)%in%"slope"] <- "slope"
table(mydata$collections.environment) 

## save phylum from mydata
class <- tail(mydata, 1) 
class <- class[,-c(1:10)]
class.df <- t(class)
class.df <- cbind(class.df, rownames(class.df))
class.df <- as.data.frame(class.df)
names(class.df) <- c("class", "species")
table(class.df$class)

##deletes rows with no collection number (Class)
mydata <- mydata[-which(is.na(mydata$collections.number)),] 

##DEFINE ALLE THE OBJECTS WE NEED
formation.name <- sort(unique(mydata$collections.formation))
datasets <- paste0("Formation", 1:length(formation.name))

Formationnames <- (formation.name)
CollpF <- c() #this will be a vector which saves the numbers of Collections per Formation
Age <- c() #Define a Vector wich collects the Age of each Formation subset
#Period <-c()
AlphaForm <- c() #assigns numbers created during the loop to a vector here called "alpha" as referring to alpha diviersities
BetaWForm <- c() # assigns numbers created during the to a vector here called "gamma" as referring to gamma diviersities
GammaForm <- c() # assigns numbers created during the to a vector here called "beta" as referring to beta diviersities
BetaSForm <- c() # a vector for additional Beta, here Sorensen or whatever
BetaJForm <- c()
BetaSimForm <- c()####NEWBETA_1####
RefForm <- c()  # a vector for adding the number of References per Formation
Maxage <- c()
Minage <- c()
Environments <- c()
Plate <- c()
abc <- paste(formation.name, "classes") #Names for files
ABC <- paste(formation.name, "abund") #Names for files
xyz <- paste(formation.name, "environments") #Names for files

Ab.Echi <- c()
Ab.Brach <- c()
Ab.Moll <- c()

Ab.Env.unkn <- c()
Ab.Env.bas <- c()
Ab.Env.shall <- c()
Ab.Env.deep <- c()
Ab.Env.marg <- c()
Ab.Env.reef <- c()
Ab.Env.slope <- c()

#Loop for creating all subsets based on Formation Name
####Creating submatrices for each formation#####
for(i in 1:length(formation.name)){
  temp <- mydata[mydata$collections.formation==formation.name[i],] #create a tempory matrice for each formation
  temp <- temp[,which(!apply(temp,2,FUN = function(x){all(x == 0)}))] #delete all species with zero occurences
  assign(datasets[i], temp)
  
  Age <- c(Age, mean(temp$collections.ma_mid)) #extract the mean age of each of those temps and saves them in the Vector "Age"
  Maxage <- c(Maxage, max(temp$collections.ma_mid))
  Minage <- c(Minage, min(temp$collections.ma_mid))
  RefForm <- c(RefForm, length(unique(temp$collections.reference_no)))
  Formation <- data.frame(formation.name)
  Plate <- c(Plate, unique(temp$collections.plate)[1])
  Environments <- c(Environments, length(unique(temp$collections.environment)))
  
  #creating dataframes to assign classes to each species which we need later...   
  temp2 <- data.frame(names(temp))
  formation.class <- merge(temp2, class.df, by.x="names.temp.", by.y = "species", all.x=TRUE)
  assign(abc[i], formation.class) 
  
  #creating dataframes for abundances of each class in each formation  
  abund.class <- as.matrix(table(formation.class$class))
  abund.class <- t(abund.class)
  assign(ABC[i], abund.class)
  
  # plot(temp$collections.environment)
  #creating dataframes for abundances of environments in each formation
  abund.envo <- as.matrix(table(temp$collections.environment))
  abund.envo <- t(abund.envo)
  sum.envo <- sum(abund.envo[,1:7])
  abund.envo.perc <- apply(abund.envo,1, function(x) {x/sum.envo*100})
  abund.envo.perc <- t(abund.envo.perc)
  assign(xyz[i], abund.envo.perc)
  
  temp <- temp[,-c(1:10)] #delete columns 1 to 10
  CollpF <- c(CollpF, nrow(temp))
  
  Alpha <- c() #assigns numbers created during the loop to a vector here called "alpha" as referring to alpha diviersities
  BetaW <- c() # assigns numbers created during the to a vector here called "gamma" as referring to gamma diviersities
  Gamma <- c() # assigns numbers created during the to a vector here called "beta" as referring to beta diviersities
  BetaS <- c() # assigns numbers..-
  BetaJ <- c()
  BetaSim <- c()   ####NEWBETA_2####
  
  for(j in 1:500)
  { ###does everything inside {} n (here 500) times###
    tempt <- t(temp)  ###transpose (t,()) yields a matrix but not a dataframe 
    tempt <- as.data.frame(tempt) ###defines matrix as dataframe which is necessary to work with it###
    Subtempt <- sample(tempt, size = 20, replace = FALSE) ###samples 30 columns of the dataframe at random###
    Subtemp <- t(Subtempt) #transpose back because, well because...
    
    NewSubtemp <- Subtemp[,which(!apply(Subtemp,2,FUN = function(x){all(x == 0)}))] #removes all species(columns) with n = 0)
    species <- specnumber(NewSubtemp) #extracts species numbers of each sample of the subset, specnumber is a funciton in "Vegan"
    
    Alpha <- c(Alpha, mean(species)) #calculates the mean (=average alpha diversity of the subset) and puts it into our created vector "alpha"
    Gamma <- c(Gamma, ncol(NewSubtemp)) #calculates the overall diversity of the subset and puts it into our created vector
    BetaW <- c(BetaW, ncol(NewSubtemp)/mean((species)-1)) #calculates beta diversity (whittaker)
    
    ###make a matrix to calculate betadiversities using betadiver of "Vegan"
    Betamat <- NewSubtemp  ####DELETE THIS[,-1] #delete first column (samples) to make "betadiv" work on the matrix#####
   
    Sorensen <- betadiver(Betamat, method="sor")
    BetaS <- c(BetaS,1-mean(Sorensen))
    
    Jaccard <- betadiver(Betamat, method="j")
    BetaJ <- c(BetaJ, 1-mean(Jaccard))
   
     ####NEWBETA_3####
    NewSubtemp[NewSubtemp>0] <- 1
    class(NewSubtemp) <- "numeric"
    Multi <- beta.multi(NewSubtemp, index.family="sor")
    SimMulti <- Multi[[1]]
    BetaSim <- c(BetaSim, mean(SimMulti))
  }
  
  Ab.Echi <- c(Ab.Echi, (abund.class[1,2]))
  Ab.Brach <- c(Ab.Brach, (abund.class[1,1]))
  Ab.Moll <- c(Ab.Moll, (abund.class[1,3]))
  
  Ab.Env.unkn <- c(Ab.Env.unkn, (abund.envo.perc[1,1]))
  Ab.Env.bas <- c(Ab.Env.bas, (abund.envo.perc[1,2]))
  Ab.Env.shall <- c(Ab.Env.shall, (abund.envo.perc[1,3]))
  Ab.Env.deep <- c(Ab.Env.deep, (abund.envo.perc[1,4]))
  Ab.Env.marg <- c(Ab.Env.marg, (abund.envo.perc[1,5]))
  Ab.Env.reef <- c(Ab.Env.reef, (abund.envo.perc[1,6]))
  Ab.Env.slope <- c(Ab.Env.slope, (abund.envo.perc[1,7]))
  
  
  ####now we collect all diversities from each run for further processing
  if(i==1){
    Alpha.Jura1 <- data.frame(Alpha)
  }else{Alpha.Jura1 <- cbind(Alpha.Jura1, Alpha)}
  
  if(i==1){
    BetaW.Jura1 <- data.frame(BetaW)
  }else{BetaW.Jura1 <- cbind(BetaW.Jura1, BetaW)}
  
  if(i==1){
    BetaJ.Jura1 <- data.frame(BetaJ)
  }else{BetaJ.Jura1 <- cbind(BetaJ.Jura1, BetaJ)}
  
  if(i==1){
    BetaS.Jura1 <- data.frame(BetaS)
  }else{BetaS.Jura1 <- cbind(BetaS.Jura1, BetaS)}
  
  if(i==1){
    Gamma.Jura1 <- data.frame(Gamma)
  }else{Gamma.Jura1 <- cbind(Gamma.Jura1, Gamma)}
  
  ####NEWBETA_4####
  if(i==1){
    BetaSim.Jura1 <- data.frame(BetaSim)
  }else{BetaSim.Jura1 <- cbind(BetaSim.Jura1, BetaSim)}
  
  ####Vectors just for mean diversities per formation  
  
  AlphaForm <- c(AlphaForm, mean(Alpha)) 
  GammaForm <- c(GammaForm, mean(Gamma))
  BetaWForm <- c(BetaWForm, mean(BetaW))
  BetaSForm <- c(BetaSForm, mean(BetaS)) 
  BetaJForm <- c(BetaJForm, mean(BetaJ))
  BetaSimForm <- c(BetaSimForm, mean(BetaSim))####NEWBETA_5####
  
  print(i) #to see the progress
  
} #End of looping

Ab.Echi <- as.data.frame(Ab.Echi)
Ab.Brach <- as.data.frame(Ab.Brach)
Ab.Moll <- as.data.frame(Ab.Moll)

colnames(Ab.Echi) <- "Echinodermata"
colnames(Ab.Brach) <- "Brachiopoda"
colnames(Ab.Moll) <- "Mollusca"

Ab.Env.unkn <- as.data.frame(Ab.Env.unkn)
Ab.Env.bas <- as.data.frame(Ab.Env.bas)
Ab.Env.shall <- as.data.frame(Ab.Env.shall)
Ab.Env.deep <- as.data.frame(Ab.Env.deep)
Ab.Env.marg <- as.data.frame(Ab.Env.marg)
Ab.Env.reef <- as.data.frame(Ab.Env.reef)
Ab.Env.slope<- as.data.frame(Ab.Env.slope)

colnames(Ab.Env.unkn) <- "unkown"
colnames(Ab.Env.bas) <- "basin"
colnames(Ab.Env.shall) <- "shallow subtidal"
colnames(Ab.Env.deep) <- "deep subtidal"
colnames(Ab.Env.marg) <- "marginal"
colnames(Ab.Env.reef) <- "reefal"
colnames(Ab.Env.slope) <- "slope"

Duration <- (Maxage - Minage)

Jura1Dat <- cbind(Age, AlphaForm, BetaSimForm,####NEWBETA_6####
                  BetaWForm, BetaSForm, BetaJForm, GammaForm, CollpF, 
                  Duration, RefForm, Formation, Environments, Ab.Env.unkn, Ab.Env.marg,
                  Ab.Env.shall, Ab.Env.reef, Ab.Env.deep, Ab.Env.slope, Ab.Env.bas,
                  Plate, Ab.Echi, Ab.Brach, Ab.Moll)

write.csv(Jura1Dat, file = "10_Jura1Dat.csv")