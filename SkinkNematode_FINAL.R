# Wog Wog skink parasite study
# Julian Resasco (jresasco@colorado.edu)
#####################
#Codes
#Site:      1 - 188
#Replicate: 1 - 6
#Size: 	   1=small, 2=medium, 3=large, (4 = cont)
#Edge:	   1=inner, 2=outer
#Topography:1=slope, 2=drainage
#Treatment: 1=fragmented, 2=controls, 3=pine matrix
#Notes about Treatment	No pine sites for 1985-1986. In 1985-86, treatment 1 was not yet fragmented.
#######################
library(pbkrtest)
library(brglm)
library(brglm2)
library(lme4)
library(car)
library(MASS)


# set working directory
setwd("/Users/julianresasco/Documents/Data/Wog Wog/Skink parasite paper")

# read in data
Data <- read.csv("HedrurisGutData.csv")

# summary statistics
tapply(Data$HWCount,Data$Treat, sum)
tapply(Data$HWPres,Data$Treat, sum)
table(Data$Treat,Data$HWPres)

tapply(Data$HWPres,Data$Treat, mean)
tapply(Data$HWPres,Data$Rep, sum)


# Color palette for figures
WogWogCols <- palette(c("#556B2F","darkolivegreen3","#CD7F32"))


#Fig2A###
#GLMMs for Hedruris wogwogensis gut content presence with Matrix (set to 3.5)

DataTreatSub <- subset(Data, TreatNum < 3.5) # for models including the matrix make this "3.5"; for analyses without the matrix use "2.5"
head(DataTreatSub)

###Make factors
DataTreatSub$TreatNum <- factor(DataTreatSub$TreatNum, levels=c(2,1,3))
DataTreatSub$Rep <- factor(DataTreatSub$Rep)

tapply(DataTreatSub$HWCount,DataTreatSub$Treat, mean)

# Null model 
m0 <- brglm(HWPres ~  1, data=DataTreatSub, family=binomial) 
m0.1 <- glm(HWPres ~  1, data=DataTreatSub, family=binomial) 
m0BC <- update(m0.1, method = "brglmFit")

# Treat model 
m1  <- brglm(HWPres ~ TreatNum , data=DataTreatSub, family=binomial) 
m1.1  <- glm(HWPres ~ TreatNum , data=DataTreatSub, family=binomial)
m1BC <- update(m1.1, method = "brglmFit")

anova(m0,m1) # brglm method
anova(m0BC,m1BC) # brglm2 method
PBmodcomp(m1BC, m0BC, nsim = 100, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
SUM <-summary(m1BC)
CI<-confint(m1BC) # CIs used

par(mfrow=c(2,2))
par(mar=c(5, 5, 2, 2) + 0.1)
plot(NA,NA,xlim=c(0,3),ylim=c(-8,8), xlab = "", xaxt='n', ylab ="Effect size of nematode \n prevalence in skinks",las=1,axes = F)
axis(1,at = 1:2, labels = c("Fragments","Matrix"),las =2)
axis(2,at =seq(from =-8, to =8, by = 2),las =2)
points(1:2,SUM$coefficients[2:3,1],pch=21,bg="white", cex=1.4)
abline(h=0,col="grey", lty=2)
abline(h=log(0.5),col="grey", lty=1)
abline(h=log(2),col="grey", lty=1)
box()

CIlo<-CI[2:3,1]
CIhi<-CI[2:3,2]
segments(1:2,CIlo,1:2,CIhi)

points(1:2,SUM$coefficients[2:3,1],pch=21,bg=WogWogCols[2:3], cex=1.5)
###

#Fig2B###
# #GLMMs for Hedruris wogwogensis gut content intensity with Matrix (set to 3.5)
DataTreatSub <- subset(Data, TreatNum < 3.5) # for models including the matrix make this "3.5"; for analyses without the matrix use "2.5"
head(DataTreatSub)

###Make factors
DataTreatSub$TreatNum <- factor(DataTreatSub$TreatNum, levels=c(2,1,3))
DataTreatSub$Rep <- factor(DataTreatSub$Rep)

tapply(DataTreatSub$HWCount,DataTreatSub$Treat, mean)
tapply(DataTreatSub$HWPres,DataTreatSub$Treat, mean)
#hist(DataTreatSub$HWCount)

# Null model 
m0 <- glm(HWCount ~  1, data=DataTreatSub, family=poisson) 
m0BC <- update(m0, method = "brglmFit")
m0nb <- glm.nb(HWCount ~  1, data=DataTreatSub) 

# Treat model 
m1 <- glm(HWCount ~  TreatNum, data=DataTreatSub, family=poisson) 
m1BC <- update(m1, method = "brglmFit")
m1nb <- glm.nb(HWCount ~  TreatNum, data=DataTreatSub) 

# test Treat vs Null
anova(m0,m1) #no bias correction Poisson
anova(m0BC,m1BC) # bias correction Poisson
anova(m0nb,m1nb) # negative binomial

# comparison of Poisson vs. NB
AIC(m1,m1nb)

#model summaries
summary(m1) # no bias correction Poisson
summary(m1BC) # bias correction Poisson
summary(m1nb) # negative binomial
glm.convert(sum.nb) # negative binomial convert

#confidence intervals ERRORS
confint(m1) # no bias correction Poisson 
confint(m1BC) # bias correction Poisson
confint(m1nb) # negative binomial

SUM <-summary(m1nb)
CI<-confint.nb(m1nb) # CIs used

plot(NA,NA,xlim=c(0,3),ylim=c(-8,8), xlab = "", xaxt='n', ylab ="Effect size of nematode \n intensity in skinks",las=1,axes = F)
axis(1,at = 1:2, labels = c("Fragments","Matrix"),las =2)
axis(2,at =seq(from =-8, to =8, by = 2),las =2)
points(1:2,SUM$coefficients[2:3,1],pch=21,bg="white", cex=1.4)
abline(h=0,col="grey", lty=2)


CIlo<-CI[2:3,1]
CIhi<-CI[2:3,2]
segments(1:2,CIlo,1:2,CIhi)
points(1:2,SUM$coefficients[2:3,1],pch=21,bg=WogWogCols[2:3], cex=1.5)
box()









##############################
# Amphipod gut models

Data2 <-read.csv("AmphiGutData.csv")
head(Data2)
# summary statistics
tapply(Data2$AmphiPres,Data2$Treat, sum)
table(Data2$Treat,Data2$AmphiPres)


#frag models
Data2TreatSub <- subset(Data2, TreatNum < 2.5) # for models including the matrix make this "3.5"; for analyses without the matrix use "2.5"
head(Data2TreatSub)
###Make factors
Data2TreatSub$TreatNum <- factor(Data2TreatSub$TreatNum, levels=c(2,1,3))
Data2TreatSub$Rep <- factor(Data2TreatSub$Rep)

tapply(Data2TreatSub$AmphiPres,Data2TreatSub$Treat, sum)


# Amphipod gut models

# Null model 
m0 <- glm(AmphiPres ~  1 , data=Data2TreatSub, family=binomial) 

# Treat model 
m1  <- glm(AmphiPres ~ TreatNum , data=Data2TreatSub, family=binomial) 

anova(m0,m1)
PBmodcomp(m1, m0, nsim = 100, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
summary(m1)
#confint(m1) # CIs used

# Size model
m2 <- glm(AmphiPres ~ TreatNum + SizeNum , data=Data2TreatSub, family=binomial) 

anova(m1,m2) 
PBmodcomp(m2, m1, nsim = 100, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
summary(m2)

# topo model
m3 <- glm(AmphiPres ~ TreatNum + TopoNum, data=Data2TreatSub, family=binomial) 

anova(m1,m3)
PBmodcomp(m3, m1, nsim = 100, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
summary(m3)

# edge model
m4 <- glm(AmphiPres ~ TreatNum + ProxNum, data=Data2TreatSub, family=binomial) 

anova(m1,m4) 
PBmodcomp(m4, m1, nsim = 100, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
summary(m4)
#confint(m4)

Anova(m4)


Data2TreatSub <- subset(Data2, TreatNum < 3.5) # for models including the matrix make this "3.5"; for analyses without the matrix use "2.5"
head(Data2TreatSub)


###Make factors
Data2TreatSub$TreatNum <- factor(Data2TreatSub$TreatNum, levels=c(2,1,3))
Data2TreatSub$Rep <- factor(Data2TreatSub$Rep)

tapply(Data2TreatSub$AmphiPres,Data2TreatSub$Treat, sum)


#GLMs for Hedruris wogwogensis gut content presence with Matrix (set to 3.5)
# Null model 
m0 <- brglm(AmphiPres ~  1, data=Data2TreatSub, family=binomial) 

# Treat model 
m1  <- brglm(AmphiPres ~ TreatNum , data=Data2TreatSub, family=binomial) 
anova(m0,m1)
PBmodcomp(m1, m0, nsim = 100, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
SUM2<-summary(m1)
CI2<-confint(m1) # CIs used


##
#Fig2C###
plot(NA,NA,xlim=c(0,3),ylim=c(-6,6), xlab = "", xaxt='n', ylab ="Effect size of amphipod \n prevalence in skinks",las=2)
axis(1,at = 1:2, labels = c("Fragments","Matrix"),las =2)
points(1:2,SUM2$coefficients[2:3,1],pch=21,bg=col, cex=1.4)
abline(h=0,col="grey", lty=2)
abline(h=log(0.5),col="grey", lty=1)
abline(h=log(2),col="grey", lty=1)

CIlo<-c(CI2[2,1],-3.245193)
CIhi<-CI2[2:3,2]
segments(1:2,CIlo,1:2,CIhi)
arrows(2, -3.245193, 2,-6.4, length = 0.1, angle = 22,
       code = 2, col = par("fg"), lty = par("lty"),
       lwd = par("lwd"))
points(1:2,SUM2$coefficients[2:3,1],pch=21,bg=WogWogCols[2:3], cex=1.5)
###




##############################
# Amphipod Pitfall Models

Data3 <-read.csv("AmphiPitData.csv")
head(Data3)

#summary statistics
sum(Data3$AmpCt)
tapply(Data3$AmpCt,Data3$Fragmentation,mean)

head(Data3)
#frag models
Data3TreatSub <- subset(Data3, Fragmentation < 2.5) # for models including the matrix make this "3.5"; for analyses without the matrix use "2.5"
head(Data3TreatSub)
###Make factors
Data3TreatSub$Fragmentation <- factor(Data3TreatSub$Fragmentation, levels=c(2,1,3))
Data3TreatSub$Blocks <- factor(Data3TreatSub$Blocks)

# Amphipod pitfall models

# Null model 
m0 <- glmer(AmpCt ~  (1|Blocks) , data=Data3TreatSub, family=poisson) 

# Treat model 
m1  <- glmer(AmpCt ~ Fragmentation + (1|Blocks) , data=Data3TreatSub, family=poisson) 

anova(m0,m1)    
summary(m1)
#confint(m1) # CIs used

# Size model
m2 <- glmer(AmpCt ~ Fragmentation + Size14 + (1|Blocks) , data=Data3TreatSub, family=poisson) 

anova(m1,m2)  
summary(m2)

# topo model
m3 <- glmer(AmpCt ~ Fragmentation + Topo + (1|Blocks) , data=Data3TreatSub, family=poisson) 

anova(m1,m3)  
summary(m3)

# edge model
m4 <- glmer(AmpCt ~ Fragmentation + Edge +(1|Blocks) , data=Data3TreatSub, family=poisson) 

anova(m1,m4)  
summary(m4)
#confint(m4)

Anova(m4)


#Amphipod pitfall Matrix models

Data3TreatSub <- subset(Data3, Fragmentation < 3.5) # for models including the matrix make this "3.5"; for analyses without the matrix use "2.5"
head(Data3TreatSub)
tapply(Data3TreatSub$AmpCt,Data3TreatSub$Fragmentation, mean)

###Make factors
Data3TreatSub$Fragmentation <- factor(Data3TreatSub$Fragmentation, levels=c(2,1,3))
Data3TreatSub$Blocks <- factor(Data3TreatSub$Blocks)


# Amphipod models
#GLMMs for Amphipod gut content presence with Matrix (set to 3.5)
# Null model 
m0 <- glmer(AmpCt ~  (1|Blocks), data=Data3TreatSub, family=poisson) 

# Treat model 
m1  <- glmer(AmpCt ~ Fragmentation +(1|Blocks) , data=Data3TreatSub, family=poisson) 
anova(m0,m1)    
summary(m1)
#confint(m1) # CIs used

#qqnorm(residuals(m1))
#plot(m1)

SUM3<-summary(m1)
CI3<-confint(m1) # CIs used
tapply(Data3TreatSub$AmpCt,Data3TreatSub$Fragmentation, mean)


##
#Fig2D###

plot(NA,NA,xlim=c(0,3),ylim=c(-3,3), xlab = "", xaxt='n', ylab ="Effect size of amphipod \n counts in pitfall traps",las=2)
axis(1,at = 1:2, labels = c("Fragments","Matrix"),las =2)
abline(h=0,col="gray", lty=2)

CIlo<-CI3[3:4,1]
CIhi<-CI3[3:4,2]
segments(1:2,CIlo,1:2,CIhi)
points(1:2,SUM3$coefficients[2:3,1],pch=21,bg=WogWogCols[2:3], cex=1.5)

#Fig2C###exp effects

plot(NA,NA,xlim=c(0,3),ylim=c(0,2), xlab = "", xaxt='n', ylab ="Effect size of amphipod \n counts in pitfall traps")
abline(h=1,col="gray", lty=2)

CIlo<- exp(CI3[3:4,1])
CIhi<- exp(CI3[3:4,2])
segments(1:2,CIlo,1:2,CIhi)
points(1:2,-exp(SUM3$coefficients[2:3,1]),pch=21,bg=WogWogCols[2:3], cex=1.5)


########
# Logistic regression

Logistic <- read.csv("LogisticRegAmp.csv")
head(Logistic)
LogiFrag <- subset(Logistic,TreatNum == 1)
LogiCont <- subset(Logistic,TreatNum == 2)
LogiMat  <- subset(Logistic,TreatNum == 3)

# AmpYr
#par(mfrow=c(2,2))
plot(log(Logistic$AmpYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiFrag$AmpYr+1), jitter(LogiFrag$HWPres,amount = 0.02), bg=WogWogCols[2],pch = 21)
points(log(LogiCont$AmpYr+1), jitter(LogiCont$HWPres,amount = 0.02), bg=WogWogCols[1],pch = 21)
points(log(LogiMat$AmpYr+1),  jitter(LogiMat$HWPres,amount = 0.02), bg=WogWogCols[3],pch = 21)

Amp<-seq(0,5000,.1)
newAmps <-data.frame(AmpYr=Amp)
predicted.prob <- predict(model,newAmps,type="resp")
lines(predicted.prob~log(Amp+1),type="l",col="black")

plot(log(Logistic$AmpYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiFrag$AmpYr+1), jitter(LogiFrag$HWPres,amount = 0.02), bg=WogWogCols[2],pch = 21)

plot(log(Logistic$AmpYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiCont$AmpYr+1), jitter(LogiCont$HWPres,amount = 0.02), bg=WogWogCols[1],pch = 21)

plot(log(Logistic$AmpYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiMat$AmpYr+1),  jitter(LogiMat$HWPres,amount = 0.02), bg=WogWogCols[3],pch = 21)


model <- glm (HWPres ~ AmpYr, data = Logistic, family = binomial)
summary(model)

model.1 <- glm (HWPres ~ AmpYr, data = LogiFrag, family = binomial)
summary(model.1)
model.2 <- glm (HWPres ~ AmpYr, data = LogiCont, family = binomial)
summary(model.2)
model.3 <- glm (HWPres ~ AmpYr, data = LogiMat, family = binomial)
summary(model.3)



######


# AmpSameYr
par(mfrow=c(2,2))
plot(log(Logistic$AmpAllYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiFrag$AmpAllYr+1), jitter(LogiFrag$HWPres,amount = 0.02), bg=WogWogCols[2],pch = 21)
points(log(LogiCont$AmpAllYr+1), jitter(LogiCont$HWPres,amount = 0.02), bg=WogWogCols[1],pch = 21)
points(log(LogiMat$AmpAllYr+1),  jitter(LogiMat$HWPres,amount = 0.02), bg=WogWogCols[3],pch = 21)

Amp<-seq(0,5000,.1)
newAmps <-data.frame(AmpAllYr=Amp)
predicted.prob <- predict(model,newAmps,type="resp")
lines(predicted.prob~log(Amp+1),type="l",col="black")

plot(log(Logistic$AmpAllYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiFrag$AmpAllYr+1), jitter(LogiFrag$HWPres,amount = 0.02), bg=WogWogCols[2],pch = 21)

plot(log(Logistic$AmpAllYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiCont$AmpAllYr+1), jitter(LogiCont$HWPres,amount = 0.02), bg=WogWogCols[1],pch = 21)

plot(log(Logistic$AmpAllYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiMat$AmpAllYr+1),  jitter(LogiMat$HWPres,amount = 0.02), bg=WogWogCols[3],pch = 21)


model <- glm (HWPres ~ AmpAllYr, data = Logistic, family = binomial)
summary(model)

model.1 <- glm (HWPres ~ AmpAllYr, data = LogiFrag, family = binomial)
summary(model.1)
model.2 <- glm (HWPres ~ AmpAllYr, data = LogiCont, family = binomial)
summary(model.2)
model.3 <- glm (HWPres ~ AmpAllYr, data = LogiMat, family = binomial)
summary(model.3)


########
# Logistic regression year lags

#AmpPit <- read.csv("All amphipod data (blocked by year) - Copy.csv")

#SameYr   <- vector(mode="numeric", length=148)
#OneYrlag <- vector(mode="numeric", length=148)
#TwoYrlag <- vector(mode="numeric", length=148)
#ThrYrlag <- vector(mode="numeric", length=148)

#i = 1
#for (i in i:148) {
#    SameYr[i] <- AmpPit[(AmpPit$year == Logistic[i,1]) & (AmpPit$pit == Logistic[i,2]),3]
#    OneYrlag[i] <- AmpPit[(AmpPit$year == Logistic[i,1]-1) & (AmpPit$pit == Logistic[i,2]),3]
#    TwoYrlag[i] <- AmpPit[(AmpPit$year == Logistic[i,1]-2) & (AmpPit$pit == Logistic[i,2]),3]
#    ThrYrlag[i] <- AmpPit[(AmpPit$year == Logistic[i,1]-3) & (AmpPit$pit == Logistic[i,2]),3]
#    i=i+1
#}

#Logistic2 <- cbind(Logistic,SameYr,OneYrlag,TwoYrlag,ThrYrlag)
#write.csv(Logistic2, "LogisticRegAmp2.csv")

# lags analysis

head(Logistic)

# lag models
OneYrlag.mod <- glm (HWPres ~ OneYrlag, data = Logistic, family = binomial)
TwoYrlag.mod <- glm (HWPres ~ TwoYrlag, data = Logistic, family = binomial)
ThrYrlag.mod <- glm (HWPres ~ ThrYrlag, data = Logistic, family = binomial)

# combined years models
Sameyr.mod <- glm (HWPres ~ SameYr, data = Logistic, family = binomial)
Last2yrs.mod <- glm (HWPres ~ Last2yrs, data = Logistic, family = binomial)
Last3yrs.mod <- glm (HWPres ~ Last3yrs, data = Logistic, family = binomial)
Last4yrs.mod <- glm (HWPres ~ Last4yrs, data = Logistic, family = binomial)

summary(Last4yrs.mod)

# all years model
AmpAllYr.mod <- glm (HWPres ~ AmpAllYr, data = Logistic, family = binomial)

summary(Last4yrs.mod)

require(AICc)
require(MuMIn)
AICc(Last2yrs.mod,Last3yrs.mod,Last4yrs.mod,Sameyr.mod)



dev.off
par(mar=c(5, 5, 2, 2) + 0.1)
# Last4yrs
#Fig 2D
model <- glm (HWPres ~ Last4yrs, data = Logistic, family = binomial)
summary(model)
plot(log(Logistic$Last4yrs+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Probability of nematode \n occurence in skink", axes = F)
axis(1,las=1)
axis(2,las=2)
box()
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiFrag$Last4yrs+1), jitter(LogiFrag$HWPres,amount = 0.02), bg=WogWogCols[2],pch = 21)
points(log(LogiCont$Last4yrs+1), jitter(LogiCont$HWPres,amount = 0.02), bg=WogWogCols[1],pch = 21)
points(log(LogiMat$Last4yrs+1),  jitter(LogiMat$HWPres,amount = 0.02), bg=WogWogCols[3],pch = 21)

Amp<-seq(0,5000,.1)
newAmps <-data.frame(Last4yrs=Amp)
predicted.prob <- predict(model,newAmps,type="resp")
lines(predicted.prob~log(Amp+1),type="l",col="black")

plot(log(Logistic$Last4yrs+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiFrag$Last4yrs+1), jitter(LogiFrag$HWPres,amount = 0.02), bg=WogWogCols[2],pch = 21)

plot(log(Logistic$Last4yrs+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiCont$Last4yrs+1), jitter(LogiCont$HWPres,amount = 0.02), bg=WogWogCols[1],pch = 21)

plot(log(Logistic$AmpAllYr+1), Logistic$HWPres,col="white",ylim=c(-0.1,1.1),
     xlab= "log(amphipod abundance + 1)", ylab = "Nematode prevalence")
abline(h=1,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
points(log(LogiMat$Last4yrs+1),  jitter(LogiMat$HWPres,amount = 0.02), bg=WogWogCols[3],pch = 21)




model.1 <- glm (HWPres ~ Last4yrs, data = LogiFrag, family = binomial)
summary(model.1)
model.2 <- glm (HWPres ~ Last4yrs, data = LogiCont, family = binomial)
summary(model.2)
model.3 <- glm (HWPres ~ Last4yrs, data = LogiMat, family = binomial)
summary(model.3)










#Amp yr

plot(Logistic$BioYear,Logistic$HWPres)
Tab <-table(Logistic$BioYear,Logistic$HWPres)
Years<-unique(Logistic$BioYear)

Tab <-cbind(Tab,sort(Years))
dev.off()
plot(Tab[,3],(Tab[,2]/(Tab[,1]+Tab[,2]))/0.55555556,ylim=c(0,1),type="b", ylab= "Relative infection", xlab="Year")



Tab2<-cbind(tapply(Logistic$AmpYr/83.640000,Logistic$BioYear,mean),sort(Years))
points(Tab2[,2],Tab2[,1],col="red",type="b")

max(Logistic$AmpYr)

plot(LogiFrag$BioYear,LogiFrag$HWPres)
Tab <-table(LogiFrag$BioYear,LogiFrag$HWPres)
Years<-unique(LogiFrag$BioYear)

Tab <-cbind(Tab,sort(Years))
dev.off()
plot(Tab[,3],(Tab[,2]/(Tab[,1]+Tab[,2]))/0.25000000,ylim=c(0,1),type="b", ylab= "Relative infection", xlab="Year")


tapply(LogiFrag$AmpYr,LogiFrag$BioYear,mean)
Tab2<-cbind(tapply(LogiFrag$AmpYr/86.100000,LogiFrag$BioYear,mean),sort(Years))
points(Tab2[,2],Tab2[,1],col="red",type="b")

max(LogiFrag$AmpYr)


