
# Models used in manuscript 'importance of kangaroo mother'

# Questions  1 (Maternal repeatability) and 2 (Covariances between behaviour and survival)

data.merged11<- read.csv(file = "data.merged11.csv", header = T)

library(MCMCglmm)

# Model I: 4-trait multivariate model####
MultiPrior8 <- list(G = list(G1 = list(V = diag(4), nu = 4)), 
                    R = list(R1=list(V=diag(1), nu=1), 
                             R2 = list(V = diag(3), nu = 3,fix=3)))
times2.2.1<- 2500                           
AverJuvTraits_Mother_MultiMCMCModel4.2<-MCMCglmm(cbind(FIDAdFem,AvFIDJuv,AvPYMov, JuvSurv) ~ trait - 1 +
                                                 at.level(x=trait,level=1):(East+ GroupSize+ TestNumber+ RelFemAlong+ Observer+ DifferenceDays)+
                                                 at.level(x=trait,level=1:2):(Year)+
                                                 at.level(x=trait,level=2):(Sex)+
                                                 at.level(x=trait,level=3):(AgeCapture+Sex+Year)+
                                                 at.level(x=trait,level=4):(Sex+Year),
                                               random = ~ us(trait):MotherID, 
                                               rcov = ~idh(at.level(x=trait, level=1)):units+us(at.level(x=trait, level=2:4)):units, family = c("gaussian","gaussian", "gaussian","categorical"), 
                                               data = data.merged11, prior = MultiPrior8, verbose = T,
                                               nitt=1300*times2.2.1,thin=1*times2.2.1,burnin=300*times2.2.1)
summary(AverJuvTraits_Mother_MultiMCMCModel4.2)
plot(AverJuvTraits_Mother_MultiMCMCModel4.2$VCV)
autocorr(AverJuvTraits_Mother_MultiMCMCModel4.2$Sol)
autocorr(AverJuvTraits_Mother_MultiMCMCModel4.2$VCV)
save(AverJuvTraits_Mother_MultiMCMCModel4.2,file="AverJuvTraits_Mother_MultiMCMCModel4.2_All2019SurvData2500x1300IttFinalModel1.rdata")
load("~/Documents/Research/Kruuk-Lab/1-Maternal-effects/Multivariate-Mother-Off-Beh-Surv-Analysis/AverJuvTraits_Mother_MultiMCMCModel4.2_All2019SurvData2500x1300IttFinalModel1.rdata")


# Maternal Repeatabilitiy estimation (between non-siblings variance)####
AverJuvTraits_Mother_MultiMCMCModel4<-AverJuvTraits_Mother_MultiMCMCModel4.2

# PY Mov####
summary(AverJuvTraits_Mother_MultiMCMCModel4)
PYMAposterior.repeatability_M2 <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "traitAvPYMov:traitAvPYMov.MotherID"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitAvPYMov:traitAvPYMov.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)2:at.level(x = trait, level = 2:4)2.units"])
hist(PYMAposterior.repeatability_M2)
posterior.mode(PYMAposterior.repeatability_M2) *100
HPDinterval(PYMAposterior.repeatability_M2) *100
save(PYMAposterior.repeatability_M2,file="PYMAposterior.repeatability_M2_RepeatabilitiesMCMC_PYM_MotherID.Rdata")

# Within-mother (or residual = between offspring/siblings variance) proportion
summary(AverJuvTraits_Mother_MultiMCMCModel4)
PYMAposterior.repeatability_O <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)2:at.level(x = trait, level = 2:4)2.units"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitAvPYMov:traitAvPYMov.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)2:at.level(x = trait, level = 2:4)2.units"])
hist(PYMAposterior.repeatability_O)
posterior.mode(PYMAposterior.repeatability_O)*100
HPDinterval(PYMAposterior.repeatability_O)*100
save(JuvFIDAposterior.repeatability_M,file="JuvFIDAposterior.repeatability_M_RepeatabilitiesMCMC_JuvFID_MotherID.Rdata")


# Juv FID####
summary(AverJuvTraits_Mother_MultiMCMCModel4)
JuvFIDAposterior.repeatability_M <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "traitAvFIDJuv:traitAvFIDJuv.MotherID"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitAvFIDJuv:traitAvFIDJuv.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)1:at.level(x = trait, level = 2:4)1.units"])
hist(JuvFIDAposterior.repeatability_M)
posterior.mode(JuvFIDAposterior.repeatability_M)*100
# mean(JuvFIDAposterior.repeatability_M)
HPDinterval(JuvFIDAposterior.repeatability_M)*100
save(JuvFIDAposterior.repeatability_M,file="JuvFIDAposterior.repeatability_M_RepeatabilitiesMCMC_JuvFID_MotherID.Rdata")

# Within-mother (or residual) proportion
summary(AverJuvTraits_Mother_MultiMCMCModel4)
JuvFIDAposterior.repeatability_O <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)1:at.level(x = trait, level = 2:4)1.units"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitAvFIDJuv:traitAvFIDJuv.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)1:at.level(x = trait, level = 2:4)1.units"])
hist(JuvFIDAposterior.repeatability_O)
posterior.mode(JuvFIDAposterior.repeatability_O)
HPDinterval(JuvFIDAposterior.repeatability_O)
save(JuvFIDAposterior.repeatability_M,file="JuvFIDAposterior.repeatability_M_RepeatabilitiesMCMC_JuvFID_MotherID.Rdata")

# Juv Surv Data Scale####
summary(AverJuvTraits_Mother_MultiMCMCModel4)
library(MCMCglmm)
library(QGglmm) #for backtranformation
RepsSurv_M2<-vector(length=1000) #did not converge
for (i in 1:1000) {
  prevalues <- predict.MCMCglmm(object = AverJuvTraits_Mother_MultiMCMCModel4, it=i, type = "terms")
  qgrSurv_M2<-QGicc(predict =prevalues,
                    var.comp  = AverJuvTraits_Mother_MultiMCMCModel4$VCV[i,"traitJuvSurv.1:traitJuvSurv.1.MotherID"],
                    var.p = AverJuvTraits_Mother_MultiMCMCModel4$VCV[i,"traitJuvSurv.1:traitJuvSurv.1.MotherID"]+
                      AverJuvTraits_Mother_MultiMCMCModel4$VCV[i,"at.level(x = trait, level = 2:4)3:at.level(x = trait, level = 2:4)3.units"], 
                    model = "binom1.logit")
  RepsSurv_M2[i]<-qgrSurv_M2$icc.obs
}
save(RepsSurv_M2,file="RepsSurv_M2_RepeatabilitiesMCMC_YoungSurvival_MotherID.Rdata")
plot(as.mcmc(RepsSurv_M2))
posterior.mode(as.mcmc(RepsSurv_M2), adjust = 1)
mean(RepsSurv_M2)
HPDinterval(as.mcmc(RepsSurv_M2))
hist(RepsSurv_M2)
mean(RepsSurv_M2)

# Within-mother (or residual) proportion
library(MCMCglmm)
library(QGglmm)
RepsSurv_J<-vector(length=1000) 
for (i in 1:1000) {
  prevalues <- predict.MCMCglmm(object = AverJuvTraits_Mother_MultiMCMCModel4, it=i, type = "terms")
  qgrSurv_J<-QGicc(predict =prevalues,
                   var.comp  = AverJuvTraits_Mother_MultiMCMCModel4$VCV[i,"at.level(x = trait, level = 2:4)3:at.level(x = trait, level = 2:4)3.units"],
                   var.p = AverJuvTraits_Mother_MultiMCMCModel4$VCV[i,"traitJuvSurv.1:traitJuvSurv.1.MotherID"]+
                     AverJuvTraits_Mother_MultiMCMCModel4$VCV[i,"at.level(x = trait, level = 2:4)3:at.level(x = trait, level = 2:4)3.units"]
                   , model = "binom1.logit")
  RepsSurv_J[i]<-qgrSurv_J$icc.obs
}
plot(as.mcmc(RepsSurv_J))
posterior.mode(as.mcmc(RepsSurv_J), adjust = 1)
mean(RepsSurv_J)
HPDinterval(as.mcmc(RepsSurv_J))
hist(RepsSurv_J)
mean(RepsSurv_J)
save(RepsSurv_J2,file="RepsSurv_J_RepeatabilitiesMCMC_YoungSurvival_YoungID.Rdata")
residua1<-(posterior.mode(as.mcmc(RepsSurv_J), adjust = 1)+posterior.mode(as.mcmc(RepsSurv_M2), adjust = 1)) - 1


# JuvSur Latent Scale####
# Maternal repeatability
summary(AverJuvTraits_Mother_MultiMCMCModel4)
JuvSurvposterior.repeatability_M <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "traitJuvSurv.1:traitJuvSurv.1.MotherID"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitJuvSurv.1:traitJuvSurv.1.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)3:at.level(x = trait, level = 2:4)3.units"])
hist(JuvSurvposterior.repeatability_M)
posterior.mode(JuvSurvposterior.repeatability_M)
mean(JuvSurvposterior.repeatability_M)
HPDinterval(JuvSurvposterior.repeatability_M)

#Within mother proportion
summary(AverJuvTraits_Mother_MultiMCMCModel4)
JuvSurvposterior.repeatability_J <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)3:at.level(x = trait, level = 2:4)3.units"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitJuvSurv.1:traitJuvSurv.1.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 2:4)3:at.level(x = trait, level = 2:4)3.units"])
hist(JuvSurvposterior.repeatability_J)
posterior.mode(JuvSurvposterior.repeatability_J)
mean(JuvSurvposterior.repeatability_J)
HPDinterval(JuvSurvposterior.repeatability_J)


# Individual Repeatability####

# Mum FID####
summary(AverJuvTraits_Mother_MultiMCMCModel4)
AdFemFIDAposterior.repeatability <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "traitFIDAdFem:traitFIDAdFem.MotherID"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitFIDAdFem:traitFIDAdFem.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 1).units"])
hist(AdFemFIDAposterior.repeatability)
posterior.mode(AdFemFIDAposterior.repeatability)
HPDinterval(AdFemFIDAposterior.repeatability)
save(AdFemFIDAposterior.repeatability,file="AdFemFIDAposterior.repeatability_RepeatabilitiesMCMC_AdFemFID_MotherID.Rdata")

#Residual proportion
AdFemFIDAposterior.repeatabilityR <- AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 1).units"]/
  (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitFIDAdFem:traitFIDAdFem.MotherID"] + 
     AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "at.level(x = trait, level = 1).units"])
hist(AdFemFIDAposterior.repeatabilityR)
posterior.mode(AdFemFIDAposterior.repeatabilityR)
HPDinterval(AdFemFIDAposterior.repeatabilityR)


##### Question 3 (importance of mother for between juvenile personality variation) and 4 (heritability)

### Model II: Bivariate model MUm & Offs FID####

# write.csv(data.merged6, file = "data.merged6.csv", row.names = F)
data.merged6<- read.csv(file = "data.merged6.csv", header = T)

library(MCMCglmm)
MultiPriorxx <- list(G = list(G1 = list(V = diag(2), n = 2),
                              G2 = list(V=diag(1), n=1)), 
                     R = list(V = diag(2), n = 2)) 
timesxx<-1000
FIDCov_MotherOffspring<-MCMCglmm(cbind(FIDAdFem,FIDJuv) ~ trait - 1 +
                                   at.level(x=trait,level=1:2):(East+ GroupSize+ TestNumber+ RelFemAlong+ Year+ DifferenceDays)+
                                   at.level(x=trait,level=1):(Observer)+
                                   at.level(x=trait,level=2):(AgeatObs+ Sex),
                                 random = ~ us(trait):MotherID+idh(at.level(x=trait, level=2)):YoungID, 
                                 rcov = ~idh(trait):units, family = c("gaussian","gaussian"), 
                                 data = data.merged6, prior = MultiPriorxx, verbose = T,
                                 nitt=1300*timesxx,thin=1*timesxx,burnin=300*timesxx)
summary(FIDCov_MotherOffspring)
plot(FIDCov_MotherOffspring$VCV)
autocorr(FIDCov_MotherOffspring$Sol)
autocorr(FIDCov_MotherOffspring$VCV)
save(FIDCov_MotherOffspring,file="FIDCov_MotherOffspring_CovarianceForQuestion31000It.rdata")
load("~/Documents/PhD Data Analyses/5 Pouch Young Response/Multivariate PY Survival Analysis - Exit&Weaning/FIDCov_MotherOffspring_CovarianceForQuestion31000It.rdata")

# MFID individual repeatability####
summary(FIDCov_MotherOffspring)
FIDCovposterior.repeatability_MFID_M <- FIDCov_MotherOffspring$VCV[, "traitFIDAdFem:traitFIDAdFem.MotherID"]/
  (FIDCov_MotherOffspring$VCV[,"traitFIDAdFem:traitFIDAdFem.MotherID"] + 
     FIDCov_MotherOffspring$VCV[, "traitFIDAdFem.units"])
hist(FIDCovposterior.repeatability_MFID_M)
posterior.mode(FIDCovposterior.repeatability_MFID_M)
mean(FIDCovposterior.repeatability_MFID_M)
HPDinterval(FIDCovposterior.repeatability_MFID_M)

# MFID Residual proportion 
FIDCovposterior.repeatability_MFID_R <- FIDCov_MotherOffspring$VCV[, "traitFIDAdFem.units"]/
  (FIDCov_MotherOffspring$VCV[,"traitFIDAdFem:traitFIDAdFem.MotherID"] + 
     FIDCov_MotherOffspring$VCV[, "traitFIDAdFem.units"])
hist(FIDCovposterior.repeatability_MFID_R)
posterior.mode(FIDCovposterior.repeatability_MFID_R)
mean(FIDCovposterior.repeatability_MFID_R)
HPDinterval(FIDCovposterior.repeatability_MFID_R)

# JFID Maternal repeatability####
FIDCovposterior.repeatability_JFID_M <- FIDCov_MotherOffspring3$VCV[, "traitFIDJuv:traitFIDJuv.MotherID"]/
  (FIDCov_MotherOffspring3$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"] + 
     FIDCov_MotherOffspring3$VCV[, "at.level(x = trait, level = 2).YoungID"]+
     FIDCov_MotherOffspring3$VCV[, "traitFIDJuv.units"])
hist(FIDCovposterior.repeatability_JFID_M)
posterior.mode(FIDCovposterior.repeatability_JFID_M)
mean(FIDCovposterior.repeatability_JFID_M)
HPDinterval(FIDCovposterior.repeatability_JFID_M)

# JFID individual repeatability
FIDCovposterior.repeatability_JFID_J <- FIDCov_MotherOffspring3$VCV[, "at.level(x = trait, level = 2).YoungID"]/
  (FIDCov_MotherOffspring3$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"] + 
     FIDCov_MotherOffspring3$VCV[, "at.level(x = trait, level = 2).YoungID"]+
     FIDCov_MotherOffspring3$VCV[, "traitFIDJuv.units"])
hist(FIDCovposterior.repeatability_JFID_J)
posterior.mode(FIDCovposterior.repeatability_JFID_J)
mean(FIDCovposterior.repeatability_JFID_J)
HPDinterval(FIDCovposterior.repeatability_JFID_J)

# JFID Residual proportion 
FIDCovposterior.repeatability_JFID_R <- FIDCov_MotherOffspring$VCV[, "traitFIDJuv.units"]/
  (FIDCov_MotherOffspring$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"] + 
     FIDCov_MotherOffspring$VCV[, "at.level(x = trait, level = 2).YoungID"]+
     FIDCov_MotherOffspring$VCV[, "traitFIDJuv.units"])
hist(FIDCovposterior.repeatability_JFID_R)
posterior.mode(FIDCovposterior.repeatability_JFID_R)
mean(FIDCovposterior.repeatability_JFID_R)
HPDinterval(FIDCovposterior.repeatability_JFID_R)


# Total repeatability####
Total_Repeat<-FIDCovposterior.repeatability_JFID_M+FIDCovposterior.repeatability_JFID_J
posterior.mode(Total_Repeat)
mean(Total_Repeat)
HPDinterval(Total_Repeat)

# Maternal proportion of total R####
MaternalPro_Total<-FIDCovposterior.repeatability_JFID_M/(FIDCovposterior.repeatability_JFID_M+FIDCovposterior.repeatability_JFID_J)
posterior.mode(MaternalPro_Total)
mean(MaternalPro_Total)
HPDinterval(MaternalPro_Total)


rJuvProportion_M<-((FIDCov_MotherOffspring3$VCV[, "traitFIDJuv:traitFIDJuv.MotherID"])/
      ((FIDCov_MotherOffspring3$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"])+(FIDCov_MotherOffspring3$VCV[,"at.level(x = trait, level = 2).YoungID"])))
mean(rJuvProportion_M) 
HPDinterval(rJuvProportion_M)

rJuvProportion_O<-((FIDCov_MotherOffspring$VCV[, "at.level(x = trait, level = 2).YoungID"])/
      ((FIDCov_MotherOffspring$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"])+(FIDCov_MotherOffspring$VCV[,"at.level(x = trait, level = 2).YoungID"])))
mean(rJuvProportion_O) 
HPDinterval(rJuvProportion_O)

rJuvProportion_Total<-((FIDCov_MotherOffspring$VCV[, "at.level(x = trait, level = 2).YoungID"])+(FIDCov_MotherOffspring$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"]))/
                     ((FIDCov_MotherOffspring$VCV[,"traitFIDJuv:traitFIDJuv.MotherID"])+(FIDCov_MotherOffspring$VCV[,"at.level(x = trait, level = 2).YoungID"])+
                        (FIDCov_MotherOffspring$VCV[,"traitFIDAdFem.units"])+ 
                        (FIDCov_MotherOffspring$VCV[,"traitFIDJuv.units"]))
mean(rJuvProportion_Total) 
HPDinterval(rJuvProportion_Total)
summary(FIDCov_MotherOffspring)


# Heritability calculation####

summary(FIDCov_MotherOffspring)
h<-((FIDCov_MotherOffspring$VCV[, "traitFIDAdFem:traitFIDJuv.MotherID"])/
      (FIDCov_MotherOffspring$VCV[,"traitFIDAdFem:traitFIDAdFem.MotherID"]))*2
mean(h) 
HPDinterval(h)

AverJuvTraits_Mother_MultiMCMCModel4

h2<-((AverJuvTraits_Mother_MultiMCMCModel4$VCV[, "traitFIDAdFem:traitAvFIDJuv.MotherID"])/
      (AverJuvTraits_Mother_MultiMCMCModel4$VCV[,"traitFIDAdFem:traitFIDAdFem.MotherID"]))*2
mean(h2) 
HPDinterval(h2)
summary(AverJuvTraits_Mother_MultiMCMCModel4)

