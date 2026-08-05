library(oSCR)
library(dplyr)
library(AICcmodavg)
library(ggplot2)
library(corrplot)
library(simplecolors)

######NE STUDY AREA MODELS########

#import capture file
NEedf = read.csv("NE_S_capth.csv")
head(NEedf)

#import trap file, 1 for each session (1 session per summer)
NEtdf1 = read.csv("NE_18_S_trap.csv")
NEtdf2 = read.csv("NE_19_S_trap.csv")
NEtdf3 = read.csv("NE_20_S_trap.csv")

#rescale coordinates
NEtdf1$X <- NEtdf1$X/1000; NEtdf1$Y <- NEtdf1$Y/1000
NEtdf2$X <- NEtdf2$X/1000; NEtdf2$Y <- NEtdf2$Y/1000
NEtdf3$X <- NEtdf3$X/1000; NEtdf3$Y <- NEtdf3$Y/1000

#import mask file with 3 covariates: human footprint index (HFI, 0-10,000), 
#forest cover (FC, 0-1), and human population density (pop, #/km2)
NEmask = read.csv("NE_mask.csv")
#rescale coordinates
NEmask$X <- NEmask$X/1000; NEmask$Y <- NEmask$Y/1000

#inspect covariates
names(NEmask)
hist(NEmask$HFI) ##ranges from 0-10000, right-skewed
hist(log(NEmask$HFI + 1)) #left-skewed, looks worse
hist(NEmask$pop) #0-400, extremely 0-inflated
hist(log(NEmask$pop+1)) #skewed towards 0, logging helps but still very 0-inflated
hist(NEmask$FC) #0-1, skewed toward 1
hist(asin(sqrt(NEmask$FC))) #helps a little

#add scaled/transformed covariates to mask file
NEmask$n.HFI <-scale(NEmask$HFI)
NEmask$n.FC <-scale(asin(sqrt(NEmask$FC)))
NEmask$n.pop <-scale(log(NEmask$pop+1))

#check for correlations among covariates
NEcovar <- select(NEmask, n.HFI, n.FC, n.pop)
NECorCovar<-cor(NEcovar)  #correlation matrix of covariates
NECorCovar
#HFI vs pop= 0.70, HFI vs FC = -0.52, pop vs FC = -0.49

##create mask dataframe for each session (same for all 3)
NEssDF <- list(data.frame(NEmask), data.frame(NEmask), data.frame(NEmask))

# create the scrFrame using data2oscr
NEsf <- data2oscr(edf = NEedf, #the EDF
                sess.col = 1, #session column
                id.col = 2, #individual column
                occ.col = 3, #occasion column
                trap.col = 4, #trap column
                tdf = list(NEtdf1,NEtdf2, NEtdf3), #list of TDFs,
                K = c(1,1,1), #makes counts of scats binary in each cell
                trapcov.names = c("Effort"), #covariate names
                tdf.sep = "/",
                ntraps = c(nrow(NEtdf1),nrow(NEtdf2),nrow(NEtdf3)) #no. traps vector
)$scrFrame

class(NEssDF) <- "ssDF"

str(NEssDF)

plot(NEsf,jit = 2)
plot(ssDF = NEssDF, scrFrame = NEsf)

##Models with binary "clog" encounters, and log Effort as covariate on p0
##HFI and pop too correlated to include in same model
##models take a while to run

NEm0 <- oSCR.fit(scrFrame=NEsf,ssDF=NEssDF,DorN="D",encmod="CLOG",
                 trimS=NEsf$mdm,
                 model=list(D ~ 1,p0 ~ log(Effort),sig ~ 1))
NEm_HFI <- oSCR.fit(scrFrame=NEsf,ssDF=NEssDF,DorN="D",encmod="CLOG",
                  trimS=NEsf$mdm,
                  model=list(D~n.HFI,p0~log(Effort),sig~1))

NEm_FC <- oSCR.fit(scrFrame=NEsf,ssDF=NEssDF,DorN="D",encmod="CLOG",
                 trimS=NEsf$mdm,
                 model=list(D~n.FC,p0~log(Effort),sig~1))

NEm_pop <- oSCR.fit(scrFrame=NEsf,ssDF=NEssDF,DorN="D",encmod="CLOG",
                model=list(D~n.pop,p0~log(Effort),sig~1))

NEm_HFI_FC <- oSCR.fit(scrFrame=NEsf,ssDF=NEssDF,DorN="D",encmod="CLOG",
                     trimS=NEsf$mdm,
                     model=list(D~n.FC+n.HFI,p0~log(Effort),sig~1))

NEm_pop_FC <- oSCR.fit(scrFrame=NEsf,ssDF=NEssDF,DorN="D",encmod="CLOG",
                     trimS=NEsf$mdm,
                     model=list(D~n.FC+n.pop,p0~log(Effort),sig~1))

#Comparing models using AIC
#fl = AIC table with a row for each model
NEfl <- fitList.oSCR(list(NEm0,NEm_HFI,NEm_FC,NEm_pop,NEm_HFI_FC,NEm_pop_FC),
                   drop.asu=T, rename=TRUE) #rename=T adds sensible model names

#makes sure the parameter names are stored as factors
for (i in 1:length(NEfl)) {NEfl[[i]]$coef.mle$param <- as.factor(NEfl[[i]]$coef.mle$param)}

saveRDS(NEfl, "NEfl.rds")

#AIC table with weight and cumulative weight added
NEms<-modSel.oSCR(NEfl)

NE.AIC.table<-as.data.frame(NEms$aic.tab) #extract AIC table as data frame

#************Model-averaged parameter estimates*************
NEma<-ma.coef(NEms)       #model averaged parameter estimates
names(NEma)             # list of names in model averaged output table
                      #ma is an annoying list of vectors that as.data.frame doesn't work for

#manually make the data frame
NEModAvg<-data.frame(
  NE.Parameter       = NEma$Parameter,
  NE.Est             = NEma$Estimate,
  NE.SE              = NEma$`Std. Error`,
  NE.RVI             = NEma$RVI,
  stringsAsFactors = FALSE) 
#add Z, p, and CIs
NEModAvg <- NEModAvg %>% 
mutate(NE.zval = NE.Est / NE.SE,      
       NE.pval = 2 * (1 - pnorm(abs(NE.zval))),
       NE.lwr  = NE.Est - 1.96 * NE.SE,
       NE.upr  = NE.Est + 1.96 * NE.SE)
write.csv(NEModAvg, "NE_ModelAvgFull.csv")

##exp of intercept should be model-avg density, but it can be higher
##this is because the exponential of the average log (manual) is less than the average of the ##exponentials (ma.coeff function)
##manual is the right method to use, because the densities should be averaged, not the log coefficients
##in this case they are nearly identical
##The ma.coef() method is more convenient for coefficients, but for actual densities, averaging on the natural scale is more accurate
##convenient, easy method:
NE.MA_Dens<-exp(NEModAvg$NE.Est[1]) #model-averaged density
NE.MA_LCI<-exp(NEModAvg$NE.lwr[1]) #LCI of density
NE.MA_UCI<-exp(NEModAvg$NE.upr[1]) #UCI of density
NE.MAresults<-cbind(NE.MA_Dens,NE.MA_LCI,NE.MA_UCI)
NE.MAresults

##Manual model-averaged Density (less convenient method):

NEAllDens <- matrix(nrow = length(NEfl), ncol = 5)  ##matrix to store results in

##extract density, SE, CIs from each model
for (i in 1:length(NEfl)) {
  model.name<-names(NEfl)[i]
  D<-exp(NEfl[[i]]$outStats[4,"mle"]) ##density estimate
  se <- NEfl[[i]]$outStats[4, "std.er"]
  D_LCI<-exp(NEfl[[i]]$outStats[4,"mle"] - (1.96*NEfl[[i]]$outStats[4,"std.er"]))
  D_UCI<-exp(NEfl[[i]]$outStats[4,"mle"] + (1.96*NEfl[[i]]$outStats[4,"std.er"]))
NEAllDens[i,]<-cbind(model.name,D,se,D_LCI,D_UCI)
  }
colnames(NEAllDens)<-c("model","NE.D","NE.se","NE.D_LCI","NE.D_UCI") #assign column names
NEAllDens<-as.data.frame(NEAllDens) #make data frame
  #density, LCI, UCI stored as characters, change to numeric
NEAllDens[2:5] <- lapply(NEAllDens[2:5], as.numeric)

  #join density estimates with AIC table, multiply model weights by density and CIs
NEDens_AIC<-left_join(NEAllDens, NE.AIC.table, by = "model") 
NEDens_AIC <- NEDens_AIC %>%
  mutate(NE.wD = D*weight, 
         NE.wD_LCI = D_LCI*weight, 
         NE.wD_UCI = D_UCI*weight)
write.csv(NEDens_AIC, "NE_DensityAICFull.csv")

##sum across models for model-averaged density and CIs
summarise_at(NEDens_AIC, c("NE.wD","NE.wD_LCI","NE.wD_UCI"), sum)

##spatial density predictions from each model

NEpredDens <- matrix(nrow = nrow(NEmask), ncol = length(NEfl))  ##matrix to store results in
colnames(NEpredDens) <- names(NEfl)

# store full get.real output 
NEpredList <- vector("list", length(NEfl))
names(NEpredList) <- names(NEfl)

##get pixel-level density for each model, takes a while
for (i in seq_along(NEfl)) {
  
  pred_i <- get.real(
    model   = NEfl[[i]],
    type    = "dens",
    newdata = NEmask,   # use mask/grid for pixel-level predictions
    d.factor = 1
  )
  
  NEpredList[[i]] <- pred_i
  NEpredDens[, i] <- pred_i$estimate
}
NEpredDens <- as.data.frame(NEpredDens)
saveRDS(NEpredList, "NEpredList.rds")

# extract and align model weights
mod_wts <- NE.AIC.table$weight
names(mod_wts) <- NE.AIC.table$model
mod_wts <- mod_wts[colnames(NEpredDens)]

# model-averaged density for each pixel
NEmask$coyD <- as.numeric(as.matrix(NEpredDens) %*% mod_wts)
write.csv(NEmask, "NEpred.csv")

#########OKANOGAN STUDY AREA MODELS#################
#import capture file
OKedf = read.csv("MV_S_capth.csv")
head(OKedf)

#import trap file, 1 for each session (1 session per summer)
OKtdf1 = read.csv("MV_18_S_trap.csv")
OKtdf2 = read.csv("MV_19_S_trap.csv")
OKtdf3 = read.csv("MV_20_S_trap.csv")

#rescale coordinates
OKtdf1$X <- OKtdf1$X/1000; OKtdf1$Y <- OKtdf1$Y/1000
OKtdf2$X <- OKtdf2$X/1000; OKtdf2$Y <- OKtdf2$Y/1000
OKtdf3$X <- OKtdf3$X/1000; OKtdf3$Y <- OKtdf3$Y/1000

#import mask file with 3 covariates: human footprint index (HFI, 0-10,000), 
#forest cover (FC, 0-1), and human population density (pop, #/km2)
OKmask = read.csv("MV_mask.csv")
#rescale coordinates
OKmask$X <- OKmask$X/1000; OKmask$Y <- OKmask$Y/1000

#inspect covariates
names(OKmask)
hist(OKmask$HFI) ##ranges from 0-10000, 0-inflated
hist(log(OKmask$HFI + 1)) #a bit skewed but better
hist(OKmask$pop) #0-200, very skewed
hist(log(OKmask$pop+1)) #skewed towards 0, logging helps but still skewed
hist(OKmask$FC) #0-1, fairly uniform
hist(asin(sqrt(OKmask$FC))) #more bell curved

##add scaled and transformed variables to mask dataframe
OKmask$n.HFI <-scale(log(OKmask$HFI + 1))
OKmask$n.FC <-scale(asin(sqrt(OKmask$FC)))
OKmask$n.pop <-scale(log(OKmask$pop+1))

#check for correlations among covariates
OKcovar <- select(OKmask, n.HFI, n.FC, n.pop)
OKCorCovar<-cor(OKcovar)  #correlation matrix of covariates
OKCorCovar
##n.pop vs n.HFI r = 0.37, pop vs FC r = -0.28, FC vs HFI = 0.11

##create mask dataframe for each session (same for all 3)
OKssDF <- list(data.frame(OKmask), data.frame(OKmask), data.frame(OKmask))

# create the scrFrame using data2oscr
OKsf <- data2oscr(edf = OKedf, #the EDF
                  sess.col = 1, #session column
                  id.col = 2, #individual column
                  occ.col = 3, #occasion column
                  trap.col = 4, #trap column
                  tdf = list(OKtdf1,OKtdf2, OKtdf3), #list of TDFs,
                  K = c(1,1,1), #makes counts of scats binary in each cell
                  trapcov.names = c("Effort"), #covariate names
                  tdf.sep = "/",
                  ntraps = c(nrow(OKtdf1),nrow(OKtdf2),nrow(OKtdf3)) #no. traps vector
)$scrFrame
class(OKssDF) = "ssDF"
str(OKssDF)

plot(OKsf,jit = 2)
plot(ssDF = OKssDF, scrFrame = OKsf)


##Models with binary "clog" encounters, and log Effort as covariate on p0
##takes a while to run 

OKm0 <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                 trimS=OKsf$mdm,
                 model=list(D ~ 1,p0 ~ log(Effort),sig ~ 1))

OKm_HFI <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                    trimS=OKsf$mdm,
                    model=list(D~n.HFI,p0~log(Effort),sig~1))

OKm_FC <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                   trimS=OKsf$mdm,
                   model=list(D~n.FC,p0~log(Effort),sig~1))

OKm_pop <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                    trimS=OKsf$mdm,
                    model=list(D~n.pop,p0~log(Effort),sig~1))

OKm_HFI_pop_FC <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                           trimS=OKsf$mdm,
                           model=list(D~n.FC+n.HFI+n.pop,p0~log(Effort),sig~1))

OKm_HFI_pop <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                        trimS=OKsf$mdm,
                        model=list(D~n.HFI+n.pop,p0~log(Effort),sig~1))

OKm_HFI_FC <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                       trimS=OKsf$mdm,
                       model=list(D~n.FC+n.HFI,p0~log(Effort),sig~1))

OKm_pop_FC <- oSCR.fit(scrFrame=OKsf,ssDF=OKssDF,DorN="D",encmod="CLOG",
                       trimS=OKsf$mdm,
                       model=list(D~n.FC+n.pop,p0~log(Effort),sig~1))

#Comparing models using AIC
#fl = AIC table with a row for each model
OKfl <- fitList.oSCR(list(OKm0,OKm_HFI,OKm_FC,OKm_pop,OKm_HFI_pop_FC,OKm_HFI_pop,
                          OKm_HFI_FC,OKm_pop_FC),
                   drop.asu=T, rename=TRUE) #rename=T adds sensible model names

#makes sure the parameter names are stored as factors
for (i in 1:length(OKfl)) {OKfl[[i]]$coef.mle$param <- as.factor(OKfl[[i]]$coef.mle$param)}

saveRDS(OKfl, "OKfl.rds")

#AIC table with weight and cumulative weight added
OKms<-modSel.oSCR(OKfl)

OK.AIC.table<-as.data.frame(OKms$aic.tab) #extract AIC table as data frame

#************Model-averaged parameter estimates*************
OKma<-ma.coef(OKms)       #model averaged parameter estimates
names(OKma)             # list of names in model averaged output table
#ma is an annoying list of vectors that as.data.frame doesn't work for

#manually make the data frame
OKModAvg<-data.frame(
  OK.Parameter       = OKma$Parameter,
  OK.Est             = OKma$Estimate,
  OK.SE              = OKma$`Std. Error`,
  OK.RVI             = OKma$RVI,
  stringsAsFactors = FALSE) 
#add Z, p, and CIs
OKModAvg <- OKModAvg %>% 
  mutate(OK.zval = OK.Est / OK.SE,      
         OK.pval = 2 * (1 - pnorm(abs(OK.zval))),
         OK.lwr  = OK.Est - 1.96 * OK.SE,
         OK.upr  = OK.Est + 1.96 * OK.SE)
write.csv(OKModAvg, "OK_ModelAvgFull.csv")

##exp of intercept should be model-avg density, but it can be higher
##this is because the exponential of the average log (manual) is less than the average of the ##exponentials (ma.coeff function)
##manual is the right method to use, because the densities should be averaged, not the log coefficients
##in this case they are nearly identical
##The ma.coef() method is more convenient for coefficients, but for actual densities, averaging on the natural scale is more accurate
##convenient, easy method:
OK.MA_Dens<-exp(OKModAvg$OK.Est[1]) #model-averaged density
OK.MA_LCI<-exp(OKModAvg$OK.lwr[1]) #LCI of density
OK.MA_UCI<-exp(OKModAvg$OK.upr[1]) #UCI of density
OK.MAresults<-cbind(OK.MA_Dens,OK.MA_LCI,OK.MA_UCI)
OK.MAresults

##Manual model-averaged Density (less convenient, easy method):

OKAllDens <- matrix(nrow = length(OKfl), ncol = 5)  ##matrix to store results in

##extract density, SE, CIs from each model
for (i in 1:length(OKfl)) {
  model.name<-names(OKfl)[i]
  D<-exp(OKfl[[i]]$outStats[4,"mle"]) ##density estimate
  se <- OKfl[[i]]$outStats[4, "std.er"]
  D_LCI<-exp(OKfl[[i]]$outStats[4,"mle"] - (1.96*OKfl[[i]]$outStats[4,"std.er"]))
  D_UCI<-exp(OKfl[[i]]$outStats[4,"mle"] + (1.96*OKfl[[i]]$outStats[4,"std.er"]))
  OKAllDens[i,]<-cbind(model.name,D,se,D_LCI,D_UCI)
}
colnames(OKAllDens)<-c("model","OK.D","OK.se","OK.D_LCI","OK.D_UCI") #assign column names
OKAllDens<-as.data.frame(OKAllDens) #make data frame
#density, LCI, UCI stored as characters, change to numeric
OKAllDens[2:5] <- lapply(OKAllDens[2:5], as.numeric)

#join density estimates with AIC table, multiply model weights by density and CIs
OKDens_AIC<-left_join(OKAllDens, OK.AIC.table, by = "model") 
OKDens_AIC <- OKDens_AIC %>%
  mutate(OK.wD = OK.D*weight, 
         OK.wD_LCI = OK.D_LCI*weight, 
         OK.wD_UCI = OK.D_UCI*weight)
write.csv(OKDens_AIC, "OK_DensityAICFull.csv")

##sum across models for model-averaged density and CIs
summarise_at(OKDens_AIC, c("OK.wD","OK.wD_LCI","OK.wD_UCI"), sum)

##spatial density predictions from each model

OKpredDens <- matrix(nrow = nrow(OKmask), ncol = length(OKfl))  ##matrix to store results in
colnames(OKpredDens) <- names(OKfl)

# store full get.real output 
OKpredList <- vector("list", length(OKfl))
names(OKpredList) <- names(OKfl)

#extracts pixel-level density for each model, takes a while
for (i in seq_along(OKfl)) {
    pred_i <- get.real(
    model   = OKfl[[i]],
    type    = "dens",
    newdata = OKmask,   # use mask/grid for pixel-level predictions
    d.factor = 1
  )
  
  OKpredList[[i]] <- pred_i
  OKpredDens[, i] <- pred_i$estimate
}
OKpredDens <- as.data.frame(OKpredDens)

# extract and align model weights
OKmod_wts <- OK.AIC.table$weight
names(OKmod_wts) <- OK.AIC.table$model
OKmod_wts <- OKmod_wts[colnames(OKpredDens)]

# model-averaged density for each pixel
OKmask$coyD <- as.numeric(as.matrix(OKpredDens) %*% OKmod_wts)
write.csv(OKmask, "OKpred.csv")

saveRDS(OKpredList, "OKpredList.rds")


