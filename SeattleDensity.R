library(oSCR)
library(dplyr)
library(AICcmodavg)
library(ggplot2)
library(corrplot)

#import capture file
SEedf = read.csv("Seattle_capt.csv")
head(SEedf)

#import trap file
SEtdf = read.csv("Seattle_traps.csv")

#rescale coordinates
SEtdf$X <- SEtdf$X/1000; SEtdf$Y <- SEtdf$Y/1000

#import mask file with covariates
SEmask = read.csv("Seattle_mask.csv")
#rescale coordinates
SEmask$X <- SEmask$X/1000; SEmask$Y <- SEmask$Y/1000

#create environmental contaminant covariate by averaging the 4 metrics
SEmask$env <- rowMeans(SEmask[,c("Lead", "PNPL", "PRMP", "TSDF", "Diesel", "Traffic", "PM25", "Pesticide")])

#inspect covariates
names(SEmask)
hist(SEmask$env) #ranges from 4-9, just scale
hist(SEmask$pop_dens)
hist(log(SEmask$pop_dens+1)) #skewed, logging helps
hist(SEmask$tree_cov/100)
hist(asin(sqrt((SEmask$tree_cov/100)))) #0-1, 0-inflated, asin transform helps a little
hist(SEmask$impervious/100)  #0-1, fairly normal
hist(asin(sqrt((SEmask$impervious/100))))  #asin transform similar, just scale
hist(SEmask$income) #fairly normal, just scale

#add scaled/transformed covariates to mask file
SEmask$n.env <-scale(SEmask$env)
SEmask$n.FC <-scale(asin(sqrt((SEmask$tree_cov/100))))
SEmask$n.pop <-scale(log(SEmask$pop_dens+1))
SEmask$n.IS <- scale(SEmask$impervious/100)
SEmask$n.inc <- scale(SEmask$income)

#check for correlations among covariates
SEcovar <- select(SEmask, n.env, n.FC, n.pop, n.IS, n.inc)
SECorCovar<-cor(SEcovar)  #correlation matrix of covariates
SECorCovar
# all are low except imperv surface and forest cover (-0.88)

##create mask dataframe 
SEssDF <- list(data.frame(SEmask))

# create the scrFrame using data2oscr
SEsf <- data2oscr(edf = SEedf, #the EDF
                sess.col = 1, #session column
                id.col = 2, #individual column
                occ.col = 3, #occasion column
                trap.col = 4, #trap column
                tdf = list(SEtdf), #list of TDFs
                K = c(1), # occasions
                trapcov.names = c("Effort"), #covariate names
                tdf.sep = "/",
                ntraps = c(nrow(SEtdf)), #no. traps vector
)$scrFrame

class(SEssDF) <- "ssDF"

str(SEssDF)

plot(SEsf,jit = 2)
plot(ssDF = SEssDF, scrFrame = SEsf)

##construct and run model set, with constraint that FC and IS cannot be in
##same model because they are too correlated

#normalized covariates
covs <- c("n.pop", "n.IS", "n.FC", "n.inc", "n.env")

make_D_formula <- function(vars) {
  if (length(vars) == 0) return(as.formula("D ~ 1"))
  as.formula(paste("D ~", paste(vars, collapse = " + ")))
}

all_subsets <- unlist(
  lapply(0:length(covs), function(k) combn(covs, k, simplify = FALSE)),
  recursive = FALSE
)

#remove models that have both IS and FC
valid_subsets <- Filter(function(v) !("n.IS" %in% v && "n.FC" %in% v), all_subsets)
D_forms <- lapply(valid_subsets, make_D_formula)

# fixed parts
p0_form  <- p0 ~ log(Effort)
sig_form <- sig ~ 1

length(D_forms)  # number of models in set (24)

fits <- vector("list", length(D_forms))
names(fits) <- sprintf("m%02d_%s",
                       seq_along(D_forms),
                       sapply(D_forms, function(f) gsub("[[:space:]]+", "", deparse(f))))

#loop to run all 24 models - this should take a while
for (i in seq_along(D_forms)) {
  fits[[i]] <- oSCR.fit(
    scrFrame = SEsf,
    ssDF     = SEssDF,
    DorN     = "D",
    encmod   = "CLOG",
    trimS    = SEsf$mdm,
    model    = list(D_forms[[i]], p0_form, sig_form)
  )
  cat("Finished", i, "of", length(D_forms), ":", names(fits)[i], "\n")
}

fl <- fitList.oSCR(fits)

#makes sure the parameter names are stored as factors
for (i in 1:length(fl)) {fl[[i]]$coef.mle$param <- as.factor(fl[[i]]$coef.mle$param)}

#AIC table with weight and cumulative weight added
ms<-modSel.oSCR(fl)

AIC.table<-as.data.frame(ms$aic.tab) #extract AIC table as data frame

saveRDS(fl, "SEfl.rds")

#************Model-averaged parameter estimates*************
ma<-ma.coef(ms)       #model averaged parameter estimates
names(ma)             # list of names in model averaged output table
                      #ma is an annoying list of vectors that as.data.frame doesn't work for

#manually make the data frame
ModAvg<-data.frame(
  Parameter       = ma$Parameter,
  Est             = ma$Estimate,
  SE              = ma$`Std. Error`,
  RVI             = ma$RVI,
  stringsAsFactors = FALSE) 
#add Z, p, and CIs
ModAvg <- ModAvg %>% 
mutate(zval = Est / SE,      
       pval = 2 * (1 - pnorm(abs(zval))),
       lwr  = Est - 1.96 * SE,
       upr  = Est + 1.96 * SE)
write.csv(ModAvg, "SE_ModelAvg.csv")

##exp of intercept should be model-avg density, but it can be higher
##this is because the exponential of the average log (manual) is less than the average of the ##exponentials (ma.coeff function)
##manual is the right method to use, because the densities should be averaged, not the log coefficients
##in this case they are nearly identical
##The ma.coef() method is more convenient for coefficients, but for actual densities, averaging on the natural scale is more accurate
##convenient, easy method:

##scale density from 250x250m to km2
cell<- ((250*250)/(1000*1000)) #cell (pixel) area in km2

MA_Dens<-exp(ModAvg$Est[1])/cell #model-averaged density in km2
MA_LCI<-exp(ModAvg$lwr[1])/cell #LCI of density
MA_UCI<-exp(ModAvg$upr[1])/cell #UCI of density
MA_results<-cbind(MA_Dens,MA_LCI,MA_UCI) 
MA_results

##Manual model-averaged Density (less convenient method):

AllDens <- matrix(nrow = length(fl), ncol = 5)  ##matrix to store results in

##extract density, SE, CIs from each model
for (i in 1:length(fl)) {
  model.name<-names(fl)[i]
  D<-exp(fl[[i]]$outStats[4,"mle"])/cell ##density estimate
  se <- fl[[i]]$outStats[4, "std.er"]
  D_LCI<-(exp(fl[[i]]$outStats[4,"mle"] - (1.96*fl[[i]]$outStats[4,"std.er"])))/cell
  D_UCI<-(exp(fl[[i]]$outStats[4,"mle"] + (1.96*fl[[i]]$outStats[4,"std.er"])))/cell
AllDens[i,]<-cbind(model.name,D,se,D_LCI,D_UCI)
  }
colnames(AllDens)<-c("model","D","se","D_LCI","D_UCI") #assign column names
AllDens<-as.data.frame(AllDens) #make data frame
  #density, LCI, UCI stored as characters, change to numeric
AllDens[2:5] <- lapply(AllDens[2:5], as.numeric)

  #join density estimates with AIC table, multiply model weights by density and CIs
Dens_AIC<-left_join(AllDens, AIC.table, by = "model") 
Dens_AIC <- Dens_AIC %>%
  mutate(wD = D*weight, 
         wD_LCI = D_LCI*weight, 
         wD_UCI = D_UCI*weight)
write.csv(Dens_AIC, "SE_DensityAIC_New.csv")

##sum across models for model-averaged density and CIs
summarise_at(Dens_AIC, c("wD","wD_LCI","wD_UCI"), sum)


#############################################################
#make predictions on the real scale to map model-averaged density

predDens <- matrix(nrow = nrow(SEmask), ncol = length(fl))  ##matrix to store results in
colnames(predDens) <- names(fl)

# store full get.real output 
predList <- vector("list", length(fl))
names(predList) <- names(fl)

##extracts pixel-level predicted density for each model in the fit list (fl)
##takes a while

for (i in seq_along(fl)) {
  
  pred_i <- get.real(
    model   = fl[[i]],
    type    = "dens",
    newdata = SEmask,   # use mask/grid for pixel-level predictions
    d.factor = 1/cell
  )
  
  predList[[i]] <- pred_i
  predDens[, i] <- pred_i$estimate
}
predDens <- as.data.frame(predDens)

saveRDS(predList, "predList.rds")

# extract and align model weights
mod_wts <- AIC.table$weight
names(mod_wts) <- AIC.table$model
mod_wts <- mod_wts[colnames(predDens)]

# model-averaged density for each pixel
SEmask$coyD <- as.numeric(as.matrix(predDens) %*% mod_wts)
write.csv(SEmask, "SEpred.csv")
