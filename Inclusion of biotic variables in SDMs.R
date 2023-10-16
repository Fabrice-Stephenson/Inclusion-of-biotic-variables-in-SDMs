##==========================      BOOTSTRAPPED BOOSTED REGRESSION TREE MODELS   =======================================##

##--------------------------  Authors: Fabrice Stephenson
##--------------------------  Aim: Example code used to explore possible influence of biotic interactions in 
##--------------------------  explaining patterns of species' distribution, abundance and explanatory power. Analysis  
##--------------------------  described in the manuscript "Inclusion of biotic variables improve predictions of 
##--------------------------  environmental niche models."
##--------------------------  Project : COME2001 and COME2101
##--------------------------  Start date : 01/05/2020
##--------------------------  End date : 1/12/2021

##======================================================================================================================##

####=======================         PACKAGES                 =========================================================####
# Load library packages
library(dismo); library(dplyr); library(gbm); library(plyr)
library(corrplot); library(snow); library(ggplot2);library(e1071)
library(devEMF); library(raster); library(rgdal); library(spatialEco)

####=======================         LOAD FILES               =========================================================####
setwd("~")
# LOAD BIOLOGICAL DATA 

load("BI.Rdata")

# A dataframe called BI with 853 rows (sampling sites) and 16 columns (Label: unique identified assigned by Kraan et al., 2015;
# Region: Estuary data was collected in; environmental variables: "Seagrass","Shallhash_cover","Chla","Phaephytin"
# "DOI", "Mud"; Biotic variables: "Ausstu_J"(juvenile Astrovenus count),"Ausstu"(adult Astrovenus count), 
# "Maclil_J" (juvenile Macomona count), "Maclil"(Adult Macomona count), "Maclil_comb"(juvenile and adult Macomona count)
# "Ausstu_comb" (juvenile and adult Austrovenus count), "Mac_patch" (70 IDW estimate of macomona count),
# "Ausstu_patch" (70 IDW estimate of Austrovenus count)


####=======================         MODEL TUNING               =======================================================####

# The code below presents the analysis for Macomona. To apply this code to Austrovenus, the biotic variables used 
# (dependant and independant) would need to be changed as described in the manuscript "Inclusion of biotic variables improve 
# predictions of environmental niche models"

env.all <- c("Seagrass","Shallhash_cover","Chla","Phaephytin","DOI", "Mud", "Region")

M1 <- gbm.step(data = BI, gbm.x = env.all, 
               gbm.y = "Maclil_comb", family = "gaussian", tree.complexity = 3,
               learning.rate = 0.001, bag.fraction = 0.6, n.folds=10, 
               max.trees = 5000, plot.main = F, tolerance.method = "auto",
               tolerance = 0.001, verbose = F, step.size = 25, silent = T)

gbm.simplify(M1, n.folds = 10, n.drops = "auto", alpha = 1, prev.stratify = T, plot = T)

####=======================         ENVIRONMENT ONLY MODEL    =======================================================####
n.boot <- 100
imp.vars <- sort(c("Seagrass", "Shellhash_cover", "Chla","Pheophytin", "DOI", "Mud"))

# Array for saving model fit metrics
deviance_mat <- array(0, c(n.boot,4))
colnames(deviance_mat) <- c("Dev.Train", "Dev.Eval", "Cor.Train", "Cor.Eval")

# Array for saving information on predictor importance
influence_mat <- array(0,c(length(imp.vars)+1,n.boot))
rownames(influence_mat) <- c("Chla", "DOI", "Mud", "Pheophytin", 'Region', 'Seagrass','Shellhash_cover')

# create environmental gradient files used for prediction of partial dependance plots
PredMins <- apply(BI[,imp.vars],2,min)
PredMaxs <- apply(BI[,imp.vars],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(imp.vars))))
names(EnvRanges) <- imp.vars
for (i in c(1:length(imp.vars))) {
     EnvRanges[,i] <- seq(PredMins[imp.vars[i]],PredMaxs[imp.vars[i]],length = 100)
}

# create 3D array for saving partial dependence plots of each env predictor from each bootstrap
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot))
dimnames(PD)[[2]] <- imp.vars
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

Region_PD <-array(0,c(3,n.boot))
dimnames(Region_PD)[[1]] <- c("Manukau", "Tauranga", "Kaipara")

# MAIN FOR LOOP
for (i in 1:n.boot){
     # split data into training and eval data (stratified by estuary = Region in DF)
     train <- ddply(BI, ~Region, function(x){ndf <- x[sample(nrow(x),nrow(x)), ]})
     eval <- BI[!BI$Label  %in% train$Label,] # Select absences based on unique identifier
     
     # RUN BRT
     M1 <- gbm.step(data = train, gbm.x = c(2:8), gbm.y = 13,
                    family = "gaussian", tree.complexity = 3,
                    learning.rate = 0.001, bag.fraction = 0.6, n.folds=10, 
                    max.trees = 5000, plot.main = F, tolerance.method = "auto",
                    tolerance = 0.001, verbose = F, step.size = 25,silent = T)
     
     # Saving partial dependence plots
     for (j in 1:length(imp.vars)){
          grid <- gbm::plot.gbm(M1, i.var = c(paste(imp.vars[[j]])), return.grid = T, type = c("response"))
          PD[,imp.vars[[j]],i] <- loess(grid$y~EnvRanges[,imp.vars[[j]]])$y
          Region_PD[,i] <- gbm::plot.gbm(M1, i.var = c("Region"), return.grid = T)$y
     }
     
     # MODEL FIT METRICS
     # Deviance explained - training data
     int.null.deviance <- M1$self.statistics$mean.null 
     int.residual.deviance <- M1$cv.statistics$deviance.mean 
     deviance_mat[i,"Dev.Train"] <- (int.null.deviance-int.residual.deviance)/int.null.deviance
     
     # Deviance explained - evaluation data
     pred <- predict.gbm(M1, eval, n.trees = M1$gbm.call$best.trees, type = "response") 
     ext.residual.deviance <- calc.deviance(eval$Maclil_comb, pred, family = "gaussian" ,calc.mean=T) 
     ext.null.deviance <- calc.deviance(eval$Maclil_comb, family = "gaussian", rep(mean(eval$Maclil_comb),nrow(eval)), calc.mean=T) 
     deviance_mat[i,"Dev.Eval"] <-(ext.null.deviance - ext.residual.deviance)/ext.null.deviance
     
     # Correlation - training data
     deviance_mat[i,"Cor.Train"] <- M1$self.statistics$correlation
     # Correlation - evaluation data
     deviance_mat[i,"Cor.Eval"] <- cor(pred, eval$Maclil_comb)
     
     # Predictor variable importance 
     M1_contrib <- as.data.frame(M1$contributions)
     env_var_ord <- M1_contrib[order(M1_contrib$var),]
     influence_mat[,i]<-env_var_ord[,"rel.inf"]
     
     print(paste("Iteration ", i, " out of ", n.boot, sep = ""))
}

# model fits
M.Fits <-round(apply(deviance_mat, 2, function(x) c(Mean = mean(x), SD = sd(x))), 2)
write.csv(M.Fits, file = "MAC_ENV.csv")

# variable importance +- SD
pred_inf <-t(round(apply(influence_mat, 1, function(x) c(Mean = mean(x), SD = sd(x))), 1))

x.label <- sort(c("Seagrass", "Shellhash cover", "Chla",
                  "Pheophytin", "LOI", "Mud"))

# plot PDs + 95 PI and save to file
emf(file = "PD_Mac_ENV.emf", emfPlus = T)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(3,3))

# change this value to fit all env range or change to have order of var importance
for (i in c(1:7)) {
     if (i < 7){
          plot(EnvRanges[,i],apply(expm1(PD[,i,]),1,mean), col = "black",type='l',
               xlab = paste(x.label[i], " (",pred_inf[imp.vars[i],1], "% ± ",pred_inf[imp.vars[i],2],")", sep = ""), 
               ylab = '',
               ylim = c(min(expm1(PD[,,])), max(expm1(PD[,,])))) #max(boot_array_EnvTran[,,])))
          # 95% PI
          UC <- na.omit(cbind(EnvRanges[,imp.vars[i]],
                              apply(expm1(PD[,i,]),1, quantile, probs= c(0.05)),
                              apply(expm1(PD[,i,]),1, quantile, probs= c(0.95))))
          polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
          rug(quantile(BI[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
     } else {boxplot(t(expm1(Region_PD)), 
                     boxwex=0.5,
                     ylim = c(min(expm1(PD[,,])), max(expm1(PD[,,]))))
          mtext(paste("Region (",pred_inf["Region",1], "% ± ",pred_inf["Region",2],")", sep = ""),
                side = 1, line = 3, cex = 0.75)
     }
}
dev.off()

####=======================         ENVIRONMENT + SPECIES MODEL    =======================================================####
n.boot <- 100
imp.vars <- sort(c("Seagrass", "Shellhash_cover", "Chla","Pheophytin", "DOI", "Mud", "Ausstu_comb"))

# Array for saving model fit metrics
deviance_mat <- array(0, c(n.boot,4))
colnames(deviance_mat) <- c("Dev.Train", "Dev.Eval", "Cor.Train", "Cor.Eval")

# Array for saving information on predictor importance
influence_mat <- array(0,c(length(imp.vars)+1,n.boot))
rownames(influence_mat) <- sort(c("Ausstu_comb", "Chla", "DOI", "Mud", "Pheophytin","Region", "Seagrass","Shellhash_cover"))

# create environmental gradient files used for prediction of partial dependance plots
PredMins <- apply(BI[,imp.vars],2,min)
PredMaxs <- apply(BI[,imp.vars],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(imp.vars))))
names(EnvRanges) <- imp.vars
for (i in c(1:length(imp.vars))) {
     EnvRanges[,i] <- seq(PredMins[imp.vars[i]],PredMaxs[imp.vars[i]],length = 100)
}

# create 3D array for saving partial dependence plots of each env predictor from each bootstrap
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot))
dimnames(PD)[[2]] <- imp.vars
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

Region_PD <-array(0,c(3,n.boot))
dimnames(Region_PD)[[1]] <- c("Manukau", "Tauranga", "Kaipara")

# MAIN FOR LOOP
for (i in 1:n.boot){
     # split data into training and eval data (stratified by estuary = Region in DF)
     train <- ddply(BI, ~Region, function(x){ndf <- x[sample(nrow(x),nrow(x)), ]})
     eval <- BI[!BI$Label  %in% train$Label,] # Select absences based on unique identifier
     
     # RUN BRT
     M1 <- gbm.step(data = train, gbm.x = c(2:8,14), gbm.y = 13,
                    family = "gaussian", tree.complexity = 3,
                    learning.rate = 0.001, bag.fraction = 0.6, n.folds=10, 
                    max.trees = 5000, plot.main = F, tolerance.method = "auto",
                    tolerance = 0.001, verbose = F, step.size = 25,silent = T)
     
     # Saving partial dependence plots
     for (j in 1:length(imp.vars)){
          grid <- gbm::plot.gbm(M1, i.var = c(paste(imp.vars[[j]])), return.grid = T, type = c("response"))
          PD[,imp.vars[[j]],i] <- loess(grid$y~EnvRanges[,imp.vars[[j]]])$y
          Region_PD[,i] <- gbm::plot.gbm(M1, i.var = c("Region"), return.grid = T)$y
     }
     
     # MODEL FIT METRICS
     # Deviance explained - training data
     int.null.deviance <- M1$self.statistics$mean.null 
     int.residual.deviance <- M1$cv.statistics$deviance.mean 
     deviance_mat[i,"Dev.Train"] <- (int.null.deviance-int.residual.deviance)/int.null.deviance
     
     # Deviance explained - evaluation data
     pred <- predict.gbm(M1, eval, n.trees = M1$gbm.call$best.trees, type = "response") 
     ext.residual.deviance <- calc.deviance(eval$Maclil_comb, pred, family = "gaussian" ,calc.mean=T) 
     ext.null.deviance <- calc.deviance(eval$Maclil_comb, family = "gaussian", rep(mean(eval$Maclil_comb),nrow(eval)), calc.mean=T) 
     deviance_mat[i,"Dev.Eval"] <-(ext.null.deviance - ext.residual.deviance)/ext.null.deviance
     
     # Correlation - training data
     deviance_mat[i,"Cor.Train"] <- M1$self.statistics$correlation
     # Correlation - evaluation data
     deviance_mat[i,"Cor.Eval"] <- cor(pred, eval$Maclil_comb)
     
     # Predictor variable importance 
     M1_contrib <- as.data.frame(M1$contributions)
     env_var_ord <- M1_contrib[order(M1_contrib$var),]
     influence_mat[,i]<-env_var_ord[,"rel.inf"]
     
     print(paste("Iteration ", i, " out of ", n.boot, sep = ""))
}

# model fits
M.Fits <-round(apply(deviance_mat, 2, function(x) c(Mean = mean(x), SD = sd(x))), 2)
write.csv(M.Fits, file = "MAC_ENV.csv")

# variable importance +- SD
pred_inf <-t(round(apply(influence_mat, 1, function(x) c(Mean = mean(x), SD = sd(x))), 1))

x.label <- sort(c("Seagrass", "Shellhash cover", "Chla",
                  "Pheophytin", "LOI", "Mud","Austrovenus"))

# plot PDs + 95 PI and save to file
emf(file = "PD_Mac_ENV_AUS.emf", emfPlus = T)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(3,3))

# change this value to fit all env range or change to have order of var importance
for (i in c(1:8)) {
     if (i < 8){
          plot(EnvRanges[,i],apply(expm1(PD[,i,]),1,mean), col = "black",type='l',
               xlab = paste(a.label[i], " (",pred_inf[imp.vars[i],1], "% ± ",pred_inf[imp.vars[i],2],")", sep = ""), 
               ylab = '',
               ylim = c(min(expm1(PD[,,])), max(expm1(PD[,,])))) #max(boot_array_EnvTran[,,])))
          # 95% PI
          UC <- na.omit(cbind(EnvRanges[,imp.vars[i]],
                              apply(expm1(PD[,i,]),1, quantile, probs= c(0.05)),
                              apply(expm1(PD[,i,]),1, quantile, probs= c(0.95))))
          polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
          rug(quantile(BI[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
     } else {boxplot(t(expm1(Region_PD)), 
                     boxwex=0.5,
                     ylim = c(min(expm1(PD[,,])), max(expm1(PD[,,]))))
          mtext(paste("Region (",pred_inf["Region",1], "% ± ",pred_inf["Region",2],")", sep = ""),
                side = 1, line = 3, cex = 0.75)
     }
}
dev.off()

####=======================         ENVIRONMENT + SPECIES + SPECIES PATCH MODEL  ======================================####
n.boot <- 100
imp.vars <- sort(c("Seagrass", "Shellhash_cover", "Chla","Pheophytin", "DOI", "Mud", "Ausstu_comb","Ausstu_patch"))

# Array for saving model fit metrics
deviance_mat <- array(0, c(n.boot,4))
colnames(deviance_mat) <- c("Dev.Train", "Dev.Eval", "Cor.Train", "Cor.Eval")

# Array for saving information on predictor importance
influence_mat <- array(0,c(length(imp.vars)+1,n.boot))
rownames(influence_mat) <- sort(c("Ausstu_comb", "Ausstu_patch","Chla", "DOI", "Mud", "Pheophytin","Region", "Seagrass",
                                  "Shellhash_cover"))

# create environmental gradient files used for prediction of partial dependence plots
PredMins <- apply(BI[,imp.vars],2,min)
PredMaxs <- apply(BI[,imp.vars],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(imp.vars))))
names(EnvRanges) <- imp.vars
for (i in c(1:length(imp.vars))) {
     EnvRanges[,i] <- seq(PredMins[imp.vars[i]],PredMaxs[imp.vars[i]],length = 100)
}

# create 3D array for saving partial dependence plots of each env predictor from each bootstrap
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot))
dimnames(PD)[[2]] <- imp.vars
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

Region_PD <-array(0,c(3,n.boot))
dimnames(Region_PD)[[1]] <- c("Manukau", "Tauranga", "Kaipara")

# MAIN FOR LOOP
for (i in 1:n.boot){
     # split data into training and eval data (stratified by estuary = Region in DF)
     train <- ddply(BI, ~Region, function(x){ndf <- x[sample(nrow(x),nrow(x)), ]})
     eval <- BI[!BI$Label  %in% train$Label,] # Select absences based on unique identifier
     
     # RUN BRT
     M1 <- gbm.step(data = train, gbm.x = c(2:8,14,16), gbm.y = 13,
                    family = "gaussian", tree.complexity = 3,
                    learning.rate = 0.001, bag.fraction = 0.6, n.folds=10, 
                    max.trees = 5000, plot.main = F, tolerance.method = "auto",
                    tolerance = 0.001, verbose = F, step.size = 25,silent = T)
     
     # Saving partial dependence plots
     for (j in 1:length(imp.vars)){
          grid <- gbm::plot.gbm(M1, i.var = c(paste(imp.vars[[j]])), return.grid = T, type = c("response"))
          PD[,imp.vars[[j]],i] <- loess(grid$y~EnvRanges[,imp.vars[[j]]])$y
          Region_PD[,i] <- gbm::plot.gbm(M1, i.var = c("Region"), return.grid = T)$y
     }
     
     # MODEL FIT METRICS
     # Deviance explained - training data
     int.null.deviance <- M1$self.statistics$mean.null 
     int.residual.deviance <- M1$cv.statistics$deviance.mean 
     deviance_mat[i,"Dev.Train"] <- (int.null.deviance-int.residual.deviance)/int.null.deviance
     
     # Deviance explained - evaluation data
     pred <- predict.gbm(M1, eval, n.trees = M1$gbm.call$best.trees, type = "response") 
     ext.residual.deviance <- calc.deviance(eval$Maclil_comb, pred, family = "gaussian" ,calc.mean=T) 
     ext.null.deviance <- calc.deviance(eval$Maclil_comb, family = "gaussian", rep(mean(eval$Maclil_comb),nrow(eval)), calc.mean=T) 
     deviance_mat[i,"Dev.Eval"] <-(ext.null.deviance - ext.residual.deviance)/ext.null.deviance
     
     # Correlation - training data
     deviance_mat[i,"Cor.Train"] <- M1$self.statistics$correlation
     # Correlation - evaluation data
     deviance_mat[i,"Cor.Eval"] <- cor(pred, eval$Maclil_comb)
     
     # Predictor variable importance 
     M1_contrib <- as.data.frame(M1$contributions)
     env_var_ord <- M1_contrib[order(M1_contrib$var),]
     influence_mat[,i]<-env_var_ord[,"rel.inf"]
     
     print(paste("Iteration ", i, " out of ", n.boot, sep = ""))
}

# model fits
M.Fits <-round(apply(deviance_mat, 2, function(x) c(Mean = mean(x), SD = sd(x))), 2)
write.csv(M.Fits, file = "MAC_ENV.csv")

# variable importance +- SD
pred_inf <-t(round(apply(influence_mat, 1, function(x) c(Mean = mean(x), SD = sd(x))), 1))

x.label <- sort(c("Seagrass", "Shellhash cover", "Chla",
                  "Pheophytin", "LOI", "Mud","Austrovenus","Austrovenus patch"))

# plot PDs + 95 PI and save to file
emf(file = "PD_MAC_ENV_AUS_PATCH.emf", emfPlus = T)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(3,3))

# change this value to fit all env range or change to have order of var importance
for (i in c(1:9)) {
     if (i < 9){
          plot(EnvRanges[,i],apply(expm1(PD[,i,]),1,mean), col = "black",type='l',
               xlab = paste(x.label[i], " (",pred_inf[imp.vars[i],1], "% ± ",pred_inf[imp.vars[i],2],")", sep = ""), 
               ylab = '',
               ylim = c(min(expm1(PD[,,])), max(expm1(PD[,,])))) #max(boot_array_EnvTran[,,])))
          # 95% PI
          UC <- na.omit(cbind(EnvRanges[,imp.vars[i]],
                              apply(expm1(PD[,i,]),1, quantile, probs= c(0.05)),
                              apply(expm1(PD[,i,]),1, quantile, probs= c(0.95))))
          polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
          rug(quantile(BI[,imp.vars[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
     } else {boxplot(t(expm1(Region_PD)), 
                     boxwex=0.5,
                     ylim = c(min(expm1(PD[,,])), max(expm1(PD[,,]))))
          mtext(paste("Region (",pred_inf["Region",1], "% ± ",pred_inf["Region",2],")", sep = ""),
                side = 1, line = 3, cex = 0.75)
     }
}
dev.off()
