


##############################################################
########################################  1.1. Data collection & Preprocessing
# install.packages("rgbif")
# install.packages("geodata")
# install.packages("terra")
# install.packages("usdm")
# install.packages("caret")
# install.packages("dismo")
# install.packages("MLmetrics“)
# install.packages("DALEX")

### Collection
# Download GBIF data from Lee et al. 2023
library(rgbif)
WCSB.occ <- occ_download_get('0428875-210914110416597') %>%  occ_download_import() 

# species data
library(terra)
DB.WCSB <- as.data.frame(WCSB.occ[c("decimalLongitude", "decimalLatitude")])
head(DB.WCSB)
WCSB.v <- vect(DB.WCSB, geom=c("decimalLongitude", "decimalLatitude"), crs="epsg:4326") 

# Collection of environment data (Env.)
library(geodata)
clim <- worldclim_global(var = 'bio', res = 10, path="GIS_data/" )

###  Preprocessing
#(1) Remove outlier (i.e. point in polar regions)
WCSB.ext <- ext(-180, 180, -60, 65) # set extent (exclude polar regions)
WCSB.v <- crop(WCSB.v, WCSB.ext)

#(2) Sampling bias (grid filtering)
r <- rast(clim) # create a SpatRaster 
r <- extend(r, ext(r)+1) # extend (expand) the extent of the SpatRaster a little
set.seed(123)
WCSB.p <- spatSample(WCSB.v, size=1, "random", strata=r) # only 1 site per grid

#(3) combined P sites and Env. data
WCSB.env <- extract(clim, WCSB.p)
values(WCSB.p) <- WCSB.env[-1] 

#(4) Remove P sites with NA data
WCSB.p <- WCSB.p[complete.cases(WCSB.env),]

#(5) Generate pesudo-absence sites with env. data
clim.masked <- mask(clim, WCSB.p, inverse = TRUE) # Area where P sites are not located 

set.seed(456)
WCSB.a <- spatSample(clim.masked, length(WCSB.p), "random", na.rm=TRUE, as.points=TRUE, ext=WCSB.ext) # Same number as 'used presence sites'

# check duplicated coordinates between P and A sites
sum(duplicated(rbind(geom(WCSB.p)[,3:4], geom(WCSB.a)[,3:4]))) # 0 : All coordinates are unique

# ### Visualization: global distribution of WCSB 
# world_map <- world(resolution = 3, path="GIS_data/" ) # World map
# par(mar=c(4,4,3,3))
# plot(world_map, col="lightyellow", border="lightgray", las=1)
# points(WCSB.a, pch=3, col="#780000", cex=0.1) # raw
# points(WCSB.p, pch=3, col="#03045e", cex=0.3) # raw

# Make Table (Data.frame)
WCSB.DB <- cbind( data.frame("OC" = c(rep(1, nrow(WCSB.p)), rep(0, nrow(WCSB.a)))), 
                  rbind( cbind(geom(WCSB.p)[,c("x", "y")], as.data.frame(WCSB.p)),
                         cbind(geom(WCSB.a)[,c("x", "y")], as.data.frame(WCSB.a))
                  ))

#(6) Solve Multicollinearity
library(usdm)
DB.env <- WCSB.DB[-c(1:3)]
ev1 <- vifcor(DB.env, th=0.9) # remove r>0.9 
t.DB.env <- exclude(DB.env, ev1)

ev2 <- vifstep(t.DB.env, th=5) # remove VIF > 5
f.DB.env <- exclude(t.DB.env, ev2) # final Env. DB

# Final DB for SDM
DB.for.SDM <- cbind(WCSB.DB[1], f.DB.env)

#(7) Data splitting (K-fold cross-validation)
# 5-folds; return values corresponding to 1 fold (1/5 = 20 %)
library(caret)
set.seed(123)
list.test.num <- createFolds(DB.for.SDM$OC, k = 5, list = T, returnTrain = FALSE) 

# Save data
save(WCSB.occ, WCSB.DB, DB.for.SDM, list.test.num, file="1.WCSB_All.RData") 
write.csv(WCSB.DB, file="1_1.DB_WCSB.csv", fileEncoding = "EUC-KR") 
write.csv(DB.for.SDM, file="1_2.DB_for_SDM.csv", fileEncoding = "EUC-KR") 

### Visualization: global distribution of WCSB  (ver 2)
world_map <- world(resolution = 3, path="GIS_data/" ) # World map
par(mar=c(4,4,3,3))
plot(world_map, col="lightyellow", border="lightgray", las=1)
grid()
points(WCSB.DB[c("x", "y")][WCSB.DB$OC == 0, ], pch=3, col="#780000", cex=0.1) # A
points(WCSB.DB[c("x", "y")][WCSB.DB$OC == 1, ], pch=3, col="#03045e", cex=0.3) # P

##############################################################
########################################  1.2. Develop the SVM
# install.packages("e1071")
load("1.WCSB_All.RData") 
library(e1071)

list.SVMs <- list()
DB.SVM.best.para <- c()
system.time({ # 3800s
for(k in 1:5){ #  k-folds
  # k=1
  # Data split
  t.train <- DB.for.SDM[-list.test.num[[k]], ] 
   
  # Grid search for optimal hyperparameter
  tuned_model <- tune.svm( OC ~ ., data = t.train, kernel = "radial", gamma = 2^(-3:1), cost = 2^(-1:1) ) # per 600s
  
  #summary(tuned_model)
  # Select best model
  list.SVMs[[k]] <- tuned_model$best.model
  
  # save model information 
  t.bp <- tuned_model$best.parameters
  t.bp$k <- k
  
  DB.SVM.best.para <-rbind(DB.SVM.best.para, t.bp)
}# for k
  
# Save data
save(list.SVMs, DB.SVM.best.para, file="2.SDM_DB.RData")  

})


##############################################################
########################################  1.3. Model evaluation 
library(e1071)
library(dismo)
library(MLmetrics)
library(caret)
library(DALEX)
load("1.WCSB_All.RData") 
load("2.SDM_DB.RData")

SVM.metrics <- c(); SVM.pred <- c(); SVM.VarImp <- c(); SVM.PDP <- c()
system.time({ # 240s 
for(k in 1:5){ 
  # Data split
  t.train<- DB.for.SDM[-list.test.num[[k]], ] 
  t.test <- DB.for.SDM[list.test.num[[k]], ] 
  t.SDM  <- list.SVMs[[k]]
  
  pred.v <- predict(t.SDM, t.test)
  eval.m <- evaluate(p=pred.v[t.test$OC == 1], a=pred.v[t.test$OC == 0])
  threhold.m <- as.numeric(threshold(eval.m)["spec_sens"])
  
  ### Evaluation metrics: Accuracy, AUROC, F1 score
  pred.OC <- ifelse(pred.v > threhold.m, 1, 0)
  
  m.AUC <- AUC(y_pred = pred.v, y_true = t.test$OC)
  m.ACC <- Accuracy(y_pred = pred.OC, y_true = t.test$OC)
  m.F1 <- F1_Score(y_true = t.test$OC, y_pred = pred.OC, positive = 1)
  
  # Create Model Explainer
  explainer.svm <- DALEX::explain(t.SDM, 
                                 data =t.test[,-1], 
                                 y = t.test[,1],
                                 label="svm")
  
  ### Permutation importance
  set.seed(123)
  t.var.imp <- model_parts(explainer.svm, N=NULL)
  t.var.imp <- t.var.imp[t.var.imp$permutation == 0, ] # extract average results
  
  ### Partial dependence plot
  t.pdp.svm  <- model_profile(explainer.svm, variable =names(t.test)[-1], type = "partial") # 12s

  plot(t.test[5])
  # Files to save
  DB.m.eval <- data.frame("k"=k, "Thr"=threhold.m, "ACC"=m.ACC, "AUC"=m.AUC, "F1"=m.F1)
  DB.m.pred <- data.frame("k"=k, "ID"=as.numeric(names(pred.v)), 
                          "Pred.value"=pred.v, "Pred.OC"=pred.OC, "Actual.OC"= pred.OC)
  DB.m.varimp <- data.frame("k"=k, "var.name"= t.var.imp$variable, "dropout.loss"= t.var.imp$dropout_loss )
  DB.m.PDP <- data.frame("k"=k, "var.name"= t.pdp.svm$agr_profiles$`_vname_`, 
                         "pdp.x"= t.pdp.svm$agr_profiles$`_x_`, "pdp.yhat"= t.pdp.svm$agr_profiles$`_yhat_`)
  
  SVM.metrics <- rbind(SVM.metrics, DB.m.eval)
  SVM.pred <- rbind(SVM.pred, DB.m.pred)
  SVM.VarImp <- rbind(SVM.VarImp, DB.m.varimp)
  SVM.PDP <- rbind(SVM.PDP, DB.m.PDP)
} # for k
  
  # Save data
  save(SVM.metrics, SVM.pred, SVM.VarImp, SVM.PDP, file="3.SDM_explain.RData")
})

plot(t.var.imp)

### Analysis (Metrics & VarImp) & Visualization (PDP)
load("3.SDM_explain.RData")
#save(SVM.metrics, SVM.pred, SVM.VarImp, SVM.PDP, file="3.SDM_explain.RData")

# Metrics
apply(SVM.metrics, 2, mean)
apply(SVM.metrics, 2, sd)

# Variable importnace
{
SVM.VarImp <- SVM.VarImp[!SVM.VarImp$var.name %in% c("_baseline_", "_full_model_"), ]

vi.mean <- aggregate(SVM.VarImp, . ~ var.name, mean) 
vi.sd <- aggregate(SVM.VarImp, . ~ var.name, sd) # SD

vi.order <- order(vi.mean$dropout.loss, decreasing = T)
vi.mean <- vi.mean[vi.order, ]
vi.sd <- vi.sd[vi.order, ]

#
colfunc <- colorRampPalette(c("#134074", "#8da9c4"))
vi.cols <- colfunc(8)

par(mar=c(4,4,3,2))
bp <- barplot(vi.mean$dropout.loss, col=vi.cols, ylim=c(0.2, 0.45), xpd=F,
              cex.names=1, las=1, ylab="Dropout Loss", xlab="", 
              names.arg=paste0("Bio", sapply(strsplit(vi.mean$var.name, "_"), function(x){ x[4] }))
              )
rect(bp-0.5, 0.2, bp+0.5, vi.mean$dropout.loss, col="NA", border="black", lwd=1.5)
arrows(x0 = bp, 
       y0 = vi.mean$dropout.loss - vi.sd$dropout.loss,
       y1 = vi.mean$dropout.loss + vi.sd$dropout.loss,
       code=3, length = 0.07, angle = 90, col="black", lwd=3)
}

# PDP
{
head(SVM.PDP)
u.var <- unique(SVM.PDP$var.name)
u.k <- unique(SVM.PDP$k)

for( i in 1:length(u.var)){ # number of variable 
  t.PDP <- SVM.PDP[SVM.PDP$var.name ==u.var[i], ] 
  par(mar=c(4,4,3,3))
  plot(t.PDP$pdp.x, t.PDP$pdp.yhat, type="n", las=1, xlab=paste0("Bio", strsplit(u.var[i], "_")[[1]][4]), ylab="Probability")
  grid()
  segments(x0=t.PDP$pdp.x, y0=par()$usr[3], y1=par()$usr[3]+(par()$usr[4]-par()$usr[3])/50, lwd=1, col="gray");box()
  for(k in 1:5){ # k-folds
    t.PDP2 <- t.PDP[t.PDP$k == k, ]
    lines(t.PDP2$pdp.x, t.PDP2$pdp.yhat, lwd=2, col="#073b4c")
  }# for k
}# for i
}

##############################################################
########################################  1.4. Prediction
library(e1071)
library(terra)

# load data (SVM model & Env. maps)
load("2.SDM_DB.RData")
cwd <- getwd()
setwd("./GIS_data/climate/wc2.1_10m/")
clim <- rast(list.files())
setwd(cwd)

# Spatial prediction
system.time({ # 660s
  svm.maps <- c()
  for(k in 1:5){
  t.map <- predict(clim, list.SVMs[[k]], core=5, na.rm=T)
  t.map <- ifel(t.map > 1, 1, t.map)
  t.map <- ifel(t.map < 0, 0, t.map)
  
  svm.maps <- c(svm.maps, t.map)
  }
  
  # stack
  svm.maps2 <- rast(svm.maps)
  
  # Mean calculation
  f.map <- app(svm.maps2, mean, cores=5, filename="4.SVM_maps_mean.tif")
})


### Visualization
library(terra)
f.map <- rast("4.SVM_maps_mean.tif")
# Global
plot(f.map, las=1)

# Europe
load("1.WCSB_All.RData") 
plot(f.map, xlim=c(-20,50), ylim=c(30, 60), las=1)
points(WCSB.DB[WCSB.DB$OC == 1, ][c("x", "y")], pch=3, col="red", cex=0.4)

# threshold
library(dismo)
load("3.SDM_explain.RData")
Ensemble.pred <- aggregate(SVM.pred, cbind(Pred.value, Actual.OC) ~ ID, mean)
en.thr <- evaluate(Ensemble.pred$Pred.value[Ensemble.pred$Actual.OC == 1], Ensemble.pred$Pred.value[Ensemble.pred$Actual.OC == 0])

f.map2 <- ifel(f.map > as.numeric(threshold(en.thr)["spec_sens"]), 1, 0)
plot(f.map2, las=1, col=c("gray40", "tomato"))

