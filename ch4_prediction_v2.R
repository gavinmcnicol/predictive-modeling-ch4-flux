# Gavin McNicol 
# Half-hourly fluxes machine learning
# This script mimics Bodesheim et al. 2018 (FLUXCOM) for CH4 prediction

# packages
library(lubridate)
library(tidyverse)
library(caret)
library(gbm)
library(Metrics)
library(ranger)
library(RColorBrewer)
library(officer)
library(svMisc)
library(tictoc)

# clear workspace
rm(list=ls())

# location/working directory 
data_loc <- "C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data"
model_output <- "C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_analysis"
modis_loc <- "C:/Users/gavin/Box/MODIS Data/AppEEARS Data"
worldclim_loc <- "C:/Users/gavin/Box/Upscaling Resources/WorldClim Data"

##
# load fluxes
fluxes <- read.csv(paste(data_loc,"/Upscaling_analysis/","flat_fluxes_v2.csv",sep=""))

# create a copy to reload from
fluxes0 <- fluxes
# fluxes <- fluxes0

# get site.names
site.names <- fluxes %>% select(ID) %>% distinct() %>%  pull()

# get metadata
metadata <- read_csv(paste(data_loc,"/Upscaling_analysis/190201_sites.csv",sep="")) %>% 
  mutate(ID = substr(site.names, 1,5),
         Class = site.class,
         Biome = site.biome,
         PF = site.pf,
         LAT = site.lat,
         LONG = site.long) %>% 
  select(ID, Class, Biome, PF, LAT, LONG)


# load geospatial predictor data 
modis <- read.csv(paste(modis_loc,"/V2.0_MODIS.csv",sep=""))
modis2 <- read.csv(paste(modis_loc,"/v2.0_MODIS_20152018.csv",sep=""))

# Convert Date into Year/DOY 
modis <- modis %>% 
  mutate(DATE = as.Date(DATE, format = "%Y-%m-%d")) %>% 
  mutate(ID = substr(ID, 1,5),
         DOY = yday(DATE),
         DATE = as_date(DATE),
         Month = as.numeric(substr(DATE, 6,7)),
         Year = as.integer(substr(DATE,1,4)),
         SRWI = b02/b05) %>% 
  select(ID, LAT, LONG ,DATE, Year, Month, DOY, LSTD, LSTN, EVI, LAI, SRWI) 

modis2 <- modis2 %>% 
  mutate(DATE = as.Date(DATE, format = "%Y-%m-%d")) %>% 
  mutate(ID = substr(ID, 1,5),
         DOY = yday(DATE),
         DATE = as_date(DATE),
         Month = as.numeric(substr(DATE, 6,7)),
         Year = as.integer(substr(DATE,1,4)),
         SRWI = b02/b05) %>% 
  select(ID, LAT, LONG ,DATE, Year, Month, DOY, LSTD, LSTN, EVI, LAI, SRWI) 

# join
modis <- modis %>% bind_rows(modis2)

# Check data coverage
modis %>% group_by(ID, Year) %>% 
  filter(ID == "BCBog" & Year > 2010) %>% 
  summarize(records1 = sum(!is.na(LSTD)),
            records2 = sum(!is.na(LSTN)),
            records3 = sum(!is.na(EVI)),
            records4 = sum(!is.na(SRWI)))

# get monthly means
modis <- modis %>% 
  group_by(ID, Year, Month) %>% 
  summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  select(ID, Year, Month, LSTD, LSTN, EVI, LAI, SRWI)

# join metadata
fluxes <- fluxes %>% left_join(metadata, by = "ID") %>% 
  mutate(index = 1:n()) %>% 
  select(index, ID, LAT, LONG, Biome, Class, PF, Year, Month, DOY, Hour, FCH4, FCH4_F, FCH4_F_ANN, TS, TA, P, VPD, GPP_DT, NEE, LE, H)

# join modis
fluxes <- fluxes %>% left_join(modis, by = c("ID","Year","Month"))

# # explore sites for gap-filled trends
# fluxes %>% 
#   filter(ID == "USMyb") %>%
#   # filter(Year == 2015 & DOY == 200 ) %>% 
#   ggplot(aes(Hour, FCH4)) +
#   geom_point(col = 'black') +
#   geom_point(aes(Hour, FCH4_F_ANN), col = 'pink', alpha = 0.3) +
#   scale_y_continuous() +
#   facet_wrap(~Year, ncol = 1) +
#   theme_bw()

# filter long gaps and fill short gaps with ANN data
FCH4_tofill <- fluxes %>% 
  filter(!is.na(FCH4_F) & is.na(FCH4)) %>%
  select(index) %>% pull()
# fill with ANN
fluxes$FCH4[FCH4_tofill] <- fluxes$FCH4_F_ANN[FCH4_tofill]
# retain only gap filled data (2.5 million half hours)
fluxes <- fluxes %>% 
  filter(!is.na(FCH4_F)) %>% 
  select(-FCH4_F, -FCH4_F_ANN)

# save again as new version
fluxes1 <- fluxes

# remove non-wetland sites
fluxes <- fluxes %>% filter(Class %in% c("Bog","Fen","Marsh","Peat plateau","Swamp","Wet tundra"))

# remove firest year of USSne (restoring)
join.sne <- fluxes %>% 
  filter(ID == "USSne" & Year == 2016)
fluxes <- fluxes %>% setdiff(join.sne)

# filter any remaining missing FCH4 rows
fluxes <- fluxes %>% filter(!is.na(FCH4))

# quick half hour plots
fluxes %>%
  filter(ID == "USBeo") %>%
  ggplot(aes(DOY, FCH4)) +
  geom_point() +
  scale_y_continuous(limits=c(0,500)) 
  # facet_wrap(~Year, ncol = 1)

# get daily means
daily <- fluxes %>% 
  group_by(ID,LAT, LONG, Biome, Class,PF,Year,Month, DOY) %>% 
  summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  select(-Hour) 

# get monthly means
monthly <- fluxes %>% 
  group_by(ID,LAT, LONG, Biome, Class,PF,Year,Month) %>% 
  summarize_all(funs(mean(., na.rm = TRUE))) %>% 
  select(-Hour, -DOY) 

# get worldclim data (static climate variables)
worldclim <- read.csv(paste(worldclim_loc, "/V2.0_sites_worldclim.csv", sep=""))
worldclim <- worldclim %>% 
  mutate(ID = substr(ID,1,5))

# temp while select is covered
worldclim <- worldclim[,c(2,6:10)]

# join worldclim data
daily <- daily %>% left_join(worldclim, by = "ID") %>% select(-index)
monthly <- monthly %>% left_join(worldclim, by = "ID") %>% select(-index)

############# save csv of v.2 monthly data
# write.csv(monthly, "C:/Users/gavin/Box/Lab Meetings/190219 Methane Session/monthly_alldata.csv")
monthly <- read.csv("C:/Users/gavin/Box/Lab Meetings/190219 Methane Session/monthly_alldata.csv")

# quick daily plots
daily %>%
  # filter(ID %in% c("USORv","JPBBY","USStJ")) %>%
  ggplot(aes(DOY, LSTD))+
  geom_point() +
  scale_y_continuous(limits=c(250,350)) +
  facet_wrap(~ID, ncol = 7)

# quick monthly plots
monthly %>%
  # filter(ID %in% c("USORv","JPBBY","USStJ")) %>%
  ggplot(aes(Month, LAI))+
  geom_point() +
  scale_y_continuous(limits=c(0,3)) +
  facet_wrap(~ID, ncol = 7)


##
# start machine learning
train <- ungroup(monthly) 

# break into label and features
train_x <- train[,]
train_y <- train[,10]$FCH4

# redefine site.names, after wetland subsetting
site.names <- train %>% select("ID") %>% distinct() %>% pull()

# create a tibble where each site name has assocaited 'close by' sites
close.sites <- list()
close.sites[[1]] <- c("BCFEN","YFBsf") #BCBog
close.sites[[2]] <- c("BCBog","YFBsf") #BCFEN
close.sites[[3]] <- c("CASCC") #CASCB
close.sites[[4]] <- c("CASCB") #CASCC
close.sites[[5]] <- c("DESfN") #placeholder (i.e. no close sites, site removed as validation site)
close.sites[[6]] <- c("DEZrk") #placeholder
close.sites[[7]] <- c("FILom") #placeholder
close.sites[[8]] <- c("FISi2") #FISi1
close.sites[[9]] <- c("FISi1") #FISi2
close.sites[[10]] <- c("JPBBY") #placeholder

close.sites[[11]] <- c("MYMLM") #MYMLM
close.sites[[12]] <- c("NZKop") #placeholder
close.sites[[13]] <- c("RUChe") #RUCh2
close.sites[[14]] <- c("RUCh2") #RUChe
close.sites[[15]] <- c("RUSAM") #placeholder
close.sites[[16]] <- c("RUVrk") #placeholder
close.sites[[17]] <- c("SEDeg") #placeholder
close.sites[[18]] <- c("SESto") #SESt1
close.sites[[19]] <- c("SESt1") #SESto
close.sites[[20]] <- c("USAtq") #placeholder

close.sites[[21]] <- c("USBes", "USNGB","USBrw") # USBeo
close.sites[[22]] <- c("USBgl") #placeholder
close.sites[[23]] <- c("USFwm") #USBms
close.sites[[24]] <- c("USBms") #USFwm
close.sites[[25]] <- c("USIcs") #placeholder
close.sites[[26]] <- c("USIvo") #placeholder
close.sites[[27]] <- c("USLos") #placeholder
close.sites[[28]] <- c("USSne","USSnd","USBi1","USTwt","USBi2") #USMyb
close.sites[[29]] <- c("USNC4") #placeholder
close.sites[[30]] <- c("USBeo","USBes","USBrw") #USNGB

close.sites[[31]] <- c("USORv", "USOWC") #placeholder
close.sites[[32]] <- c("USWPT", "USCRT", "USORv") # USOWC
close.sites[[33]] <- c("USMyb","USSnd","USBi1","USBi2") #USSne
close.sites[[34]] <- c("USStj") #placeholder
close.sites[[35]] <- c("USSnd","USBi1","USTw4","USTwt","USBi2") #USTw1
close.sites[[36]] <- c("USSnd","USBi1","USTw1","USTwt","USBi2") #USTw4
close.sites[[37]] <- c("USOWC", "USCRT", "USORv") #USWPT

close.sites.list <- list()
close.sites.list <- as_tibble(cbind(sites = as.character(site.names), proximate = close.sites))

# create folds without proximate sites
folds_train <- list()
for (i in 1:length(site.names)) {
  folds_train[[i]] <- ungroup(train_x) %>%
    mutate(IDX = 1:n()) %>%
    filter(!ID %in% c(close.sites.list$proximate[[i]], as.character(site.names[i]))) %>%
    select("IDX") %>% pull()
}


## impute missing predictor values 
# select all predictor variables that can be used to preprocess predictors
train_x <- train[,c("PF","TA","TS","P","GPP_DT", "NEE", "LE","H", "VPD",
                    "Bio1","Bio3","Bio4","Bio12","Bio15",
                    "LAI","SRWI","EVI")]

# preprocess
pp <- preProcess(train_x, method = c("bagImpute"))

train_x_pp <- predict(pp, train_x)
train_x_pp$ID <- train$ID


## now subset rf predictors and train model
# select only predictors
train_x <- train_x_pp[,c("GPP_DT","LE","H","NEE","TS", # site level
                         "Bio1","Bio3","Bio4","Bio12","Bio15",  # static climate
                         "LAI","SRWI","EVI")] # modis (monthly)

# fill NAs in original training set
train[,c("PF","TA","TS","P","GPP_DT", "NEE", "LE","H", "VPD",
         "Bio1","Bio3","Bio4","Bio12","Bio15",
         "LAI","SRWI","EVI")] <- 
  train_x_pp[,c("PF","TA","TS","P","GPP_DT", "NEE", "LE","H", "VPD",
                "Bio1","Bio3","Bio4","Bio12","Bio15",
                "LAI","SRWI","EVI")]


# quick monthly plots
train %>%
  # filter(ID %in% c("USORv","JPBBY","USStJ")) %>%
  ggplot(aes(DOY, EVI))+
  geom_point() +
  scale_y_continuous(limits=c(0,0.75)) +
  facet_wrap(~ID, ncol = 7)

## set up lists (for rf)
tgrid <- list()
myControl <- list()

## Create tune-grid
tgrid <- expand.grid(
  .mtry = c(2:13),
  .splitrule = "variance", 
  .min.node.size = c(2,5,10,20)
)

## Create trainControl object
myControl <- trainControl(
  method = "oob",
  classProbs = FALSE,
  verboseIter = TRUE,
  savePredictions = TRUE,
  index = folds_train
)

## train rf on folds
rf_model <- list()
for (i in 1:length(site.names)){
  rf_model[[i]] <- train(
    x = train_x[folds_train[[i]],], 
    y = train_y[folds_train[[i]]],
    method = 'ranger',
    trControl = myControl,
    tuneGrid =tgrid,
    num.trees = 300
  )
  print(i)
}

## save/output model structure
saveRDS(rf_model, paste(model_output,"/V2.0/190214_rf_wetlands_monthly_no_LST.rds",sep=""))

## gbm - not working great 
# ## set up lists (gbm)
# tgrid <- list()
# myControl <- list()
# 
# 
# ## Create tune-grid
# tgrid <- expand.grid(
#   .interaction.depth = c(8,10,12),
#   .n.trees = 500, 
#   .shrinkage = c(0.05,0.03),
#   .n.minobsinnode = 10
# )
# 
# ## Create trainControl object
# myControl <- trainControl(
#   method = "cv",
#   number = 10
# )
# 
# # set metric
# metric <- "RMSE"
# 
# ## train gbm on folds
# gbm_model <- list()
# for (i in 1:length(site.names)){
#   gbm_model[[i]] <- train(
#     train_x[folds_train[[i]],],
#     train_y[folds_train[[i]]], 
#     method = 'gbm',
#     trControl = myControl,
#     tuneGrid = tgrid,
#      metric = metric
#   )
#   print(i)
# }
# 
# ## save/output model structure
# saveRDS(gbm_model, paste(model_output,"/190213_gbm_models.rds",sep=""))



# load older models
# rf_model <- readRDS(paste(model_output,"/V2.0/190214_rf_wetlands_monthly_no_LST.rds",sep=""))



# look at r2 for all models
x <- c()
for (i in 1:length(rf_model)){
 x[1] <- max(rf_model[[i]]$results$Rsquared)
}

summary(x)

# get all predictions
rf.pred <- list()
for (i in 1:length(rf_model)) {
  rf.pred[[i]] <- ungroup(train) %>% 
    filter(ID == site.names[i]) %>%   
    mutate(FCH4P = predict(rf_model[[i]], .))
}
rf.pred.all <- bind_rows(rf.pred)

### plotting predictions
# plot all site-years
rf.pred.all %>% 
  filter(ID %in% c("BCBog","USBeo","FISi1","USSne","RUChe","USMyb","USORv","USOWC","NZKop")) %>% 
  ggplot(aes(Month, FCH4))+
  geom_point(size = 2, col = "grey") +
  # scale_y_continuous(limits = c(0,300))+
  geom_point(aes(Month, FCH4P), col = 'orange', size = 2, alpha = 0.8)+
  facet_wrap(~ID, ncol = 3, scale = 'free') +
theme_bw() +
  theme(panel.border = element_blank(), 
        axis.title=element_text(size=14), axis.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x= 'Month', y = expression(CH[4]*' Flux (nmol m'^{-2}*' s'^{-1}*')')) +
  theme(strip.text = element_text(face="bold", size=8),
        strip.background = element_rect(fill=cols, colour=cols,size=1))

# global R2 
ungroup(rf.pred.all) %>% 
  group_by(Class) %>% 
  filter(!ID %in% c("USORv", "USOWC")) %>% 
  summarize(Rsquared = summary(lm(FCH4P ~ FCH4))$adj.r.squared,
            NSE = 1 - sum((FCH4 - FCH4P)^2) / sum((FCH4 - mean(FCH4))^2),
            Mean_res = mean(FCH4 - FCH4P),
            SD_res = sd(FCH4 - FCH4P),
            MedianO = median(FCH4),
            sdO = sd(FCH4),
            MedianP = median(FCH4P),
            sdP = sd(FCH4P),
            samples = n()) %>% 
  arrange(desc(Rsquared)) %>%   View()
  filter(Rsquared < 0.3) %>%
  select(sdP) %>% 
  pull() %>% median()


  rf.pred.all %>% 
  group_by(ID, Year, Month) %>% 
    summarize(Biome = Biome[1]) %>% 
  select(Biome) %>% 
  pull() %>% 
  as.factor() %>% 
  summary()

  # all site pred vs. obs
rf.pred.all %>% 
  # filter(Class == "Marsh") %>% 
  filter(!ID %in% c("USORv", "USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P)) +
  geom_point() +
  scale_x_continuous(limits = c(0,150)) +
  scale_y_continuous(limits = c(0,150)) +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) +
  facet_wrap(~ID, ncol = 6)

# by site predictions
rf.pred.all %>% 
  filter(ID == "JPBBY") %>% 
  ggplot(aes(DOY, FCH4))+
  geom_point() +
  geom_point(aes(DOY, FCH4P), col = 'green')+
  facet_wrap(~ID, ncol = 6) 

# summaies by month
summary.ch4 <- rf.pred.all %>%
  group_by(ID, Month) %>%
  summarize(Class = Class[1],
            Biome = Biome[1],
            MeanO = mean(FCH4),
            MeanP = mean(FCH4P),
            MedianO = median(FCH4),
            MedianP = median(FCH4P),
            Rsquared = cor(FCH4, FCH4P)^2,
            NSE = 1 - sum((FCH4 - FCH4P)^2) / sum((FCH4 - mean(FCH4))^2),
            Mean_res = mean(FCH4 - FCH4P),
            SD_res = sd(FCH4 - FCH4P),
            samp = n()) %>%
  arrange(desc(NSE)) 
View(summary.ch4)

# all months pred vs obs
summary.ch4 %>% 
  filter(!ID %in% c("USOWC", "USORv")) %>% 
  ggplot(aes(MedianO, MedianP, col  =Biome), na.rm = TRUE) + 
  # geom_smooth(method = "lm", se = F, na.rm = TRUE)+
  geom_point(alpha = 0.8, na.rm = TRUE, size = 2) + 
  # geom_point(alpha = 0.9, na.rm = TRUE, aes(MeanO, MeanP, col = Biome)) +
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.major = element_blank(), 
        panel.background = element_blank()) +
  scale_x_continuous(limits = c(1,300)) +
  scale_y_continuous(limits = c(1,300)) +
  labs(x = "Validation Data", y = "RF Predictions") +
  geom_abline(slope = 1, intercept = log10(1),
                            na.rm = FALSE, show.legend = NA)

# pred vs. obs by biome (removing OWORv and OWC)
rf.pred.all %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P, col = Biome)) +
  geom_point(alpha = 0.6) + 
  scale_x_log10(limits = c(1,800)) +
  scale_y_log10(limits = c(1,800)) +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) +
  facet_wrap(~Biome)

# pred vs. obs by class (removing OWORv and OWC)
rf.pred.all %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P, col = Biome)) +
  geom_point(alpha = 0.6) + 
  scale_x_log10(limits = c(1,800)) +
  scale_y_log10(limits = c(1,800)) +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) +
  facet_wrap(~Class)

# plot performance by site-year
plot <- list()
plot.pred.obs <- list()
for (i in 1:length(rf_model)){
  plot[[i]] <- rf.pred.all %>%
    filter(ID == site.names[i]) %>% 
    ggplot(aes(Month, FCH4, col = "Observed")) + 
    geom_point() + 
    geom_point(aes(Month, FCH4P, col = "Predicted")) +
    facet_wrap(~ID + Year, ncol = 6) +
    scale_y_continuous(limits = c(0,600)) +
    theme(legend.position = "bottom")
  
  plot.pred.obs[[i]] <- rf.pred[[i]] %>% 
    ggplot(aes(FCH4, FCH4P, col = Year)) + 
    geom_point() + 
    facet_wrap(~ID, ncol = 3, scale = "free")
}

#  all data, pred vs. obs with 1:1 line
rf.pred.all %>% 
  # group_by(ID) %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P), na.rm = TRUE) + 
  geom_smooth(method = "lm", se = F, na.rm = TRUE)+
  geom_point(alpha = 0.5, na.rm = TRUE) + 
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.major = element_blank(), 
        panel.background = element_blank()) +
  scale_x_log10(limits = c(1,1000)) +
  scale_y_log10(limits = c(1,1000)) +
  labs(x = "Validation Data", y = "RF Predictions") +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) 

# plot histograms of pred. vs. obs by Class 
rf.pred.all %>% 
  group_by(Class) %>% 
  ggplot(aes(FCH4, fill = 'class'),  na.rm = TRUE) +
  geom_histogram(fill = 'black', alpha = 0.8) +
  geom_histogram(aes(FCH4P), alpha = 0.5, na.rm = TRUE) +
  facet_wrap(~Class) +
  scale_x_log10(limits = c(1,1000)) +
  theme(legend.position = 'bottom')


# create list of names
biome_names <- list(
  '1'="Tundra",
  '2'="Boreal",
  '3'="Temperate",
  '4'="(Sub)Tropical"
)

# create a labeller function
biome_labeller <- function(variable,value){
  return(biome_names[value])
}

# plot histograms of pred. vs. obs by Biome
rf.pred.all %>% 
  mutate(Biome = as.factor(Biome)) %>% 
  mutate(BiomeRank = recode_factor(Biome, Tundra = "1", Boreal = "2", Temperate = "3", Tropical = "4")) %>% 
  group_by(Biome) %>% 
  ggplot(aes(FCH4),  na.rm = TRUE) +
  geom_histogram(fill = 'BLACK', alpha = 0.8) +
  geom_histogram(aes(FCH4P, fill = 'red'), alpha = 0.5, na.rm = TRUE) +
  facet_wrap(~BiomeRank, ncol = 1, labeller = biome_labeller) +
  scale_x_log10(limits = c(1,1000)) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        legend.position = 'none',
        axis.title=element_text(size=16), axis.text=element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x=expression(CH[4]*' Flux (nmol m'^{-2}*' s'^{-1}*')'), y = 'Count (Site-Months)') +
  theme(strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill=cols, colour=cols,size=1))


# plot pred versus temperature
rf.pred.all %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(TS, FCH4)) + 
  geom_point(alpha = 0.3, na.rm = TRUE) + 
  theme_minimal() +
  scale_x_continuous(limits = c(-10,30)) +
  scale_y_continuous(limits = c(0,750)) +
  geom_point(aes(TS, FCH4P), na.rm= TRUE, col = 'blue', alpha = 0.6) 



#######  IF trained on daily fluxes
# aggregate daily by month
rf.pred.month <- rf.pred.all %>% 
  group_by(ID, Year, Month) %>% 
  summarize_all(funs(mean(., na.rm = TRUE)))

# plot all site-years
rf.pred.month %>% 
  # filter(ID == "DEZrk") %>% 
  ggplot(aes(DOY, FCH4))+
  geom_point() +
  # scale_y_continuous(limits = c(0,300))+
  geom_point(aes(DOY, FCH4P), col = 'green')+
  facet_wrap(~ID, ncol = 6, scale = 'free') 

# global R2 
ungroup(rf.pred.month) %>% 
  # group_by(ID) %>% 
  summarize(Rsquared = summary(lm(FCH4P ~ FCH4))$adj.r.squared,
            NSE = 1 - sum((FCH4 - FCH4P)^2) / sum((FCH4 - mean(FCH4))^2),
            Mean_res = mean(FCH4 - FCH4P),
            SD_res = sd(FCH4 - FCH4P),
            MedianO = median(FCH4),
            sdO = sd(FCH4),
            MedianP = median(FCH4P),
            sdP = sd(FCH4P),
            samples = n()) %>% 
  arrange(desc(Rsquared)) %>% 
  View()

# all site pred vs. obs
rf.pred.month %>% 
  # filter(Class == "Marsh") %>% 
  ggplot(aes(FCH4, FCH4P)) +
  geom_point() +
  scale_x_continuous(limits = c(0,150)) +
  scale_y_continuous(limits = c(0,150)) +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) +
  facet_wrap(~ID, ncol = 6)

# by site predictions
rf.pred.month %>% 
  filter(ID == "JPBBY") %>% 
  ggplot(aes(DOY, FCH4))+
  geom_point() +
  geom_point(aes(DOY, FCH4P), col = 'green')+
  facet_wrap(~ID, ncol = 6) 

# summaies by month
summary.ch4 <- rf.pred.month %>%
  group_by(ID, Month) %>%
  summarize(Class = Class[1],
            Biome = Biome[1],
            MeanO = mean(FCH4),
            MeanP = mean(FCH4P),
            MedianO = median(FCH4),
            MedianP = median(FCH4P),
            Rsquared = cor(FCH4, FCH4P)^2,
            NSE = 1 - sum((FCH4 - FCH4P)^2) / sum((FCH4 - mean(FCH4))^2),
            Mean_res = mean(FCH4 - FCH4P),
            SD_res = sd(FCH4 - FCH4P),
            samp = n()) %>%
  arrange(desc(NSE)) 

# all months pred vs obs
summary.ch4 %>% 
  ggplot(aes(MedianO, MedianP, col  =Biome), na.rm = TRUE) + 
  # geom_smooth(method = "lm", se = F, na.rm = TRUE)+
  geom_point(alpha = 0.8, na.rm = TRUE, size = 2) + 
  # geom_point(alpha = 0.9, na.rm = TRUE, aes(MeanO, MeanP, col = Biome)) +
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.major = element_blank(), 
        panel.background = element_blank()) +
  scale_x_log10(limits = c(1,1000)) +
  scale_y_log10(limits = c(1,1000)) +
  labs(x = "Validation Data", y = "RF Predictions") +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA)

# pred vs. obs by biome (removing OWORv and OWC)
rf.pred.month %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P, col = Biome)) +
  geom_point(alpha = 0.6) + 
  scale_x_log10(limits = c(1,800)) +
  scale_y_log10(limits = c(1,800)) +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) +
  facet_wrap(~Biome)

# pred vs. obs by class (removing OWORv and OWC)
rf.pred.month %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P, col = Biome)) +
  geom_point(alpha = 0.6) + 
  scale_x_log10(limits = c(1,800)) +
  scale_y_log10(limits = c(1,800)) +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) +
  facet_wrap(~Class)

# plot performance by site-year
plot <- list()
plot.pred.obs <- list()
for (i in 1:length(rf_model)){
  plot[[i]] <- rf.pred.month %>%
    filter(ID == site.names[i]) %>% 
    ggplot(aes(Month, FCH4, col = "Observed")) + 
    geom_point() + 
    geom_point(aes(Month, FCH4P, col = "Predicted")) +
    facet_wrap(~ID + Year, ncol = 6) +
    scale_y_continuous(limits = c(0,600)) +
    theme(legend.position = "bottom")
  
  plot.pred.obs[[i]] <- rf.pred[[i]] %>% 
    ggplot(aes(FCH4, FCH4P, col = Year)) + 
    geom_point() + 
    facet_wrap(~ID, ncol = 3, scale = "free")
}

#  all data, pred vs. obs with 1:1 line
rf.pred.month %>% 
  # group_by(ID) %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(FCH4, FCH4P), na.rm = TRUE) + 
  geom_smooth(method = "lm", se = F, na.rm = TRUE)+
  geom_point(alpha = 0.5, na.rm = TRUE) + 
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.major = element_blank(), 
        panel.background = element_blank()) +
  scale_x_log10(limits = c(1,1000)) +
  scale_y_log10(limits = c(1,1000)) +
  labs(x = "Validation Data", y = "RF Predictions") +
  geom_abline(slope = 1, intercept = log10(1),
              na.rm = FALSE, show.legend = NA) 

# plot histograms of pred. vs. obs by Class 
rf.pred.month %>% 
  group_by(Class) %>% 
  ggplot(aes(FCH4, fill = 'class'),  na.rm = TRUE) +
  geom_histogram(fill = 'black', alpha = 0.8) +
  geom_histogram(aes(FCH4P), alpha = 0.5, na.rm = TRUE) +
  facet_wrap(~Class) +
  scale_x_log10(limits = c(1,1000)) +
  theme(legend.position = 'bottom')

# plot histograms of pred. vs. obs by Biome
rf.pred.month %>% 
  group_by(Biome) %>% 
  ggplot(aes(FCH4, fill = 'class'),  na.rm = TRUE) +
  geom_histogram(fill = 'black', alpha = 0.8) +
  geom_histogram(aes(FCH4P), alpha = 0.5, na.rm = TRUE) +
  facet_wrap(~Biome) +
  scale_x_log10(limits = c(1,1000)) +
  theme(legend.position = 'bottom')

# plot pred versus temperature
rf.pred.month %>% 
  filter(!ID %in% c("USORv","USOWC")) %>% 
  ggplot(aes(TS, FCH4)) + 
  geom_point(alpha = 0.3, na.rm = TRUE) + 
  theme_minimal() +
  scale_x_continuous(limits = c(-10,30)) +
  scale_y_continuous(limits = c(0,750)) +
  geom_point(aes(TS, FCH4P), na.rm= TRUE, col = 'blue', alpha = 0.6) 
