################################################################################
## Title: Multivariate Bayesian spatio-temporal P-spline models to analyse    ##
##        crimes against women                                                ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T.- Ugarte, M.D.                            ##
##                                                                            ##
## doi:                                                                       ##
##                                                                            ##
################################################################################
## Fit univariate spatio-temporal P-spline models                             ##
################################################################################
rm(list=ls())

## libraries
library(INLA); library(splines); library(spdep)

## Folder to save results
if(!file.exists("resul")) {dir.create("resul")}

## data loading
load("data_Psplines.RData")

## Variables' names (crimes)
crimes <- c("rape", "assault", "cruelty", "kidnapping")


## The strategy to use for the approximations
## One of 'gaussian', 'simplified.laplace' (default), 'laplace' or 'adaptive'
strategy <- "simplified.laplace"


## Function to run several univariate Bayesian spatio-temporal P-spline models
source("./functions/ps_st.R")

################################################################################
##  1) Fitting univariate spatio-temporal P-spline models                     ##
################################################################################

##################################################
## 1.1) Fitting model:                          ##
##      Spatial prior: RW1/RW2                  ##
##      Temporal prior: RW1/RW2                 ##
##      Spatio-temporal random effect: Type I   ##
##################################################

RW1.RW1.T1 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes),
                    centered=TRUE, strategy=strategy,
                    prior.spatial=1,
                    prior.temporal=1,
                    prior.interaction=1)

RW1.RW2.T1 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes),
                    centered=TRUE, strategy=strategy,
                    prior.spatial=1,
                    prior.temporal=2,
                    prior.interaction=1)

RW2.RW1.T1 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes),
                    centered=TRUE, strategy=strategy,
                    prior.spatial=2,
                    prior.temporal=1,
                    prior.interaction=1)

RW2.RW2.T1 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes),
                    centered=TRUE, strategy=strategy,
                    prior.spatial=2,
                    prior.temporal=2,
                    prior.interaction=1)

## save
resulta.T1 <- list(RW1.RW1.T1=RW1.RW1.T1, RW1.RW2.T1=RW1.RW2.T1,
                   RW2.RW1.T1=RW2.RW1.T1, RW2.RW2.T1=RW2.RW2.T1)

save(resulta.T1, file= "./resul/resulta_univariate_T1.RData")


##################################################
## 1.2) Fitting model:                          ##
##      Spatial prior: RW1/RW2                  ##
##      Temporal prior: RW1/RW2                 ##
##      Spatio-temporal random effect: Type II  ##
##################################################

RW1.RW1.T2 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=1, prior.temporal=1, prior.interaction=2,
                    strategy=strategy)

RW1.RW2.T2 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=1, prior.temporal=2, prior.interaction=2,
                    strategy=strategy)

RW2.RW1.T2 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=2, prior.temporal=1, prior.interaction=2, 
                    strategy=strategy)

RW2.RW2.T2 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=2, prior.temporal=2, prior.interaction=2,
                    strategy=strategy)

## save
resulta.T2 <- list(RW1.RW1.T2=RW1.RW1.T2, RW1.RW2.T2=RW1.RW2.T2,
                   RW2.RW1.T2=RW2.RW1.T2, RW2.RW2.T2=RW2.RW2.T2)

save(resulta.T2, file= "./resul/resulta_univariate_T2.RData")


##################################################
## 1.3) Fitting model:                          ##
##      Spatial prior: RW1/RW2                  ##
##      Temporal prior: RW1/RW2                 ##
##      Spatio-temporal random effect: Type III ##
##################################################

RW1.RW1.T3 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=1, prior.temporal=1, prior.interaction=3, 
                    strategy=strategy)

RW1.RW2.T3 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=1, prior.temporal=2, prior.interaction=3, 
                    strategy=strategy)

RW2.RW1.T3 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=2, prior.temporal=1, prior.interaction=3, 
                    strategy=strategy)

RW2.RW2.T3 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=2, prior.temporal=2, prior.interaction=3, 
                    strategy=strategy)

## save
resulta.T3 <- list(RW1.RW1.T3=RW1.RW1.T3, RW1.RW2.T3=RW1.RW2.T3,
                   RW2.RW1.T3=RW2.RW1.T3, RW2.RW2.T3=RW2.RW2.T3)

save(resulta.T3, file= "./resul/resulta_univariate_T3.RData")

##################################################
## 1.4) Fitting model:                          ##
##      Spatial prior: RW1/RW2                  ##
##      Temporal prior: RW1/RW2                 ##
##      Spatio-temporal random effect: Type IV  ##
##################################################

RW1.RW1.T4 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=1, prior.temporal=1, prior.interaction=4, 
                    strategy=strategy)

RW1.RW2.T4 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=1, prior.temporal=2, prior.interaction=4, 
                    strategy=strategy)

RW2.RW1.T4 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=2, prior.temporal=1, prior.interaction=4, 
                    strategy=strategy)

RW2.RW2.T4 <- ps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                    crimes=crimes, Expcrimes=paste0("e_",crimes), centered=TRUE,
                    prior.spatial=2, prior.temporal=2, prior.interaction=4, 
                    strategy=strategy)

## save
resulta.T4 <- list(RW1.RW1.T4=RW1.RW1.T4, RW1.RW2.T4=RW1.RW2.T4,
                   RW2.RW1.T4=RW2.RW1.T4, RW2.RW2.T4=RW2.RW2.T4)

save(resulta.T4, file= "./resul/resulta_univariate_T4.RData")

################################################################################
################################################################################