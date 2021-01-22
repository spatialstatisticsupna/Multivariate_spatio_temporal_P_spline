################################################################################
## Title: Multivariate Bayesian spatio-temporal P-spline models to analyse    ##
##        crimes against women                                                ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T.- Ugarte, M.D.                            ##
##                                                                            ##
## doi:                                                                       ##
##                                                                            ##
################################################################################
## Fit multivariate spatio-temporal P-spline models                           ##
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
## One of 'gaussian' (RECOMMENDED), 'simplified.laplace' (default), 'laplace' or 'adaptive'
strategy <- "simplified.laplace"


## Function to run several multivariate Bayesian spatio-temporal P-spline models
source("./functions/mps_st.R")

## Function to fit several multivariate Bayesian spatio-temporal P-spline models
source("functions/generic_mps.R")

################################################################################
## 1) Fitting multivariate spatio-temporal P-spline models                    ##
##    with temporal correlations                                              ##
################################################################################

##################################################
## 1.1) Fitting model:                          ##
##      Spatial prior: RW1                      ##
##      Temporal prior: RW1                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW1.RW1.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                     crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=1,
                     prior.temporal=1,
                     prior.interaction=2,
                     temp.corre=TRUE)

##################################################
## 1.2) Fitting model:                          ##
##      Spatial prior: RW1                      ##
##      Temporal prior: RW2                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW1.RW2.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                     crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=1,
                     prior.temporal=2,
                     prior.interaction=2,
                     temp.corre=TRUE)

##################################################
## 1.3) Fitting model:                          ##
##      Spatial prior: RW2                      ##
##      Temporal prior: RW1                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW2.RW1.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                     crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=2,
                     prior.temporal=1,
                     prior.interaction=2,
                     temp.corre=TRUE)

##################################################
## 1.4) Fitting model:                          ##
##      Spatial prior: RW2                      ##
##      Temporal prior: RW2                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW2.RW2.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                     crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=2,
                     prior.temporal=2,
                     prior.interaction=2,
                     temp.corre=TRUE)


##################################################
## Save                                         ##
##################################################
resulta.with <- list(RW1.RW1.T2=RW1.RW1.T2, RW1.RW2.T2=RW1.RW2.T2,
                     RW2.RW1.T2=RW2.RW1.T2, RW2.RW2.T2=RW2.RW2.T2)

save(resulta.with, file= "./resul/resulta_with.RData")

################################################################################
## 2) Fitting multivariate spatio-temporal P-spline models                    ##
##    without temporal correlations                                           ##
################################################################################

##################################################
## 2.1) Fitting model:                          ##
##      Spatial prior: RW1                      ##
##      Temporal prior: RW1                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW1.RW1.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year",
                     crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=1,
                     prior.temporal=1,
                     prior.interaction=2,
                     temp.corre=FALSE)

##################################################
## 2.2) Fitting model:                          ##
##      Spatial prior: RW1                      ##
##      Temporal prior: RW2                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW1.RW2.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year", crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=1,
                     prior.temporal=2,
                     prior.interaction=2,
                     temp.corre=FALSE)

##################################################
## 2.3) Fitting model:                          ##
##      Spatial prior: RW2                      ##
##      Temporal prior: RW1                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW2.RW1.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year", crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=2,
                     prior.temporal=1,
                     prior.interaction=2,
                     temp.corre=FALSE)

##################################################
## 2.4) Fitting model:                          ##
##      Spatial prior: RW2                      ##
##      Temporal prior: RW2                     ##
##      Spatio-temporal random effect: Type II  ##
##################################################
RW2.RW2.T2 <- mps_st(carto=carto, data=data, ID.area="ID_area", ID.year="ID_year", crimes=crimes, Expcrimes=paste0("e_",crimes),
                     centered=TRUE, strategy=strategy,
                     prior.spatial=2,
                     prior.temporal=2,
                     prior.interaction=2,
                     temp.corre=FALSE,
                     initial=TRUE)

##################################################
## Save                                         ##
##################################################
resulta.without <- list(RW1.RW1.T2=RW1.RW1.T2, RW1.RW2.T2=RW1.RW2.T2,
                        RW2.RW1.T2=RW2.RW1.T2, RW2.RW2.T2=RW2.RW2.T2)

save(resulta.without, file= "./resul/resulta_without.RData")
################################################################################
################################################################################