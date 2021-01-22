# Functions

The [mps_st.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/functions/mps_st.R) file contains ```mps_st``` function.

The [ps_st.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/functions/ps_st.R) file contains ```ps_st``` function.

This [generic_mps.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/functions/generic_mps.R) file contains the ```generic_mps_spat``` and ```generic_mps_temp``` functions, necessary to fit multivariate P-spline models using INLA.


# Function ```mps_st```

## Description
Fit a multivariate Bayesian spatio-temporal P-spline models.

## Usage
```
mps_st(carto=NULL, data=NULL, ID.area=NULL, ID.year=NULL, crimes=NULL, Expcrimes=NULL,
       prior.spatial=1, prior.temporal=1, prior.interaction=2, temp.corre=TRUE,
       order.B=3, k.long=NULL, k.lati=NULL, k.time=NULL, centered=TRUE, initial=TRUE, 
       strategy="simplified.laplace", rerun=FALSE )
```

## Arguments

| Argument | Description |
|:---	| :--- | 
| ```carto```	| object of class ```SpatialPolygonsDataFrame``` or ```sf```. This object must contain at least the variable with the identifiers of the spatial areal units specified in  the argument ```ID.area```. | 
| ```data```	| object of class data.frame. This object must contain at least the	target variables of interest specified in the arguments ```ID.area```, ```ID.year```, ```crimes``` and ```Expcrimes```. |
| ```ID.area``` | character; name of the variable which contains the IDs of spatial areal units. |
| ```ID.year``` | character; name of the variable which contains the IDs of temporal units. |
| ```crimes``` | character; name of the variables that contain the observed number of crimes for each area and time unit. |
| ```Expcrimes```	| character; name of the variables that contain the expected number of crimes for each area and time unit. |
| ```prior.spatial```	| one of either 1 (default, for 'RW1') or 2 ('RW2'), which specifies the prior distribution considered on the coefficients of the spatial B-splines. |
| ```prior.temporal``` | one of either 1 (default, for 'RW1') or 2 ('RW2'), which specifies the prior distribution considered on the coefficients of the temporal B-splines. |
| ```prior.interaction``` | prior distributions for the espatio-temporal random effects,  one of either 1 (default, for 'Type I'), 2 ('Type II'), 3 ('Type III') or 4 ('Type IV'). |
| ```temp.corre``` | logical value (default TRUE); if FALSE, a models without correlation between the coefficients of the temporal P-splines is fitted. |
| ```order.B```	| order B-splines (default 3). |
| ```k.long```	| Number of internal intervals (x1: longitude), if NULL k.long = min(length(unique(x1))/4, 40). |
| ```k.lati```	| Number of internal intervals (x2: latitude), if NULL k.long = min(length(unique(x2))/4, 40). |
| ```k.time```	| Number of internal intervals (x3: time), if NULL k.long = min(length(unique(x3))/4, 40). |
| ```centered``` | logical value (default TRUE); TRUE to center the smooth functions, FALSE to center the coefficients. |
| ```strategy``` | one of either 'gaussian', 'simplified.laplace' (default), 'laplace' or 'adaptive', which specifies the approximation strategy considered |
| ```initial```	| logical value (default TRUE); if TRUE, the initial values ​​of the hyperparameters are estimates using standardized incidence rates; if FALSE, the initial values are 1. |
| ```rerun``` | logical value (default FALSE); if TRUE, this function will take the result in object, and rerun 'inla' again. |



# Function ```ps_st```

## Description
Fit a univariate Bayesian spatio-temporal P-spline models.

## Usage
```
ps_st(carto=NULL, data=NULL, ID.area=NULL, ID.year=NULL, crimes=NULL, Expcrimes=NULL,
      prior.spatial=1, prior.temporal=1, prior.interaction=2,
      order.B=3, k.long=NULL, k.lati=NULL, k.time=NULL, centered=TRUE,
      strategy="simplified.laplace", rerun=FALSE )
```

## Arguments

| Argument | Description |
|:---	| :--- | 
| ```carto```	| object of class ```SpatialPolygonsDataFrame``` or ```sf```. This object must contain at least the variable with the identifiers of the spatial areal units specified in  the argument ```ID.area```. | 
| ```data```	| object of class data.frame. This object must contain at least the	target variables of interest specified in the arguments ```ID.area```, ```ID.year```, ```crimes``` and ```Expcrimes```. |
| ```ID.area``` | character; name of the variable which contains the IDs of spatial areal units. |
| ```ID.year``` | character; name of the variable which contains the IDs of temporal units. |
| ```crimes``` | character; name of the variables that contain the observed number of crimes for each area and time unit. |
| ```Expcrimes```	| character; name of the variables that contain the expected number of crimes for each area and time unit. |
| ```prior.spatial``` | one of either 1 (default, for 'RW1') or 2 ('RW2'), which specifies the prior distribution considered on the coefficients of the spatial B-splines. |
| ```prior.temporal``` | one of either 1 (default, for 'RW1') or 2 ('RW2'), which specifies the prior distribution considered on the coefficients of the temporal B-splines. |
| ```prior.interaction``` | prior distributions for the espatio-temporal random effects,  one of either 1 (default, for 'Type I'), 2 ('Type II'), 3 ('Type III') or 4 ('Type IV'). |
| ```order.B```	| order B-splines (default 3). |
| ```k.long``` | Number of internal intervals (x1: longitude), if NULL k.long = min(length(unique(x1))/4, 40). |
| ```k.lati``` | Number of internal intervals (x2: latitude), if NULL k.long = min(length(unique(x2))/4, 40). |
| ```k.time``` | Number of internal intervals (x3: time), if NULL k.long = min(length(unique(x3))/4, 40). |
| ```centered``` | logical value (default TRUE); TRUE to center the smooth functions, FALSE to center the coefficients. |
| ```strategy``` | one of either 'gaussian', 'simplified.laplace' (default), 'laplace' or 'adaptive', which specifies the approximation strategy considered |
| ```rerun``` | logical value (default FALSE); if TRUE, this function will take the result in object, and rerun 'inla' again. |



# Functions ```generic_mps_spat``` and ```generic_mps_temp```

## Description
Implementation of spatial (```generic_mps_spat```) and temporal (```generic_mps_temp```generic_mps_temp')  multivariate P-spline models using the ```rgeneric``` model of INLA.

## Usage
```
generic_mps_spat(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                         "log.prior", "quit"), 
                 theta = NULL)

generic_mps_temp(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                         "log.prior", "quit"),
                 theta = NULL)                         
```
## Arguments

| Argument | Description |
|:---	| :--- |
| cmd	| Internall functions used by the ```rgeneric``` model to define the latent effect. | 
| theta	| Vector of hyperparameters. |