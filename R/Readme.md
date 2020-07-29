# R code

This folder contains the necessary R functions to fit the multivariate spatio-temporal areal models described in Vicente et al. (2020), and to reproduce the results.

The [data_Psplines.RData](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/data_Psplines.RData) file contains the following R objects:

- ```data```: contains the data set used. It is a dataframe with the following variables,
	- **dist**: Districts
	- **state**: State (Maharashtra)
	- **year**: Year (2001:2013)
	- **rape**: Observed number of Rape
	- **assault**: Observed number of Assault or criminal force to woman with intent to outrage her modesty
	- **cruelty**: Observed number of Cruelty by husband or relatives of husband
	- **kidnapping**: Observed number of Kidnapping and abduction of women
	- **e_rape**: Expected number of Rape
	- **e_assault**: Expected number of Assault or criminal force to woman with intent to outrage her modesty
	- **e_cruelty**: Expected number of Cruelty by husband or relatives of husband
	- **e_kidnapping**: Expected number of Kidnapping and abduction of women
	- **ID_area**: Area Identifiers (Districts)
	- **ID_year**: Time Identifiers (Year)


- ```carto```: SpatialPolygonDataFrame object with the cartography of the 34 districts of Maharashtra



The [run_multi_Centered_Psplines_models.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_multi_Centered_Psplines_models.R) file allow you to adjust the Multivariate P-spline models using INLA.

The [run_multi_Centered_Psplines_models_without_temp_corr.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_multi_Centered_Psplines_models_without_temp_corr.R) file allow you to adjust the Multivariate P-spline models (without temporal correlations) using INLA.


The [run_multi_Psplines_models.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_multi_Psplines_models.R) file allow you to adjust the Multivariate P-spline models with centered coefficients using INLA.

The [run_multi_Psplines_models_without_temp_corr.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_multi_Psplines_models_without_temp_corr.R) file allow you to adjust the Multivariate P-spline models (without temporal correlations) with centered coefficients using INLA.


The [run_Mmodels.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_Mmodels.R) file allow you to adjust the fixed-effect and random-effects M-models with a iCAR model for the spatial and temporal random effect, using INLA.

The [run_multi_Mmodels_Centered_Psplines_models.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_multi_Mmodels_Centered_Psplines_models.R) file allow you to adjust the fixed-effect/random-effects M-models with a iCAR model for the spatial random effect and P-spline models for the temporal random effect, using INLA.





The [run_univ_Centered_Psplines_models.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_univ_Centered_Psplines_models.R) file allow you to adjust the univariate P-spline models using INLA.

The [run_univ_Psplines_models.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_univ_Psplines_models.R) file allow you to adjust the univariate P-spline models with centered coefficients using INLA.

The [run_univ_iCAR.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_univ_iCAR.R) file allow you to adjust the univariate iCAR models using INLA.




The [functions](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/functions) folder contains the necessary functions to fit Multivariate P-spline models and M-models using INLA.

The [fbplot](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/fbplot) folder contains a modification of the 'fbplot' function of the 'fda' library.



The file [reproduce_paper_Psplines.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/reproduce_paper_Psplines.R) permit to reproduce the results given in the paper.



 
