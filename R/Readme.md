# R code

This folder contains the necessary R functions to fit the multivariate spatio-temporal areal models described in Vicente et al. (2021), and to reproduce the results.

The ```data_Psplines.RData``` file contains the following R objects:

- ```data```: contains the data set used. It is a dataframe with the following variables,
	- **dist**: Districts
	- **state**: States (Maharashtra)
	- **year**: Years (2001:2013)
	- **rape**: Observed number of Rapes
	- **assault**: Observed number of Assaults or criminal force to woman with intent to outrage her modesty
	- **cruelty**: Observed number of cases of the crime "Cruelty by husband or relatives of husband"
	- **kidnapping**: Observed number of cases of the crime "Kidnapping and abduction of women"
	- **e_rape**: Expected number of Rapes
	- **e_assault**: Expected number of Assaults or criminal force to woman with intent to outrage her modesty
	- **e_cruelty**: Expected number of of cases of the crime "Cruelty by husband or relatives of husband"
	- **e_kidnapping**: Expected number of of cases of the crime "Kidnapping and abduction of women"
	- **ID_area**: Area Identifiers (Districts)
	- **ID_year**: Time Identifiers (Year)


- ```carto```: SpatialPolygonDataFrame object with the cartography of the 34 districts of Maharashtra


The [run_mps_st.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_mps_st.R) file allows you to fit multivariate Bayesian spatio-temporal P-spline models using INLA.

The [run_ps_st.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_ps_st.R) file allows you to fit univariate Bayesian spatio-temporal P-spline models using INLA.


The [functions](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/functions) folder contains the necessary functions to fit univariate and multivariate P-spline models using INLA.

The [fbplot](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/fbplot) folder contains a modification of the 'fbplot' function of the 'fda' library.


The file [reproduce_results_paper.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/reproduce_results_paper.R) permits to reproduce the results given in the paper.

**IMPORTANT** 
You need to run [run_mps_st.R](https://github.com/spatialstatisticsupna/Multivariate_spatio_temporal_P_spline/blob/master/R/run_mps_st.R) with the 'simplified.laplace' approximation to reproduce the results of the paper.
