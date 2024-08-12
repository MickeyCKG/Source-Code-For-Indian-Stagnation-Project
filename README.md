# Source-Code-For-Indian-Stagnation-Project
Here are the MATLAB modules we developed to load and spatially interpolate the CMIP6 daily and monthly output. Please see our relevant research paper in Nature Communications (https://doi.org/10.1038/s41467-024-51462-y).

%% We included two Matlab modules (scripts) and three associated functions:

(1) Step1_Read_CMIP6_Daily.m (the module to read and regrid CMIP6 daily variables);

(2) Step1_Read_CMIP6_Monthly.m (the module to read and regrid CMIP6 monthly variables);

(3) Read_CMIP6_daily.m (the function to read CMIP6 daily variables, that is called in Step1_Read_CMIP6_Daily.m);

(4) Read_CMIP6_monthly.m (the function to read CMIP6 monthly variables, that is called in Step1_Read_CMIP6_Monthly.m);

(5) Regrid_Coef.m (the function to regrid CMIP6 variables with flux-conserving method, that is called both in Step1_Read_CMIP6_Daily.m and Step1_Read_CMIP6_Monthly.m);

%% Before use these scripts, please :

(1) Download CMIP6 data from https://aims2.llnl.gov/search/cmip6/;

(2) Using the 'GFDL_ESM4' model and the 'historical' scenario as an example, place the daily or monthly model output in './GFDL_ESM4/historical/'. 'GFDL_ESM4' and 'historical' should be accordingly assigned to the variables 'model_name' and 'scenario_name', respectively, in the Matlab modules.



