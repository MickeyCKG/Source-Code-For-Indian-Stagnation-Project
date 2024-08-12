# Source-Code-For-Indian-Stagnation-Project
Here are the MATLAB modules we developed to load and spatially interpolate the CMIP6 daily and monthly output. Please see our relevant research paper in Nature Communications (https://doi.org/10.1038/s41467-024-51462-y).

We included two Matlab modules (scripts) and three associated functions:

(1) Step1_Read_CMIP6_Daily.m (the module to read and regrid CMIP6 daily variables);

(2) Step1_Read_CMIP6_Monthly.m (the module to read and regrid CMIP6 monthly variables);

(3) Read_CMIP6_daily.m (the function to read CMIP6 daily variables, that is called in Step1_Read_CMIP6_Daily.m);

(4) Read_CMIP6_monthly.m (the function to read CMIP6 monthly variables, that is called in Step1_Read_CMIP6_Monthly.m);

(5) Regrid_Coef.m (the function to regrid CMIP6 variables with flux-conserving method, that is called both in Step1_Read_CMIP6_Daily.m and Step1_Read_CMIP6_Monthly.m);
