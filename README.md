# YangtzeRiverFlood
This Matlab code repository contains two folders: One is FloodAnalysis, which contains four main scripts for detecting, validating, and quantify flood events from lacustrine sedimentary data: 
(1) flood_analysis_main.m is the main script for detecting flood events from a lake sedimentary record and validate against the instrumental record; 
(2) h_threshold.m is used to optimize the threshold for identifying flood events; 
(3) reg_bootstrap.m is used to bootstrap the linear regression model for reconstructing flood magnitude; 
(4) flood_synthesis.m is used to generate a composite flood record from individual flood records. 

The other is ModEAnalysis, which contains four main scripts for analyzing and visualize the Modern Era Reanalysis (ModE-RA) products: 
(1) CONVxmonth.m is used to calculate the annual maximum consecutive 3-month mean vertically integrated moisture flux and convergence; 
(2) Rxmonth.m is used to calculate the annual maximum consecutive 3-month total precipitation; 
(3) WPSH_main.m is used to calculate the western Pacific subtropical high (WPSH) strength index through an EOF analysis of the annual maximum consecutive 3-month mean 850-hPa wind stream function; 
(4) ITCZ_main.m is used to calculate the intertropical convergence zone (ITCZ) latitudinal index from the annual maximum consecutive 3-month mean convergence of 850-hPa horizontal winds.
