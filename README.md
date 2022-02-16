# Fire-Prediction-Western-US
BurnArea_GAM_Obs_Forecasting_Ziter_validation_PDSImod_PCA_higherPeakSWE.R

Calculate best-fit and leave-one-year-out skill for 100-member ensemble using PDSI drought conditions with other antecedent climate conditions over a domain using the 150mm peak SWE threshold 

BurnArea_GAM_Obs_Forecasting_Ziter_validation_PDSImod_PCA.R

Calculate best-fit and leave-one-year-out skill for 100-member ensemble using PDSI drought conditions with other antecedent climate conditions over a domain using the 100mm peak SWE threshold 

BurnArea_GAM_Obs_Forecasting_Ziter_validation_SWEImod_PCA_higherPeakSWE.R

Calculate best-fit and leave-one-year-out skill for 100-member ensemble using snow drought conditions with other antecedent climate conditions over a domain using the 150mm peak SWE threshold 

BurnArea_GAM_Obs_Forecasting_Ziter_validation_SWEImod_PCA_lowerPeakSWE.R

Calculate best-fit and leave-one-year-out skill for 100-member ensemble using snow drought conditions with other antecedent climate conditions over a domain using the 50mm peak SWE threshold 

BurnArea_GAM_Obs_Forecasting_Ziter_validation_SWEImod_PCA_SummerVars.R

Calculate best-fit and leave-one-year-out skill for 100-member ensemble using PDSI drought conditions with other antecedent and in-year fire season climate conditions over a domain using the 100mm peak SWE threshold 

BurnArea_GAM_Obs_Forecasting_Ziter_validation_SWEImod_PCA.R

Calculate best-fit and leave-one-year-out skill for 100-member ensemble using snow drought conditions with other antecedent climate conditions over a domain using the 100mm peak SWE threshold 

Convert_MTBS_to_Grid.m

Bring MTBS fire shapes to MODIS burned area 500-m grid

Decrease_April1_SWE_NoahMP.m

Create state/restart files for Noah-MP simulation that impose a 30% increase and decrease to April 1 SWE from a reference Noah-MP simulation 

export_MCD64A1_timeseries.m

create monthly burned area files from MCD64A1 observed burned area

extract_modelout_1layer_daily_DecrSWE_cheyenne.m

extract daily time series of outputs from the Noah-MP simulation corresponding with a 30% decrease in April 1 SWE

extract_modelout_1layer_daily_IncrSWE_cheyenne.m

extract daily time series of outputs from the Noah-MP simulation corresponding with a 30% increase in April 1 SWE

extract_modelout_1layer_daily_reference_cheyenne.m

extract daily time series of outputs from the baseline/reference Noah-MP simulation 

Fire_Frequency_By_Month.m

Create boxplots plot burned area by each month, where boxplot variability is from variability across different years

Get_Covariates_and_BA_forecast_1984_2020_peakSWE_higher.m

Create annual time series of burned area and climate predictors from 1984-2020 for each elevation bin to be used in generalized additive model construction. Resulting outputs correspond with the domain using the 150 mm peak SWE threshold.

Get_Covariates_and_BA_forecast_1984_2020_peakSWE_lower.m

Create annual time series of burned area and climate predictors from 1984-2020 for each elevation bin to be used in generalized additive model construction. Resulting outputs correspond with the domain using the 50 mm peak SWE threshold.

Get_Covariates_and_BA_forecast_1984_2020.m

Create annual time series of burned area and climate predictors from 1984-2020 for each elevation bin to be used in generalized additive model construction. Resulting outputs correspond with the domain using the 100 mm peak SWE threshold.

Get_Monthly_PDSI_from_PRISM.m

Calculate monthly PDSI from PRISM precipitation and temperature

Get_UA_Annual_Peak_SWE.m

Get annual peak SWE from the UA SWE product

Get_UA_Total_Spring_SWE.m

Get cumulative SWE from the UA SWE product (used in SWEI calculation performed by Get_Covariates_and_BA_forecast_1984_2020*.m)

nearestneighbour.m

Matlab function to nearest neighbor match two spatial grids

NN_NLDAS_to_BAgrid.m

Nearest neighbor match NLDAS-2 grid to the MODIS 500 m grid

NN_UA_and_PRISM_to_BAgrid.m

Nearest neighbor match UA-SWE grid to the MODIS 500 m grid

Plot_Domain_for_Obs_Analysis.m

Plot burned fraction spatial distribution and scatter plot of domain burned area vs. total western US burned area used in the study domain figure.

Plot_NoahMP_April1SWE_Experiment.m

Plot time series from the Noah-MP experimentation analyzing how perturbations in April 1 SWE effect summer drought conditions

plot_Obs_forecast_BA_1984_2020_PDSImod_PCA_higherSWE.m

Plot time series and scatter plots comparing simulated and observed burned area for the region corresponding with a 150mm peak SWE threshold. Simulations are from GAMs that use PDSI drought conditions with other antecedent predictors.

plot_Obs_forecast_BA_1984_2020_PDSImod_PCA.m

Plot time series and scatter plots comparing simulated and observed burned area for the region corresponding with a 100mm peak SWE threshold. Simulations are from GAMs that use PDSI drought conditions with other antecedent predictors.

plot_Obs_forecast_BA_1984_2020_SWEImod_PCA_HigherPeakSWE.m

Plot time series and scatter plots comparing simulated and observed burned area for the region corresponding with a 150mm peak SWE threshold. Simulations are from GAMs that use snow drought conditions with other antecedent predictors.

plot_Obs_forecast_BA_1984_2020_SWEImod_PCA_LowerPeakSWE.m

Plot time series and scatter plots comparing simulated and observed burned area for the region corresponding with a 50mm peak SWE threshold. Simulations are from GAMs that use snow drought conditions with other antecedent predictors.

plot_Obs_forecast_BA_1984_2020_SWEImod_PCA_SummerVars.m

Plot time series and scatter plots comparing simulated and observed burned area for the region corresponding with a 100mm peak SWE threshold. Simulations are from GAMs that use snow drought conditions with other antecedent and in-year fire season predictors.

plot_Obs_forecast_BA_1984_2020_SWEImod_PCA.m

Plot time series and scatter plots comparing simulated and observed burned area for the region corresponding with a 100mm peak SWE threshold. Simulations are from GAMs that use snow drought conditions with other antecedent predictors.

Plot_Predictor_Sensitivities_PCA_ensemble_higherPeakSWE.m

Plot importance of antecedent predictor figures corresponding to the domain using the 150 mm peak SWE threshold. 

Plot_Predictor_Sensitivities_PCA_ensemble_lowerPeakSWE.m

Plot importance of antecedent predictor figures corresponding to the domain using the 50 mm peak SWE threshold. 

Plot_Predictor_Sensitivities_PCA_ensemble.m

Plot importance of antecedent predictor figures corresponding to the domain using the 100 mm peak SWE threshold. 

plot_scatter_of_covariates.m

plot scatter plots of predictors and burned area on an annual time scale

Report_Npredictor_Table.m

Plot Taylor skill score based on leave-one-year-out cross validation for 100-model ensemble for models using 2-11 predictors. Used to select optimal number of predictors.

Select_BestMod_Zbins_ObsAnalysis_PDSI_higherPeakSWE.R

Select best 100-member ensemble for GAMs using PDSI drought conditions with other antecedent predictors corresponding to the domain using the 150 mm peak SWE threshold. 

Select_BestMod_Zbins_ObsAnalysis_PDSI_PCA_skillTable.R

Calculate and export leave-one-year-out Taylor skill scores for 100-member GAM ensembles using 2-11 predicotrs. Used in Report_Npredictor_Table.m. Domain corresponds with the 100 mm peak SWE threshold. Models use PDSI drought conditions with other antecedent predictors.

Select_BestMod_Zbins_ObsAnalysis_PDSI_PCA.R

Select best 100-member ensemble for GAMs using PDSI drought conditions with other antecedent predictors corresponding to the domain using the 100 mm peak SWE threshold. 

Select_BestMod_Zbins_ObsAnalysis_SnowDrought_higherPeakSWE.R

Select best 100-member ensemble for GAMs using snow drought conditions with other antecedent predictors corresponding to the domain using the 150 mm peak SWE threshold. 

Select_BestMod_Zbins_ObsAnalysis_SnowDrought_lowerPeakSWE.R

Select best 100-member ensemble for GAMs using snow drought conditions with other antecedent predictors corresponding to the domain using the 50 mm peak SWE threshold. 

Select_BestMod_Zbins_ObsAnalysis_SnowDrought_PCA_skillTable.R

Calculate and export leave-one-year-out Taylor skill scores for 100-member GAM ensembles using 2-11 predicotrs. Used in Report_Npredictor_Table.m. Domain corresponds with the 100 mm peak SWE threshold. Models use snow drought conditions with other antecedent predictors.

Select_BestMod_Zbins_ObsAnalysis_SnowDrought_PCA_SummerVars.R

Select best 100-member ensemble for GAMs using snow drought conditions with other antecedent and year-of fire season predictors corresponding to the domain using the 100 mm peak SWE threshold. 


Select_BestMod_Zbins_ObsAnalysis_SnowDrought_PCA.R

Select best 100-member ensemble for GAMs using snow drought conditions with other antecedent predictors corresponding to the domain using the 100 mm peak SWE threshold. 

SnowDrought_BestMod_Zbins_ObsAnalysis_PredictorSensitivity_PCA_higherPeakSWE.R

Perform leave-one-column-out (LOCO) analysis for GAMs using snow drought conditions with other antecedent predictors over the domain using the 150 mm peak SWE threshold. 

SnowDrought_BestMod_Zbins_ObsAnalysis_PredictorSensitivity_PCA_lowerPeakSWE.R

Perform leave-one-column-out (LOCO) analysis for GAMs using snow drought conditions with other antecedent predictors over the domain using the 50 mm peak SWE threshold. 

SnowDrought_BestMod_Zbins_ObsAnalysis_PredictorSensitivity_PCA.R

Perform leave-one-column-out (LOCO) analysis for GAMs using snow drought conditions with other antecedent predictors over the domain using the 100 mm peak SWE threshold. 


