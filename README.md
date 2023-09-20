### order of analysis for method comparison ###

#for trapping data 
prelim_traps.R cleans raw data from Access database and produces file with all capture recapture records in CleanTrapData.csv and site location info in TrapSiteInfo.csv
#for count data
prlim_count.R cleans raw counts from original survey123 collected in 2021 and 2022 to be used by prelim_counts_updated.R
prelim_counts_updated.R takes cleaned 2021-2022 count data from prelim_counts.R and cleans new survey123 collected in 2023. Combines and cleans it. Produces clean counts formatted
as countsALL_longdf.csv with data not-summarized (by either method (visual or auditory) or type (transect or point), countsXmethodALL.csv with data summarized by method, 
and countALLinfo.csv showing the site info for sites with counts. 
#for predictor data
ClimPrep.R takes PRISM and weather station data for 30 years and summarizes it for temporal analysis with long-term trapping data. 
HabPrepRAP.R converts 30 x 30m rasters to 300 m stacks
HabCovExtractRAP.R matches habitat data to all sites surveyed. 
Data4JAnalysis.R takes habitat data extracted in HabCovExtractRAP.R for all sites surveyed during 2021-2023 and weather station data and links it to clean counts and clean trapping data.
Data4JAGs.R readies clean count and predictors for analysis in JAGs.
CountsBayesAnalysis.R takes data from Data4JAGs.R and runs an N-mixture model. 
CountsBayesPlotting.R takes model output and plots results.
ModelEval_CountsJAGs.R does model fit for counts analysis.
