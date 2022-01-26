# strain_growth
Supporting code for "Dehydration tolerance in epidemic versus nonepidemic MRSA demonstrated by isothermal microcalorimetry", authored by: Valérie O. Baede a, Mehri Tavakol a, Margreet C. Vos a, Gwenan M. Knight b+, Willem J.B. van Wamel a+# on behalf of the MACOTRA study group. 

## Aim of code 
To extract key parameters from time series of bacterial growth 

## Aim of related paper 
To determine the impact of dessication on bacterial survival. 

To do this we used this code to extract key predictive parameters (e.g. time to maximal exponential growth) and then calculate a linear relationship between initial inoculum and time to peak to predict surviving bacterial levels. 

## FILE structure
- "data" holds the raw data for this project
- "output" holds the key parameter and plot outputs
- "proof_of_principle" holds the data and then the output for the comparison of optical density and calorimeter data
- "plots/figures" holds the paper figues 
- "plots/exp_growth" holds the exponential growth analysis 
- "plots/final_data_split_highlighted/" shows the final curves of all datasets used with the cut point highlighted (red point)
- "plots/linear_fit/" holds the summary output of linear model fit to the data 
- "plots/output_fit" holds the individual fit to each dataset

## CODE structure

### Code to data clean: "1_"
Run 1_data_cleaning.R 

This generates all the data sources for this paper by standardising the varying timeseries data in the data folder. Inspection of the growth curves for potential contamination leads to some data being removed. This can be seen in the plots generated in the individual R code for each data set (data cleaning 1). 

### Extract parameters: "2_": all strains

Run 2_analysis.R to extract all parameters needed for this analysis. 

This needs files 
- "grofit_functions.R" which takes the needed functions from the now unsupported GROFIT package. 
- "functions_for_heat_curves.R" which has the find peaks function and the main function for extracting the parameters. Here ODD is used to characterise any oddities in the curves. Most of this functionality is not used in this paper but is used to characterise differences between the curves (e.g. wide peaks). 

Returns two files 
- output/YOURNAMECODE_all_time_series_fit_params.csv": time series + whether odd or not + parameter values      
- output/YOURNAMECODE_all_model_fit_params.csv: summary parameters for each time series dataset  

### Extract parameters: "2_": non-MACOTRA strans
Run 2_analysis_non_mactora.R to get curves / analysis for the second figure of the paper based on the OD vs calorimeter analysis. This can only be run once the parameters have been extracted (in above "2_" code). 

### Proof of principle:
Run proof_of_principle.R to perform the comparison analysis of OD vs. calorimeter data using the same "non-macotra" strains above. 

Returns plots to "proof_of_principle" 

### Perform analysis: "3_": exponential growth cutoff determination 
In order to perform the analysis, the first step is to determine the threshold that is needed to determine how variable the exponential growth can be when the top 5% of variability is removed. This is determined by running 3_exponential_growth_variation.R. 

Returns 
- plots of individual strain exponential growth variation to plots/exp_growth.

### Perform analysis: "3_": 
Run 3_clean_exp_logred.R 

This needs files
- function_linear_model.R which calculates the linear model fit to the extracted data
