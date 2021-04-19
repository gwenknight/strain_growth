# strain_growth
Supporting code for Baede et al. (2021)

## Aim of code 
To extract key parameters from time series of bacterial growth 

## Aim of related paper 
To determine the impact of dessication on bacterial survival. 

To do this we used this code to extract key predictive parameters (e.g. time to maximal exponential growth) and then calculate a linear relationship between initial inoculum and time to peak to predict surviving bacterial levels. 

## FILE structure
- "data" holds the raw data for this project
- "output" holds the key parameter and plot outputs

## CODE structure

### Code to data clean: "1_"
Run 1_data_cleaning.R 

This generates all the data sources for this paper by standardising the varying timeseries data in the data folder. Inspection of the growth curves for potential contamination leads to some data being removed. This can be seen in the plots generated in the individual R code for each data set (data cleaning 1). 

## Extract parameters: "2_"
Run 2_analysis.R

This needs files 
"grofit_functions.R" which takes the needed functions from the now unsupported GROFIT package. 
"functions_for_heat_curves.R" which has the find peaks function and the main function for extracting the parameters. Here ODD is used to characterise any oddities in the curves. Most of this functionality is not used in this paper but is used to characterise differences between the curves (e.g. wide peaks). 

Returns two files 
output/YOURNAMECODE_all_time_series_fit_params.csv": time series + whether odd or not + parameter values      
output/YOURNAMECODE_all_model_fit_params.csv: summary parameters for each time series dataset  

Returns plots
plot/shoulder_curves: just those strains with a shoulder highlights the point at which they are cut (first cut doesn't always give good time to peak)


### "3_"
3_a_clean_label.R cleans / labels datasets or replicates to be removed: outputs param_labelled.csv, and plots up to cut. Explores exponential growth



# Use exploring_exponential_growth.R to find cutoff

3_b_fit_linear.R checks exponential and fits linear
