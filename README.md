# strain_growth
Strain growth analysis of calorimeter data, pre and post dessication 


### "1_"
Code to data clean

## "2_"
Code to characterise strains 
- extract exponential growth
- time to peak 
- characterise if odd peak / width / shoulder

First 2_analysis.R for all strains

Returns two files 
output/cut_all_time_series_fit_params.csv": time series + whether odd or not + new cut exponential growth     
output/cut_all_model_fit_params.csv: summary parameters for each dataset

Returns plots
plot/shoulder_curves: just those strains with a shoulder highlights the point at which they are cut (first cut doesn't always give good time to peak)


### "3_"
3_a_clean_label.R cleans / labels datasets or replicates to be removed: outputs param_labelled.csv, and plots up to cut. Explores exponential growth



# Use exploring_exponential_growth.R to find cutoff

3_b_fit_linear.R checks exponential and fits linear
