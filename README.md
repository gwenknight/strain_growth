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
Then 2_analysis_cut.R for those that have a shoulder

Returns two files 
output/cut_all_ddm.csv: time series + whether odd or not + any "cut" analysis due to a shoulder with new cut exponential growth     
output/cut_all_param.csv: just summary parameters for each dataset

### "3_"
3_a_clean_label.R cleans / labels datasets or replicates to be removed: outputs param_labelled.csv
