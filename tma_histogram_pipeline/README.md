# README for Histogram Plot Script

## Overview

This R script is designed to generate two histogram plots based on user-defined parameters. It accepts a CSV file input containing data related to Nucleus and Cytoplasm DAB OD mean values, processes the data using various statistical methods, and generates histograms with Weibull distribution fitting. Additionally, the script can filter data based on further criteria (such as Nucleus Max Caliper and Cell DAB OD mean) to produce another histogram with specified cutoffs.

## Requirements

The script requires the following R libraries:
- `data.table`
- `tidyverse`
- `fitdistrplus`

Ensure the required packages are installed before running the script.



## Commands:

Rscript -e 'install.packages("tidyverse", repos="http://cran.rstudio.com/")'

Rscript -e 'library(tidyverse)'

Rscript -e 'install.packages("data.table", repos="http://cran.rstudio.com/")'

Rscript -e 'library(data.table)'

Rscript -e' install.packages("fitdistrplus",repos="http://cran.rstudio.com/")'

Rscript -e 'library(fitdistrplus)'


## Usage

Place the CSV file(s) you want to analyze in the same working directory as the script.

Run the script using the following command:

Rscript tma_histogram_pipeline.R  `input_file` `bin_size` `col1` `col1_cutoff` `col2` `col2_cutoff` `col3` `col3_cutoff` `col4` `col4_cutoff`

`input_file`: The name of the CSV file (e.g., "14895_TMA1.2401b_Breast_core_G-1_S9.6.csv").
`bin_size`: Number of bins for the histograms.
`col1`: Column for "Nucleus: DAB OD mean".
`col1_cutoff`: write 0 to delete all negative values. 
`col2`: Column for "Cytoplasm: DAB OD mean".
`col2_cutoff`: write 0 to delete all negative values
`col3`: Column for "Nucleus: Max caliper" (optional).
`col3_cutoff`: Cutoff value for col3 (e.g., 25). #all values below cutoff will be taken for further analysis.
`col4`: Column for "Cell: DAB OD mean" (optional).
`col4_cutoff`: Cutoff value for col4 (e.g., 0.18) # all values above this cutoff will be taken for further analysis.


Example code:

Rscript tma_histogram_pipeline.R 14895_TMA1.2401b_Breast_core_G-1_S9.6.csv 75 "Nucleus: DAB OD mean" 0 "Cytoplasm: DAB OD mean" 0

The input file is 14895_TMA1.2401b_Breast_core_G-1_S9.6.csv.
The bin size is set to 75.
The script first filters rows where "Nucleus: DAB OD mean" and "Cytoplasm: DAB OD mean" are greater than 0. Meaning all positive values in a given column is retained.

Or you can apply further filtering with "Nucleus: Max caliper ` 25" and "Cell: DAB OD mean ` 0.18".

Rscript tma_histogram_pipeline.R 14895_TMA1.2401b_Breast_core_G-1_S9.6.csv 75 "Nucleus: DAB OD mean" 0 "Cytoplasm: DAB OD mean" 0 "Nucleus: Max caliper" 25 "Cell: DAB OD mean" 0.18

## Script Logic
**1.Reading the Input Data: **
The script reads the CSV file and selects specific columns based on the provided column names.
Print all the parameters entered by users.
Filters rows where both "Nucleus: DAB OD mean" and "Cytoplasm: DAB OD mean" have positive values. Save both positive and negative values dataset separately in working directory.

**2.Perform and print statistical Calculations:**
-Performs a Kolmogorov-Smirnov (KS) test.
-Calculates ANOVA test
-Calculates the Euclidean distance between these values.

**3.First part of analysis:**
Create histogram with probability density function for both absolute and normalized frequency. Save all graph input and final graphs in the working directory.
Next, create histogram and fits Weibull distributions for both absolute and normalized frequency.
Save all graph input and final graphs in the working directory

**4.Second part of analysis (Optional):**
If cutoffs for "Nucleus: Max caliper" and/or "Cell: DAB OD mean" are provided, the script further filters the data based on these additional criteria. Steps 1,2 and 3 repeated again with filtered dataset. Generates another histogram with Weibull distribution fitting for the filtered data.
.

## Output files details:

**1. CSV files: **
a)	TMA with positive cells.csv: Contains rows where "Nucleus: DAB OD mean" and "Cytoplasm: DAB OD mean" are both greater than 0.

b)	TMA with negative cells.csv: Contains rows where these conditions are not met.

c)	graphinput no cutoff`filename`.csv: Data used to generate the first histogram.

d)	TMA matched with user-specified parameter cutoffs.csv: Data filtered by additional cutoff conditions.

e)	graphinput after cutoff`filename`.csv: Data used to generate the second histogram.

f) All statiscal tables basic stats and for Weibull parameters.


**2. Histogram Plots: **
Four histogram plots will be saved for each part of the histogram analysis. For the first part of the analysis (without additional cutoffs), the following names are used:

 - probability density graph no cutoff `filename`.tiff
 
 - probability density graph normalized no cutoff`filename`.tiff
 
 - 	Weibull graph no cutoff`filename`.tiff
 
 - 	Weibull graph normalized no cutoff`filename`.tiff
 

If additional cutoff criteria are provided, the second part of analysis will execute with following names:

 - probability density graph after cutoff `filename`.tiff
 
 - probability density graph normalized after cutoff`filename`.tiff
 
 - Weibull graph after cutoff`filename`.tiff
 
 - Weibull graph normalized after cutoff`filename`.tiff



## Notes
Ensure the CSV file contains the columns required for analysis.
The script generates a log in the terminal, detailing the number of positive and negative cells, statistical tests results, and confirmation of saved files.


## Contact

For any questions or issues, contact Jyoti Devendra Adala.


