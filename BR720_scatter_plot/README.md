#README: Comparative analysis for S9.6 and RNase III+H treatment effects

## Overview

This project analyzes the effects of RNase III and S9.6 treatment on cellular data derived from BR720 datasets. The primary objectives include:

1. Data preprocessing and cleaning.
2. Subsetting by specific cores.
3. Generating comparative statistical analyses and visualizations.
4. Exploring correlations using scatter plots and regression models.
5. Visualizing density distributions with 3D surface plots.


## Requirements

The script requires the following R libraries:
 - `data.table`
 - `tidyverse`
 - `plotly`
 - `MASS`

Ensure the required packages are installed before running the script.

## Input Data

-BR720_Rnase treated_IDC only_OSK_MERGED_new.csv
-BR720_S9.6 only_IDC only_OSK_MERGED_new.csv
Place these files in the working directory.

## Steps Performed

1. **Data Loading and Preprocessing:**
   - Load datasets using fread from data.table.-
   -Select columns of interest: Image, Name, Class, TMA core, Centroid X µm, Centroid Y µm, Hematoxylin: Nucleus: Mean, DAB: Nucleus: Mean.
   -Filter rows where Hematoxylin: Nucleus: Mean and DAB: Nucleus: Mean are greater than zero.

2. **Subsetting by Cores:**
   - Define specific cores (e.g., "C-8", "H-9").
   - Subset data for each core.
   - Save processed data for all classes and positive-class-only rows into separate .csv files.

3. **Statistical Analysis and Visualization:**
   - For each core, generate scatter plots with regression lines.
   - Perform correlation tests (Pearson) and fit linear models for:All class data, Positive class-only data, Comparison of RNase III and S9.6 treatments.

4. **3D Surface Plots:**
   - Visualize density distributions using plotly for each core and treatment group.


## Outputs

1. **Processed Datasets:**
-BR720_RNase_treated_all_class_[CORE].csv
-BR720_RNase_treated_no_negative_class_[CORE].csv
-BR720_S9.6_only_all_class_[CORE].csv
-BR720_S9.6_only_no_negative_class_[CORE].csv

2. **Graph:**
- BR720_RNase_treated_all_class_[CORE]_scatter_plot.tiff
-BR720_S9.6_only_all_class_[CORE]_scatter_plot.tiff

## How to Run

1. Place the input dataset in the working directory.
2. Run the script using R.
3. The outputs will be saved in the working directory.

## Notes

-Adjust the scale_x_log10 and scale_y_log10 parameters in plots if log scaling is not required.
-Modify the list of cores in cores to analyze different cores.

## Contact

For any questions or issues, contact Jyoti Devendra Adala.
