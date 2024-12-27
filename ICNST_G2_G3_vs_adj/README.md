# README: S9.6 Signal Intensity Analysis in ICNST

## Overview

This project analyzes the S9.6 signal intensity across different grades in invasive carcinoma of no special type (ICNST). The script processes the dataset, normalizes it, and visualizes the results to compare nucleus signal intensity distributions across three categories: "Adjacent normal tissue," "Grade 2," and "Grade 3."

## Requirements

The script requires the following R libraries:
- `data.table`
- `tidyverse`

Ensure the required packages are installed before running the script.

## Input Data

The input dataset, `BC081120f_MERGED_SORTED_No Negative_OSK.csv`, must contain the following columns:
- `Nucleus: DAB OD mean`: Nucleus signal intensity values.
- `Pathology diagnosis`: Tissue pathology descriptions.
- `Grade`: Cancer grade information.
- `TMA core`: Tissue microarray core identifiers.

## Steps Performed

1. **Data Cleaning:**
   - Removes rows with negative nucleus signal intensity values.

2. **Filtering:**
   - Categorizes the data into:
     - Adjacent normal tissue types.
     - Grade 2 cancer tissue.
     - Grade 3 cancer tissue.

3. **Binning and Aggregation:**
   - Bins the nucleus signal intensity into 75 equal intervals.
   - Calculates mean signal intensity and frequency for each bin.

4. **Normalization:**
   - Computes normalized average frequencies across tissue microarray cores for each category.

5. **Visualization:**
   - Creates a density plot of signal intensities and overlay mean values by bin.

6. **Statistical Analysis:**
   - (Optional) Performs ANOVA to test for significant differences in mean signal intensities across categories.

## Outputs

1. **Processed Datasets:**
   - `BC08112_Fig7_S9.6_ICNST_grade2_vs_3_vs_adj_density_input.csv`: Intermediate dataset used for graph plotting.
   - `BC08112_Fig7_S9.6_ICNST_grade2_vs_3_vs_adj_geom_point.csv`: Dataset with binned data and normalized frequencies.

2. **Graph:**
   - `BC08112_Fig7_S9.6_ICNST_grade2_vs_3_vs_adj_core_wise_cell_avg_per_bin_normalised.tiff`: Density plot comparing signal intensities across categories.

## How to Run

1. Place the input dataset in the working directory.
2. Run the script using R.
3. The outputs will be saved in the working directory.

## Notes

- Update the input file name in `fread` if the dataset name changes.
- Ensure the dataset contains the required columns with accurate data types.
- Adjust the plot aesthetics as needed.

## Contact

For any questions or issues, contact Jyoti Devendra Adala.
