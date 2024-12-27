# README:Analysis of S9.6 and RNase III+H Signal Intensities Across Stages in ICNST

## Overview

This script analyzes the distribution of signal intensities (DAB: Nucleus: Mean) for two experimental conditions—S9.6 and RNase III+H—across various Stages (Normal, Stage 1, Stage 2, and Stage 3) in ICNST. The outputs include density plots, statistical summaries, and files for further analysis.

## Requirements

The script requires the following R libraries:
- `data.table`
- `tidyverse`

Ensure the required packages are installed before running the script.

## Input Data

-BR720_IDC_only_RNase_treated_osk_No_negative.csv: Data for RNase III+H-treated samples (Stage 1, 2, 3).
-BR720_IDC_only_S9.6_positive_osk_No_negative.csv: Data for S9.6-treated samples (Stage 1, 2, 3).
-BR720_RNase_treated_N_No_Negative.csv: Data for RNase III+H-treated normal samples.
-BR720_S9.6_N_No_negative.csv: Data for S9.6-treated normal samples.

## Steps Performed

1. **Data Cleaning:**
   - Removes rows with negative nucleus signal intensity values.

2. **Filtering:**
   - Filter datasets into separate Stage-specific subsets (Stage 1, Stage 2, Stage 3, and Normal).
   - Add identifiers (e.g., "Stage 1") to facilitate group comparisons

3. **Data Transformation:**
   -Combine Datasets: Merge subsets into one dataset for combined analysis.
   -Binning: Determine the range of DAB: Nucleus: Mean.Bin values into 75 equal intervals and calculate midpoints for plotting.

4. **Normalization:**
   - Computes normalized average frequencies across tissue microarray cores for each category.

5. **Visualization:**
   - Creates a density plot of signal intensities and overlay mean values by bin.

6. **Statistical Analysis:**
   - Performs ANOVA to test for significant differences in DAB Nucleus mean signal intensities across categories.
   -Tukey HSD:Perform pairwise comparisons to identify specific group differences

## Outputs

1. **Processed Datasets:**
  -BR720_fig8a_Stage1_vs_2_vs_3_vs_normal_RNase_density_input.csv: RNase III+H dataset after preprocessing.

  -BR720_fig8a_Stage1_vs_2_vs_3_vs_normal_RNase_geom_point.csv: RNase III+H dataset ready for plotting.
  -BR720_fig8a_Stage1_vs_2_vs_3_vs_normal_S9.6_density_input.csv: S9.6 dataset after preprocessing.
  -BR720_fig8a_Stage1_vs_2_vs_3_vs_normal_S9.6_geom_point.csv: S9.6 dataset ready for plotting.

2. **Graph:**
  - BR720_fig8a_ICNST_Stage1_vs_2_vs_3_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff: RNase III+H density plot (bandwidth=0.05).
-BR720_fig8a_ICNST_Stage1_vs_2_vs_3_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff: RNase III+H density plot (bandwidth=0.02).
- BR720_fig8a_ICNST_Stage1_vs_2_vs_3_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff: RNase III+H density plot (bandwidth=0.05).
-BR720_fig8a_ICNST_Stage1_vs_2_vs_3_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff: RNase III+H density plot (bandwidth=0.02).

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


