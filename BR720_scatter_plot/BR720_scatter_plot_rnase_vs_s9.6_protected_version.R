# ------------------------------------------------------------------------------
# This code is part of paper: RNA:DNA hybrids determine tissue pathological heterogeneity and aggressiveness in clinically diagnosed cancers. 
# Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
# For updates and contributions, visit : https://github.com/adalaj

# Purpose: This R code is designed to analyze and visualize the effects of RNase treatment and S9.6 treatment on a given dataset. 
#It processes data, performs filtering, generates scatter plots, computes correlations, fits linear models, and creates 3D surface plots for selected treatment groups and core samples
# ------------------------------------------------------------------------------


#JA: load required libraries
library(data.table)
library(tidyverse)
library(plotly)
library(MASS)


#JA: Step1: Load dataset and preprocessing
rnase_tted<- fread("BR720_Rnase treated_IDC only_OSK_MERGED_new.csv")

#JA:select desired columns
rnase_tted <- rnase_tted %>% dplyr::select(Image, Name, Class, `TMA core`, `Centroid X µm`, `Centroid Y µm`, `Hematoxylin: Nucleus: Mean`, `DAB: Nucleus: Mean`)
#JA:685694 in rnase treated


#JA:filtering rows where both Hematoxylin and DAB nucleus greater than zero
filtered_rnase_tted<- rnase_tted %>% filter(`Hematoxylin: Nucleus: Mean` >0, #JA:previously i missed this step before giving data to ojashpi 
                                            `DAB: Nucleus: Mean`>0)
#JA:685624  


  
not_tted<- fread("BR720_S9.6 only_IDC only_OSK_MERGED_new.csv")
not_tted <- not_tted %>% dplyr::select(Image, Name, Class, `TMA core`, `Centroid X µm`, `Centroid Y µm`, `Hematoxylin: Nucleus: Mean`, `DAB: Nucleus: Mean`)
#JA:825712 in not treated that has only S9.6


#JA:filtering rows where both Hematoxylin and DAB nucleus greater than zero
filtered_not_tted<- not_tted %>% filter(`Hematoxylin: Nucleus: Mean` >0, 
                                            `DAB: Nucleus: Mean`>0)
#JA:825644  


#JA: Step2: Subsetting by Cores

cores<- list("C-8", "H-9")
for ( i in cores) {
  core_tted<- filtered_rnase_tted %>% filter(`TMA core` == i)
  fwrite(core_tted, paste("BR720_RNase_treated_all_class_", i, ".csv", sep = ""))
  
  #JA:removes all class negatives and keep only class positives
  core_tted_positive_class_only<- core_tted %>% filter(`Class` !="Negative")
  print(paste("Total number of positive class cells in RNase treated(BR720) =", nrow(core_tted_positive_class_only)))
  fwrite(core_tted_positive_class_only, paste("BR720_RNase_treated_no_negative_class_", i, ".csv", sep = ""))
  
  #JA:repeat this code for no treatment dataset as well
  core_not_tted<- filtered_not_tted %>% filter(`TMA core` == i)
  fwrite(core_not_tted, paste("BR720_S9.6_only_all_class_", i, ".csv", sep = ""))
 
  #JA:removes all class negatives and keep only class positives
  core_not_tted_positive_class_only<- core_not_tted %>% filter(`Class` !="Negative")
  print(paste("Total number of positive class cells in S9.6 only(BR720) =", nrow(core_not_tted_positive_class_only)))
  fwrite(core_not_tted_positive_class_only, paste("BR720_S9.6_only_no_negative_class_", i, ".csv", sep = ""))
  }






#JA:Step3: Visualization, statiscal analysis and comparative analysis


#JA:1) RNAase treatment all class
#JA:2) RNase only treated no negative class
#JA:3) S9.6 treatment all class
#JA:4) S9.6 only no negative class

cores<- list ("C-8", "H-9")
for (i in cores){
  #JA:prepare filename for the current core
  print(i)
  filename1<- paste("BR720_RNase_treated_all_class_", i, ".csv", sep = "")
  rnase_tted<- fread(filename1, sep = ",", header = TRUE)
  print(filename1)
  
  #JA:generate plot
  plot1<- ggplot(rnase_tted, aes(x = `Hematoxylin: Nucleus: Mean`, y = `DAB: Nucleus: Mean`))+
    geom_point(color= "red") +                   # Scatter plot
    geom_smooth(method = "lm", se = FALSE, color= "black") +  # Regression line
    labs(
      title = paste("Core:",i, "RNase III + H  treatment (All Class)"), 
      x = "Hematoxylin: Nucleus: Mean (log10)", 
      y = "DAB:Nucleus: Mean (log10)")+
    
    scale_x_log10(limits = c(0.05, 1),breaks= scales::log_breaks(base = 10, n= 10))+ #remove this if you dont want log scale
    scale_y_log10(limits = c(0.01,1.00),breaks= c(0.01, 0.10, 1.00))+ #remove this if you dont want log scale
    theme(text = element_text(size = 50), 
          panel.grid = element_blank(),
          axis.line = element_line(color= "black"),
          panel.background = element_blank())
  
  ggsave(paste("BR720_RNase_treated_all_class_",i,"_scatter_plot.tiff", sep = ""),
         plot= plot1, width = 20, height = 10, dpi = 600)
  
  
  
  #JA:Print details to the log
  
  print(paste("Total number of cells in RNase treated(BR720)=", nrow(rnase_tted)))
  
  #JA:perform coorelation test
  correlation_test_filename1<- cor.test(
    rnase_tted$`Hematoxylin: Nucleus: Mean`, 
    rnase_tted$`DAB: Nucleus: Mean`, 
    method = "pearson")
  
  print(correlation_test_filename1)
  
  # JA:Fit linear model and print summary
  
  lm_model1 <- lm(`DAB: Nucleus: Mean` ~ `Hematoxylin: Nucleus: Mean`, data = rnase_tted)
  print(summary(lm_model1))
  
  print(paste("Transforming the variables to log scale for ", filename1, sep = ""))
  rnase_tted$log_hematoxylin_nucleus_mean <- log10(rnase_tted$`Hematoxylin: Nucleus: Mean`)
  rnase_tted$log_DAB_nucleus_mean <- log10(rnase_tted$`DAB: Nucleus: Mean`)
  
  #JA:perform coorelation test for log transformed 
  log_transformed_correlation_test_filename1<- cor.test(
    rnase_tted$log_hematoxylin_nucleus_mean, 
    rnase_tted$log_DAB_nucleus_mean, 
    method = "pearson")
  
  print(log_transformed_correlation_test_filename1)
  log_transformed_lm_model1 <- lm(log_DAB_nucleus_mean ~ log_hematoxylin_nucleus_mean, data = rnase_tted)
  print(summary(log_transformed_lm_model1))
  
  fwrite(rnase_tted, filename1)
  
  
  print("Make similar graph for rnase treated in no negative class input, plot2")
  
  filename2<- paste("BR720_RNase_treated_no_negative_class_", i, ".csv", sep = "")
  rnase_tted_pos<- fread(filename2, sep = ",", header = TRUE)
  print(filename2)
  
  
  #JA:generate plot
  plot2<- ggplot(rnase_tted_pos, aes(x = `Hematoxylin: Nucleus: Mean`, y = `DAB: Nucleus: Mean`))+
    geom_point(color= "red") +                   # Scatter plot
    geom_smooth(method = "lm", se = FALSE, color= "black") +  # Regression line
    labs(
      title = paste("Core:",i, "RNaseIII + H treatment(No Negative Class)"), 
      x = "Hematoxylin: Nucleus: Mean (log10)", 
      y = "DAB: Nucleus: Mean (log10)")+
    
    scale_x_log10(limits = c(0.05, 1),breaks= scales::log_breaks(base = 10, n= 10))+ #remove this if you dont want log scale
    scale_y_log10(limits = c(0.01,1.00),breaks= c(0.01, 0.10, 1.00))+ #remove this if you dont want log scale
    theme(text = element_text(size = 45), 
          panel.grid = element_blank(),
          axis.line = element_line(color= "black"),
          panel.background = element_blank())
  
  ggsave(paste("BR720_RNase_treated_no_negative_class_",i,"_scatter_plot.tiff", sep = ""),
         plot= plot2, width = 20, height = 10, dpi = 600)
  
  
  #JA:Print details to the log
  
  print(paste("Total number of positive class cells in RNase treated(BR720)=", nrow(rnase_tted_pos)))
  
  #JA:perform coorelation test
  correlation_test_filename2<- cor.test(
    rnase_tted_pos$`Hematoxylin: Nucleus: Mean`, 
    rnase_tted_pos$`DAB: Nucleus: Mean`, 
    method = "pearson")
  
  print( correlation_test_filename2)
  
  #JA: Fit linear model and print summary
  
  lm_model2 <- lm(`DAB: Nucleus: Mean` ~ `Hematoxylin: Nucleus: Mean`, data = rnase_tted_pos)
  print(summary(lm_model2))  
  
  
  
  print("Make similar graph for S9.6 all class input, plot3")
  
  filename3<- paste("BR720_S9.6_only_all_class_", i, ".csv", sep = "")
  not_tted<- fread(filename3, sep = ",", header = TRUE)
  print(filename3)
  
  
  #JA:generate plot
  plot3<- ggplot(not_tted, aes(x = `Hematoxylin: Nucleus: Mean`, y = `DAB: Nucleus: Mean`))+
    geom_point(color= "blue") +                   # Scatter plot
    geom_smooth(method = "lm", se = FALSE, color= "black") +  # Regression line
    labs(
      title = paste("Core:",i, "S9.6 only (All Class)"), 
      x = "Hematoxylin: Nucleus: Mean (log10)", 
      y = "DAB: Nucleus: Mean (log10)")+
    
    scale_x_log10(limits = c(0.05, 1),breaks= scales::log_breaks(base = 10, n= 10))+ #remove this if you dont want log scale
    scale_y_log10(limits = c(0.01,1.00),breaks= c(0.01, 0.10, 1.00))+ #remove this if you dont want log scale
    theme(text = element_text(size = 50), 
          panel.grid = element_blank(),
          axis.line = element_line(color= "black"),
          panel.background = element_blank())
  
  ggsave(paste("BR720_S9.6_only_all_class_",i,"_scatter_plot.tiff", sep = ""),
         plot= plot3, width = 20, height = 10, dpi = 600)
  
  
  #JA:Print details to the log
  
  print(paste("Total number of cells in S9.6 only(BR720)=", nrow(not_tted)))
  
  #JA:perform coorelation test
  correlation_test_filename3<- cor.test(
    not_tted$`Hematoxylin: Nucleus: Mean`, 
    not_tted$`DAB: Nucleus: Mean`, 
    method = "pearson")
  
  print( correlation_test_filename3)
  
  #JA:Fit linear model and print summary
  
  lm_model3 <- lm(`DAB: Nucleus: Mean` ~ `Hematoxylin: Nucleus: Mean`, data = not_tted)
  print(summary(lm_model3))
  
  print("Make similar graph for S9.6 negative class input, plot4")
  
  filename4<- paste("BR720_S9.6_only_no_negative_class_", i, ".csv", sep = "")
  not_tted_pos<- fread(filename4, sep = ",", header = TRUE)
  print(filename4)
  
  
  #JA:generate plot
  plot4<- ggplot(not_tted_pos, aes(x = `Hematoxylin: Nucleus: Mean`, y = `DAB: Nucleus: Mean`))+
    geom_point(color= "blue") +                   # Scatter plot
    geom_smooth(method = "lm", se = FALSE, color= "black") +  # Regression line
    labs(
      title = paste("Core:",i, "S9.6 only (No Negative Class)"), 
      x = "Hematoxylin: Nucleus: Mean (log10)", 
      y = "DAB: Nucleus: Mean (log10)")+
    
    scale_x_log10(limits = c(0.05, 1),breaks= scales::log_breaks(base = 10, n= 10))+ #remove this if you dont want log scale
    scale_y_log10(limits = c(0.01,1.00),breaks= c(0.01, 0.10, 1.00))+ #remove this if you dont want log scale
    theme(text = element_text(size = 50), 
          panel.grid = element_blank(),
          axis.line = element_line(color= "black"),
          panel.background = element_blank())
  
  ggsave(paste("BR720_S9.6_only_no_negative_class_",i,"_scatter_plot.tiff", sep = ""),
         plot= plot4, width = 20, height = 10, dpi = 600)
  
  
  #JA:Print details to the log
  
  print(paste("Total number of positive class cells in S9.6 only(BR720)=", nrow(not_tted_pos)))
  
  #JA:perform coorelation test
  correlation_test_filename4<- cor.test(
    not_tted_pos$`Hematoxylin: Nucleus: Mean`, 
    not_tted_pos$`DAB: Nucleus: Mean`, 
    method = "pearson")
  
  print( correlation_test_filename4)
  
  #JA:Fit linear model and print summary
  
  lm_model4 <- lm(`DAB: Nucleus: Mean` ~ `Hematoxylin: Nucleus: Mean`, data = not_tted_pos)
  print(summary(lm_model4))
  
}


#JA: Create 3D surface plot for all four plots above
cores<- list ("C-8", "H-9")

for (i in cores){
  #JA:prepare filename for the current core
  print(i)
  filename1<- paste("BR720_RNase_treated_all_class_", i, ".csv", sep = "")
  rnase_tted<- fread(filename1, sep = ",", header = TRUE)
  print(filename1)
  print("vs")
  
  filename3<- paste("BR720_S9.6_only_all_class_", i, ".csv", sep = "")
  not_tted<- fread(filename3, sep = ",", header = TRUE)
  print(filename3)
  
  #JA:make z axis
  z_rnase_tted <- with(rnase_tted, kde2d(`Hematoxylin: Nucleus: Mean`, `DAB: Nucleus: Mean`, n = 50))
  z_not_tted <- with(not_tted, kde2d(`Hematoxylin: Nucleus: Mean`, `DAB: Nucleus: Mean`, n = 50))
  
  
  #JA:generate plot
  fig1<- plot_ly()
  fig1 <- fig1 %>% add_surface(
    x = z_rnase_tted$x,
    y = z_rnase_tted$y,
    z = z_rnase_tted$z,
    opacity = 0.98,
    surfacecolor = z_rnase_tted$z,  #JA: This uses the 'z' values to color the surface
    colorscale = "Reds",    #JA: Blue color scale
    showscale = TRUE,        #JA: Show color scale on the side
    name = "RNase III + H treated"    #JA: Legend name for C8 treated
  )
  
  #JA: Add surface for C8 No Treatment, with red color scale and legend
  fig1 <- fig1 %>% add_surface(
    x = z_not_tted$x,
    y = z_not_tted$y,
    z = z_not_tted$z,
    opacity = 0.98,
    surfacecolor = z_not_tted$z,  #JA: This uses the 'z' values to color the surface
    colorscale = "Blues",    #JA: Red color scale
    showscale = TRUE,       #JA: Show color scale on the side
    name = "S9.6 only"  #JA: Legend name for C8 no treatment
  )
  #JA: Customize layout with title, subtitle, and font size for all text
  
  fig1 <- fig1 %>% layout(
    title = paste("Core:", i, " ",":3D Surface Plot for S9.6 only (blue) and RNase III+ H Treatment (red) in all cells", sep=""),
    scene = list(
      xaxis = list(
        title = "Hematoxylin: Nucleus: Mean",
        #JA:range = c(0, 0.3),  #JA: Limit x-axis range
        titlefont = list(size = 25),  #JA: Font size for x-axis title
        tickfont = list(size = 18)  #JA: Font size for x-axis tick labels
      ),
      yaxis = list(
        title = "DAB: Nucleus: Mean",
        titlefont = list(size = 18),  #JA: Font size for y-axis title
        tickfont = list(size = 18)  #JA: Font size for y-axis tick labels
      ),
      zaxis = list(
        title = "Frequency",
        titlefont = list(size = 18),  #JA: Font size for z-axis title
        tickfont = list(size = 18)  #JA: Font size for z-axis tick labels
      )
    ),
    legend = list(
      font = list(size = 20)  #JA: Font size for legend labels
    ),
    titlefont = list(size = 15)  #JA: Font size for the title
  )
  
  fig1
  htmlwidgets::saveWidget(fig1, paste("3D_surface_plot_fig5a_",i ,".html", sep = ""))
  
  #JA: Repeat the same with no negative class
  filename2<- paste("BR720_RNase_treated_no_negative_class_", i, ".csv", sep = "")
  rnase_tted_pos<- fread(filename2, sep = ",", header = TRUE)
  print(filename2)
  print("vs")
  
  filename4<- paste("BR720_S9.6_only_no_negative_class_", i, ".csv", sep = "")
  not_tted_pos<- fread(filename4, sep = ",", header = TRUE)
  print(filename4)
  
  #JA:#JA:make z axis
  z_rnase_tted_pos <- with( rnase_tted_pos, kde2d(`Hematoxylin: Nucleus: Mean`, `DAB: Nucleus: Mean`, n = 50))
  z_not_tted_pos <- with(not_tted_pos, kde2d(`Hematoxylin: Nucleus: Mean`, `DAB: Nucleus: Mean`, n = 50))
  
  
  #JA:generate plot
  fig2<- plot_ly()
  fig2 <- fig2 %>% add_surface(
    x = z_rnase_tted_pos$x,
    y = z_rnase_tted_pos$y,
    z = z_rnase_tted_pos$z,
    opacity = 0.98,
    surfacecolor = z_rnase_tted_pos$z,  #JA: This uses the 'z' values to color the surface
    colorscale = "Reds",    #JA: Blue color scale
    showscale = TRUE,        #JA: Show color scale on the side
    name = "RNase III + H treated"    #JA: Legend name for C8 treated
  )
  
  #JA: Add surface for C8 No Treatment, with red color scale and legend
  fig2 <- fig2 %>% add_surface(
    x = z_not_tted_pos$x,
    y = z_not_tted_pos$y,
    z = z_not_tted_pos$z,
    opacity = 0.98,
    surfacecolor = z_not_tted_pos$z,  #JA: This uses the 'z' values to color the surface
    colorscale = "Blues",    #JA: Red color scale
    showscale = TRUE,       #JA: Show color scale on the side
    name = "S9.6 only"  #JA: Legend name for C8 no treatment
  )
  #JA: Customize layout with title, subtitle, and font size for all text
  
  fig2 <- fig2 %>% layout(
    title = paste("Core:", i, " ",":3D Surface Plot for S9.6 only (blue) and RNase III+ H Treatment (red) in no negative cells", sep=""),
    scene = list(
      xaxis = list(
        title = "Hematoxylin: Nucleus: Mean",
        #JA:range = c(0, 0.3),  #JA: Limit x-axis range
        titlefont = list(size = 25),  #JA: Font size for x-axis title
        tickfont = list(size = 18)  #JA: Font size for x-axis tick labels
      ),
      yaxis = list(
        title = "DAB: Nucleus: Mean",
        titlefont = list(size = 18),  #JA: Font size for y-axis title
        tickfont = list(size = 18)  #JA: Font size for y-axis tick labels
      ),
      zaxis = list(
        title = "Frequency",
        titlefont = list(size = 18),  #JA: Font size for z-axis title
        tickfont = list(size = 18)  #JA: Font size for z-axis tick labels
      )
    ),
    legend = list(
      font = list(size = 20)  #JA: Font size for legend labels
    ),
    titlefont = list(size = 15)  #JA: Font size for the title
  )
  
  fig2
  htmlwidgets::saveWidget(fig2, paste("3D_surface_plot_fig5b_no_negative",i ,".html", sep = ""))
  
  
}

