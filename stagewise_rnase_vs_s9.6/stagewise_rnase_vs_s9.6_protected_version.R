# ------------------------------------------------------------------------------
# This code is part of paper: RNA:DNA hybrids determine tissue pathological heterogeneity and aggressiveness in clinically diagnosed cancers. 
# Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
# For updates and contributions, visit : https://github.com/adalaj

# Purpose: This script analyzes the distribution of signal intensities (DAB: Nucleus: Mean) for two experimental conditions—S9.6 and 
#RNase III+H—across various stages (Normal, Stage 1, Stage 2, and Stage 3) in ICNST
#The outputs include density plots, statistical summaries, and files for further analysis.
# ------------------------------------------------------------------------------


#JA: load required libraries
library(data.table)
library(tidyverse)

#JA: Step1: Load dataset and preprocessing

rnase <- fread("BR720_IDC_only_RNase_treated_osk_No_negative.csv", sep = ",", header = TRUE)
#JA:169708

#JA:select desired columns
unique(rnase$Stage)
#JA:[1] "I"    "IIA"  "IIB"  "IIIB" "IIIA"

#JA:select desired columns
rnase<- rnase %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                         "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Stage")

#JA:removing all negative values from fig8b input
rnase<- rnase %>% filter(`DAB: Nucleus: Mean` > 0) #JA: there is no negative values in this dataset



#JA:filtering for stage wise
stage1 <- rnase %>% filter(Stage == "I")
#JA:2383

stage2<- rnase %>% filter(Stage %in% c("IIA", "IIB"))
#JA:129302

stage3<- rnase %>% filter(Stage %in% c("IIIA", "IIIB")) 
#JA:38023


#JA:inset normal file for RNase treated:
normal<- fread("BR720_RNase treated_N_No Negative.csv", header = TRUE, sep = ",")
range(normal$`DAB: Nucleus: Mean`) #JA: no negative values
#JA:[1] 0.0342 0.9654

normal <- normal %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                            "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Stage")
#JA:260

#JA:Setting a column as identifier
stage1$Identifier <- "Stage 1"
stage2$Identifier <- "Stage 2"
stage3$Identifier <- "Stage 3"
normal$Identifier <- "Normal"

#JA: Step2: Data Transformation
#JA:Combing all three dataset as one for analyzing them together
fig8b_input<- rbind(stage1, stage2, stage3, normal)
#JA:2383+129302+38023+260 = 169968

#JA:First find the range for binning
range<- range(fig8b_input$`DAB: Nucleus: Mean`)
#JA:[1] 0.0248  1.1394

a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA:add bin range in a new column

fig8b_input <- fig8b_input %>%
  dplyr::mutate(Bin = cut(`DAB: Nucleus: Mean`, breaks = bin_breaks_75, include.lowest = TRUE))
fig8b_input$Identifier <- factor(fig8b_input$Identifier, levels = c("Normal", "Stage 1", "Stage 2", "Stage 3"))

fwrite(fig8b_input, "BR720_fig8b_stage1_vs_2_vs_3_vs_normal_RNase_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig8b_graph_input <- fig8b_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`DAB: Nucleus: Mean`, na.rm = TRUE),2), #JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA:Displaying bin midpoints
fig8b_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig8b_graph_input$Bin)
fig8b_graph_input <- fig8b_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig8b_graph_input$lower_limit<- as.numeric(fig8b_graph_input$lower_limit)
fig8b_graph_input$upper_limit<- as.numeric(fig8b_graph_input$upper_limit)
fig8b_graph_input<- fig8b_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA:rearrange the columns in graoh input
fig8b_graph_input<- fig8b_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#JA:for that find how many cores we have in each identifier set

normal_core<- length(unique(normal$`TMA core`)) #JA:1
stage1_core <- length(unique(stage1$`TMA core`)) #JA:2
stage2_core<- length(unique(stage2$`TMA core`)) #JA:39
stage3_core<- length(unique(stage3$`TMA core`)) #JA:15

#JA: Step3: Normalization:
#JA: calculating average cells count compared between core wise in each bin
fig8b_graph_input <- fig8b_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Normal" ~ Frequency / normal_core,
    Identifier == "Stage 1" ~ Frequency / stage1_core,
    Identifier == "Stage 2" ~ Frequency / stage2_core,
    Identifier == "Stage 3" ~ Frequency / stage3_core,
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "stage2", or "stage3"  
  ))

#JA:#JA:calculate normalized average frequency                                                                                               
fig8b_graph_input<- fig8b_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Normal" ~ average_frequency / sum(fig8b_graph_input$average_frequency[fig8b_graph_input$Identifier == "Normal"]),
    Identifier == "Stage 1" ~ average_frequency / sum(fig8b_graph_input$average_frequency[fig8b_graph_input$Identifier == "Stage 1"]),
    Identifier == "Stage 2" ~ average_frequency / sum(fig8b_graph_input$average_frequency[fig8b_graph_input$Identifier == "Stage 2"]),
    Identifier == "Stage 3" ~ average_frequency / sum(fig8b_graph_input$average_frequency[fig8b_graph_input$Identifier == "Stage 3"]),
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "stage2", or "stage3"
  ))

fig8b_graph_input$Identifier <- factor(fig8b_graph_input$Identifier, levels = c("Normal", "Stage 1", "Stage 2", "Stage 3"))
fwrite(fig8b_graph_input, "BR720_fig8b_stage1_vs_2_vs_3_vs_normal_RNase_geom_point.csv")

#JA: Step4: Visualization
max_value <- round(max(fig8b_graph_input$Bin_Midpoint), 1)



fig8b_graph<- ggplot()+
  geom_density(
    #JA: density plot for "S9.6 in normal"
    data = fig8b_input %>% filter(Identifier == "Normal"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_graph_input$norm_average_frequency[fig8b_graph_input$Identifier == "Normal"]) * diff(bin_breaks_75)[1],
      color = "Normal"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+ 
  
  #JA:plot for Stage 1
  
  geom_density(
    data = fig8b_input %>% filter(Identifier == "Stage 1"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_graph_input$norm_average_frequency[fig8b_graph_input$Identifier == "Stage 1"]) * diff(bin_breaks_75)[1],
      color = "Stage 1"
    ),
    size = 1,#JA: just increase thickness of the lines
    kernel = "gaussian",
    bw = 0.05)+ 
  
  #JA:plot for stage2
  geom_density(
    data = fig8b_input %>% filter(Identifier == "Stage 2"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_graph_input$norm_average_frequency[fig8b_graph_input$Identifier == "Stage 2"]) * diff(bin_breaks_75)[1],
      color = "Stage 2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+ 
  
  #JA:plot for stage3
  geom_density(
    data = fig8b_input %>% filter(Identifier == "Stage 3"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_graph_input$norm_average_frequency[fig8b_graph_input$Identifier == "Stage 3"]) * diff(bin_breaks_75)[1],
      color = "Stage 3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+ 
  
  
  geom_point(
    data = fig8b_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Normal" = "blue", "Stage 1" = "darkolivegreen", "Stage 2" = "red", "Stage 3" = "black")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different stages in ICNST",
    subtitle ="RNase III+ RNase H",
    x = "DAB intensity mean per nucleus",
    y = "Normalized Frequency",
    color = "Legends"
  ) +
  
  #JA:scale_x_log10()+
  #JA:scale_y_log10()+
  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), #JA:facet wrap is not useful because of lot of tma_cores
        legend.position = "top")

#JA:change to 0.05bw
ggsave("BR720_fig8b_ICNST_stage1_vs_2_vs_3_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff", 
       plot = fig8b_graph, width=15,height=10, dpi=600)

#JA:change to 0.02bw
ggsave("BR720_fig8b_ICNST_stage1_vs_2_vs_3_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff", 
       plot = fig8b_graph, width=15,height=10, dpi=600)




#JA:Step5: Statistical analysis:
anova_result <- aov(`DAB: Nucleus: Mean` ~ Identifier, data = fig8b_input)
summary (anova_result)
#JA:                Df Sum Sq Mean Sq F value Pr(>F)    
#JA:Identifier       3    205   68.44    2480 <2e-16 ***
#JA: Residuals   169964   4691    0.03                   
#JA:---
#JA:  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_result <- TukeyHSD(anova_result) #JA:This output is from the Tukey HSD (Honest Significant Difference) test, which shows the results of pairwise comparisons between the means of the groups in your ANOVA model. Here's a breakdown of each part:

print(tukey_result)
#JA:Tukey multiple comparisons of means
#JA:95% family-wise confidence level

#JA:Fit: aov(formula = `DAB: Nucleus: Mean` ~ Identifier, data = fig8b_input)

#JA:$Identifier
#JA:diff          lwr         upr     p adj
#JA:Stage 1-Normal   0.02453931 -0.003335312  0.05241392 0.1071309
#JA:Stage 2-Normal   0.13456197  0.108067297  0.16105664 0.0000000
#JA:Stage 3-Normal   0.05587880  0.029320388  0.08243722 0.0000004
#JA:Stage 2-Stage 1  0.11002266  0.101199741  0.11884559 0.0000000
#JA:Stage 3-Stage 1  0.03133950  0.022326967  0.04035203 0.0000000
#JA:Stage 3-Stage 2 -0.07868317 -0.081172962 -0.07619337 0.0000000




setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/Ojashpi/TMA_data/TMA_final_figures/fig5")
S9.6 <- fread("BR720_IDC_only_S9.6_positive_osk_No_negative.csv", sep = ",", header = TRUE)
#JA:269381


unique(S9.6$Stage)
#JA:[1] "I"    "IIA"  "IIB"  "IIIB" "IIIA"


S9.6<- S9.6 %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                       "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Stage")

#JA:removing all negative values from fig8b_S9.6 input
S9.6<- S9.6 %>% filter(`DAB: Nucleus: Mean` > 0) #JA: there is no negative values in this dataset



#JA:filtering for grade wise
stage1 <- S9.6 %>% filter(Stage == "I")
#JA:5563

stage2<- S9.6 %>% filter(Stage %in% c("IIA", "IIB"))
#JA:211143

stage3<- S9.6 %>% filter(Stage %in% c("IIIA", "IIIB")) 
#JA:52675

normal<- fread("BR720_S9.6_N_No negative.csv", header = TRUE, sep = ",")
range(normal$`DAB: Nucleus: Mean`) #JA: no negative values
#JA:[1] 0.0409 0.7994

normal <- normal %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                            "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Stage")
#JA:789

#JA:Setting a column as identifier
stage1$Identifier <- "Stage 1"
stage2$Identifier <- "Stage 2"
stage3$Identifier <- "Stage 3"
normal$Identifier <- "Normal"


#JA:Combing all three dataset as one for analyzing them together
fig8b_S9.6_input<- rbind(stage1, stage2, stage3, normal)
#JA:5563+21143+52675+789 = 270170

#JA:First find the range for binning
range<- range(fig8b_S9.6_input$`DAB: Nucleus: Mean`)
#JA:[1]  0.0119 1.3881

a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA:add bin range in a new column

fig8b_S9.6_input <- fig8b_S9.6_input %>%
  dplyr::mutate(Bin = cut(`DAB: Nucleus: Mean`, breaks = bin_breaks_75, include.lowest = TRUE))
fig8b_S9.6_input$Identifier <- factor(fig8b_S9.6_input$Identifier, levels = c("Normal", "Stage 1", "Stage 2", "Stage 3"))

fwrite(fig8b_S9.6_input, "BR720_fig8b_stage1_vs_2_vs_3_vs_normal_S9.6_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig8b_S9.6_graph_input <- fig8b_S9.6_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`DAB: Nucleus: Mean`, na.rm = TRUE),2), #JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA:Displaying bin midpoints
fig8b_S9.6_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig8b_S9.6_graph_input$Bin)
fig8b_S9.6_graph_input <- fig8b_S9.6_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig8b_S9.6_graph_input$lower_limit<- as.numeric(fig8b_S9.6_graph_input$lower_limit)
fig8b_S9.6_graph_input$upper_limit<- as.numeric(fig8b_S9.6_graph_input$upper_limit)
fig8b_S9.6_graph_input<- fig8b_S9.6_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA:rearrange the columns in graoh input
fig8b_S9.6_graph_input<- fig8b_S9.6_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#JA:for that find how many cores we have in each identifier set

normal_core<- length(unique(normal$`TMA core`)) #JA:1
stage1_core <- length(unique(stage1$`TMA core`)) #JA:2
stage2_core<- length(unique(stage2$`TMA core`)) #JA:39
stage3_core<- length(unique(stage3$`TMA core`)) #JA:15

#JA: calculating average cells count compared between core wise in each bin
fig8b_S9.6_graph_input <- fig8b_S9.6_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Normal" ~ Frequency / normal_core,
    Identifier == "Stage 1" ~ Frequency / stage1_core,
    Identifier == "Stage 2" ~ Frequency / stage2_core,
    Identifier == "Stage 3" ~ Frequency / stage3_core,
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "stage2", or "stage3"  
  ))

#JA:#JA:calculate normalized average frequency                                                                                               
fig8b_S9.6_graph_input<- fig8b_S9.6_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Normal" ~ average_frequency / sum(fig8b_S9.6_graph_input$average_frequency[fig8b_S9.6_graph_input$Identifier == "Normal"]),
    Identifier == "Stage 1" ~ average_frequency / sum(fig8b_S9.6_graph_input$average_frequency[fig8b_S9.6_graph_input$Identifier == "Stage 1"]),
    Identifier == "Stage 2" ~ average_frequency / sum(fig8b_S9.6_graph_input$average_frequency[fig8b_S9.6_graph_input$Identifier == "Stage 2"]),
    Identifier == "Stage 3" ~ average_frequency / sum(fig8b_S9.6_graph_input$average_frequency[fig8b_S9.6_graph_input$Identifier == "Stage 3"]),
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "stage2", or "stage3"
  ))

fig8b_S9.6_graph_input$Identifier <- factor(fig8b_S9.6_graph_input$Identifier, levels = c("Normal", "Stage 1", "Stage 2", "Stage 3"))
fwrite(fig8b_S9.6_graph_input, "BR720_fig8b_stage1_vs_2_vs_3_vs_normal_S9.6_geom_point.csv")

max_value <- round(max(fig8b_S9.6_graph_input$Bin_Midpoint), 1)



fig8b_S9.6_graph<- ggplot()+
  geom_density(
    #JA: density plot for "S9.6 in normal"
    data = fig8b_S9.6_input %>% filter(Identifier == "Normal"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_S9.6_graph_input$norm_average_frequency[fig8b_S9.6_graph_input$Identifier == "Normal"]) * diff(bin_breaks_75)[1],
      color = "Normal"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+ #JA: changed to 0.01, 0.02, 0.03 and saved as BR720_fig8b_ICNST_stage123_normal_0.01bw.tiff
   
  
  #JA:plot for Stage 1
  
  geom_density(
    data = fig8b_S9.6_input %>% filter(Identifier == "Stage 1"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_S9.6_graph_input$norm_average_frequency[fig8b_S9.6_graph_input$Identifier == "Stage 1"]) * diff(bin_breaks_75)[1],
      color = "Stage 1"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for stage2
  geom_density(
    data = fig8b_S9.6_input %>% filter(Identifier == "Stage 2"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_S9.6_graph_input$norm_average_frequency[fig8b_S9.6_graph_input$Identifier == "Stage 2"]) * diff(bin_breaks_75)[1],
      color = "Stage 2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for stage3
  geom_density(
    data = fig8b_S9.6_input %>% filter(Identifier == "Stage 3"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8b_S9.6_graph_input$norm_average_frequency[fig8b_S9.6_graph_input$Identifier == "Stage 3"]) * diff(bin_breaks_75)[1],
      color = "Stage 3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  
  geom_point(
    data = fig8b_S9.6_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Normal" = "blue", "Stage 1" = "darkolivegreen", "Stage 2" = "red", "Stage 3" = "black")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different stages in ICNST",
    subtitle ="S9.6 Only",
    x = "DAB intensity mean per nucleus",
    y = "Normalized Frequency",
    color = "Legends"
  ) +
  
  #JA:scale_x_log10()+
  #JA:scale_y_log10()+
  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), #JA:facet wrap is not useful because of lot of tma_cores
        legend.position = "top")


ggsave("BR720_fig8b_ICNST_stage1_vs_2_vs_3_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff",
       plot = fig8b_S9.6_graph, width=15,height=10, dpi=600)
 
#JA:change bw 0.02
ggsave("BR720_fig8b_ICNST_stage1_vs_2_vs_3_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff",
       plot = fig8b_S9.6_graph, width=15,height=10, dpi=600)

#JA:Statistics:
anova_result <- aov(`DAB: Nucleus: Mean` ~ Identifier, data = fig8b_S9.6_input)
summary (anova_result)
#JA:Df Sum Sq Mean Sq F value Pr(>F)    
#JA:Identifier       3    265   88.34    2545 <2e-16 ***
#JA: Residuals   270166   9377    0.03                   
#JA:---
#JA:Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_result <- TukeyHSD(anova_result) #JA:This output is from the Tukey HSD (Honest Significant Difference) test, which shows the results of pairwise comparisons between the means of the groups in your ANOVA model. Here's a breakdown of each part:

print(tukey_result)
#JA:Tukey multiple comparisons of means
#JA:95% family-wise confidence level

#JA:Fit: aov(formula = `DAB: Nucleus: Mean` ~ Identifier, data = fig8b_S9.6_input)

#JA:$Identifier
#JA:diff          lwr          upr     p adj
#JA:Stage 1-Normal  -0.008914034 -0.027121824  0.009293755 0.5899692
#JA:Stage 2-Normal   0.082013532  0.064942238  0.099084825 0.0000000
#JA:Stage 3-Normal   0.008149211 -0.009017416  0.025315837 0.6143753
#JA:Stage 2-Stage 1  0.090927566  0.084426455  0.097428677 0.0000000
#JA:Stage 3-Stage 1  0.017063245  0.010315768  0.023810722 0.0000000
#JA:Stage 3-Stage 2 -0.073864321 -0.076195398 -0.071533244 0.0000000



