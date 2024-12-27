# ------------------------------------------------------------------------------
#This code is part of paper: RNA:DNA hybrids determine tissue pathological heterogeneity and aggressiveness in clinically diagnosed cancers. 
#Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
#For updates and contributions, visit : https://github.com/adalaj

#Purpose: This script analyzes the distribution of signal intensities (DAB: Nucleus: Mean) for two experimental conditions—S9.6 and 
#RNase III+H—across various tumors (Normal, T1, T2, T3, and T4) in ICNST
#The outputs include density plots, statistical summaries, and files for further analysis.
# ------------------------------------------------------------------------------


#JA: load required libraries
library(data.table)
library(tidyverse)

#JA: Step1: Load dataset and preprocessing

rnase <- fread("BR720_IDC_only_RNase_treated_osk_No_negative.csv", sep = ",", header = TRUE)
#JA:169708

unique(rnase$TNM)
#JA: [1] "T1N0M0" "T1N1M0" "T2N1M0" "T2N0M0" "T4N0M0" "T3N0M0" "T4N2M0"
#JA:[8] "T2N2M0" "T4N1M0" "T3M0N0"

#JA:select desired columns
rnase<- rnase %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                         "DAB: Nucleus: Mean" ,"Pathology diagnosis", "TNM")

#JA:removing all negative values from fig8C input
rnase<- rnase %>% filter(`DAB: Nucleus: Mean` > 0)#JA: there is no negative values in this dataset



#JA:filtering for stage wise
tumor1 <- rnase %>% filter(TNM %in% c("T1N0M0", "T1N1M0"))
#JA:4604

tumor2<- rnase %>% filter(TNM %in% c("T2N1M0", "T2N0M0","T2N2M0"))
#JA:102174

tumor3<- rnase %>% filter(TNM %in% c("T3N0M0", "T3M0N0")) 
#JA:29486

tumor4<- rnase %>% filter(TNM %in% c("T4N0M0", "T4N2M0", "T4N1M0")) 
#JA:33444


#JA:inset normal file for RNase treated:9,473
normal<- fread("BR720_RNase treated_N_No Negative.csv", header = TRUE, sep = ",")
range(normal$`DAB: Nucleus: Mean`)#JA:JA: no negative values
#JA:[1] 0.0342 0.9654

normal <- normal %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                            "DAB: Nucleus: Mean" ,"Pathology diagnosis", "TNM")
#JA:260

#JA:Setting a column as identifier
tumor1$Identifier <- "T1"
tumor2$Identifier <- "T2"
tumor3$Identifier <- "T3"
tumor4$Identifier <- "T4"
normal$Identifier <- "Normal"

#JA: Step2: Data Transformation
#JA:Combing all three dataset as one for analyzing them together
fig8c_input<- rbind(tumor1, tumor2, tumor3, tumor4, normal)
#JA:4604+102174+29486+33444+260 = 169968

#JA:First find the range for binning
range<- range(fig8c_input$`DAB: Nucleus: Mean`)
#JA:[1] 0.0248  1.1394

a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA:add bin range in a new column

fig8c_input <- fig8c_input %>%
  dplyr::mutate(Bin = cut(`DAB: Nucleus: Mean`, breaks = bin_breaks_75, include.lowest = TRUE))
fig8c_input$Identifier <- factor(fig8c_input$Identifier, levels = c("Normal", "T1", "T2", "T3", "T4"))

fwrite(fig8c_input, "BR720_fig8c_tumor1_vs_2_vs_3_vs_4_vs_normal_RNase_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig8c_graph_input <- fig8c_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`DAB: Nucleus: Mean`, na.rm = TRUE),2),#JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA:Displaying bin midpoints
fig8c_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig8c_graph_input$Bin)
fig8c_graph_input <- fig8c_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig8c_graph_input$lower_limit<- as.numeric(fig8c_graph_input$lower_limit)
fig8c_graph_input$upper_limit<- as.numeric(fig8c_graph_input$upper_limit)
fig8c_graph_input<- fig8c_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA:rearrange the columns in graoh input
fig8c_graph_input<- fig8c_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#JA:for that find how many cores we have in each identifier set

normal_core<- length(unique(normal$`TMA core`))#JA:1
tumor1_core <- length(unique(tumor1$`TMA core`))#JA:4
tumor2_core<- length(unique(tumor2$`TMA core`))#JA:30
tumor3_core<- length(unique(tumor3$`TMA core`))#JA:9
tumor4_core <- length(unique(tumor4$`TMA core`))#JA:13


#JA: Step3: Normalization:
#JA: calculating average cells count compared between core wise in each bin
fig8c_graph_input <- fig8c_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Normal" ~ Frequency / normal_core,
    Identifier == "T1" ~ Frequency / tumor1_core,
    Identifier == "T2" ~ Frequency / tumor2_core,
    Identifier == "T3" ~ Frequency / tumor3_core,
    Identifier == "T4" ~ Frequency / tumor4_core,
    TRUE ~ NA_real_#JA:JA: Handles cases where Identifier do not match to either "Adjacent",  "tumor2", "tumor3"  or "tumor4"
  ))

#JA:#JA:calculate normalized average frequency                                                                                               
fig8c_graph_input<- fig8c_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Normal" ~ average_frequency / sum(fig8c_graph_input$average_frequency[fig8c_graph_input$Identifier == "Normal"]),
    Identifier == "T1" ~ average_frequency / sum(fig8c_graph_input$average_frequency[fig8c_graph_input$Identifier == "T1"]),
    Identifier == "T2" ~ average_frequency / sum(fig8c_graph_input$average_frequency[fig8c_graph_input$Identifier == "T2"]),
    Identifier == "T3" ~ average_frequency / sum(fig8c_graph_input$average_frequency[fig8c_graph_input$Identifier == "T3"]),
    Identifier == "T4" ~ average_frequency / sum(fig8c_graph_input$average_frequency[fig8c_graph_input$Identifier == "T4"]),
    TRUE ~ NA_real_#JA:JA: Handles cases where Identifier do not match to either "Adjacent",  "tumor2", or "tumor3"
  ))

fig8c_graph_input$Identifier <- factor(fig8c_graph_input$Identifier, levels = c("Normal", "T1", "T2", "T3", "T4"))
fwrite(fig8c_graph_input, "BR720_fig8c_tumor1_vs_2_vs_3_vs_4_vs_normal_RNase_geom_point.csv")

#JA: Step4: Visualization
max_value <- round(max(fig8c_graph_input$Bin_Midpoint), 1)



fig8c_graph<- ggplot()+
  geom_density(
 #JA: density plot for "S9.6 in normal"
    data = fig8c_input %>% filter(Identifier == "Normal"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_graph_input$norm_average_frequency[fig8c_graph_input$Identifier == "Normal"]) * diff(bin_breaks_75)[1],
      color = "Normal"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for T1
  
  geom_density(
    data = fig8c_input %>% filter(Identifier == "T1"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_graph_input$norm_average_frequency[fig8c_graph_input$Identifier == "T1"]) * diff(bin_breaks_75)[1],
      color = "T1"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for tumor2
  geom_density(
    data = fig8c_input %>% filter(Identifier == "T2"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_graph_input$norm_average_frequency[fig8c_graph_input$Identifier == "T2"]) * diff(bin_breaks_75)[1],
      color = "T2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for tumor3
  geom_density(
    data = fig8c_input %>% filter(Identifier == "T3"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_graph_input$norm_average_frequency[fig8c_graph_input$Identifier == "T3"]) * diff(bin_breaks_75)[1],
      color = "T3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for tumor4
  geom_density(
  data = fig8c_input %>% filter(Identifier == "T4"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_graph_input$norm_average_frequency[fig8c_graph_input$Identifier == "T4"]) * diff(bin_breaks_75)[1],
      color = "T4"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  geom_point(
    data = fig8c_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Normal" = "blue", "T1" = "darkolivegreen", "T2" = "red", "T3" = "black", "T4" = "darkorange")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different tumor size in ICNST",
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
        panel.background = element_blank(),#JA:facet wrap is not useful because of lot of tma_cores
        legend.position = "top")

#JA:change bw to 0.02
ggsave("BR720_fig8c_ICNST_tumor1_vs_2_vs_3_vs_4_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff", 
       plot = fig8c_graph, width=15,height=10, dpi=600)

#JA:change bw to 0.05
ggsave("BR720_fig8c_ICNST_tumor1_vs_2_vs_3_vs_4_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff", 
       plot = fig8c_graph, width=15,height=10, dpi=600)


#JA:Step5: Statistical analysis:
anova_result <- aov(`DAB: Nucleus: Mean` ~ Identifier, data = fig8c_input)
summary (anova_result)
#JA:                Df    Sum Sq   Mean Sq F value    Pr(>F)    
#JA:Identifier       4    182      45.38    1636      <2e-16 ***
#JA:  Residuals   169963   4714    0.03                   
#JA:---
#JA:Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_result <- TukeyHSD(anova_result)#JA:This output is from the Tukey HSD (Honest Significant Difference) test, which shows the results of pairwise comparisons between the means of the groups in your ANOVA model. Here's a breakdown of each part:

print(tukey_result)
#JA:Tukey multiple comparisons of means
#JA:95% family-wise confidence level

#JA:Fit: aov(formula = `DAB: Nucleus: Mean` ~ Identifier, data = fig8c_input)

#JA:$Identifier
#JA:Identifier
#JA:diff          lwr          upr     p adj
#JA:T1-Normal  0.044583611  0.015624187  0.073543035 0.0002589
#JA:T2-Normal  0.134231832  0.106021210  0.162442453 0.0000000
#JA:T3-Normal  0.127483886  0.099185143  0.155782629 0.0000000
#JA:T4-Normal  0.056902066  0.028617963  0.085186169 0.0000004
#JA:T2-T1      0.089648221  0.082803582  0.096492859 0.0000000
#JA:T3-T1      0.082900276  0.075701057  0.090099494 0.0000000
#JA:T4-T1      0.012318455  0.005177002  0.019459909 0.0000250
#JA:T3-T2     -0.006747945 -0.009751228 -0.003744662 0.0000000
#JA:T4-T2     -0.077329765 -0.080191813 -0.074467718 0.0000000
#JA:T4-T3     -0.070581820 -0.074211007 -0.066952634 0.0000000



#JA:for normal please use BR720_S9.6_N_No negative.csv
S9.6 <- fread("BR720_IDC_only_S9.6_positive_osk_No_negative.csv", sep = ",", header = TRUE)
#JA:269381


unique(S9.6$TNM)
#JA: [1] "T1N0M0" "T1N1M0" "T2N1M0" "T2N0M0" "T4N0M0" "T3N0M0" "T4N2M0"
#JA:[8] "T2N2M0" "T4N1M0" "T3M0N0"


S9.6<- S9.6 %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                       "DAB: Nucleus: Mean" ,"Pathology diagnosis", "TNM")

#JA:removing all negative values from fig8c_S9.6 input
S9.6<- S9.6 %>% filter(`DAB: Nucleus: Mean` > 0)#JA: there is no negative values in this dataset



#JA:filtering for grade wise
tumor1 <- S9.6 %>% filter(TNM %in% c("T1N0M0", "T1N1M0"))
#JA:9473

tumor2<- S9.6 %>% filter(TNM %in% c("T2N1M0", "T2N0M0","T2N2M0"))
#JA:167531

tumor3<- S9.6 %>% filter(TNM %in% c("T3N0M0", "T3M0N0")) 
#JA:47056

tumor4<- S9.6 %>% filter(TNM %in% c("T4N0M0", "T4N2M0", "T4N1M0")) 
#JA:45321


normal<- fread("BR720_S9.6_N_No negative.csv", header = TRUE, sep = ",")
range(normal$`DAB: Nucleus: Mean`)#JA: no negative values
#JA:[1]  0.0409 0.7994


normal <- normal %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                            "DAB: Nucleus: Mean" ,"Pathology diagnosis", "TNM")
#JA:789

#JA:Setting a column as identifier
tumor1$Identifier <- "T1"
tumor2$Identifier <- "T2"
tumor3$Identifier <- "T3"
tumor4$Identifier <- "T4"
normal$Identifier <- "Normal"


#JA:Combing all three dataset as one for analyzing them together
fig8c_S9.6_input<- rbind(tumor1, tumor2, tumor3, tumor4, normal)
#JA:9473+167531+47056+45321+789 = 270170

#JA:First find the range for binning
range<- range(fig8c_S9.6_input$`DAB: Nucleus: Mean`)
#JA:[1]  0.0119 1.3881

a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA:add bin range in a new column

fig8c_S9.6_input <- fig8c_S9.6_input %>%
  dplyr::mutate(Bin = cut(`DAB: Nucleus: Mean`, breaks = bin_breaks_75, include.lowest = TRUE))
fig8c_S9.6_input$Identifier <- factor(fig8c_S9.6_input$Identifier, levels = c("Normal", "T1", "T2", "T3", "T4"))

fwrite(fig8c_S9.6_input, "BR720_fig8c_tumor1_vs_2_vs_3_vs_4_vs_normal_S9.6_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig8c_S9.6_graph_input <- fig8c_S9.6_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`DAB: Nucleus: Mean`, na.rm = TRUE),2),#JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA:Displaying bin midpoints
fig8c_S9.6_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig8c_S9.6_graph_input$Bin)
fig8c_S9.6_graph_input <- fig8c_S9.6_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig8c_S9.6_graph_input$lower_limit<- as.numeric(fig8c_S9.6_graph_input$lower_limit)
fig8c_S9.6_graph_input$upper_limit<- as.numeric(fig8c_S9.6_graph_input$upper_limit)
fig8c_S9.6_graph_input<- fig8c_S9.6_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA:rearrange the columns in graoh input
fig8c_S9.6_graph_input<- fig8c_S9.6_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#JA:for that find how many cores we have in each identifier set

normal_core<- length(unique(normal$`TMA core`))#JA:1
tumor1_core <- length(unique(tumor1$`TMA core`))#JA:4
tumor2_core<- length(unique(tumor2$`TMA core`))#JA:30
tumor3_core<- length(unique(tumor3$`TMA core`))#JA:9
tumor4_core<- length(unique(tumor4$`TMA core`))#JA:13

#JA: calculating average cells count compared between core wise in each bin
fig8c_S9.6_graph_input <- fig8c_S9.6_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Normal" ~ Frequency / normal_core,
    Identifier == "T1" ~ Frequency / tumor1_core,
    Identifier == "T2" ~ Frequency / tumor2_core,
    Identifier == "T3" ~ Frequency / tumor3_core,
    Identifier == "T4" ~ Frequency / tumor4_core,
    TRUE ~ NA_real_#JA:JA: Handles cases where Identifier do not match to either "Adjacent",  "tumor2", "tumor3"  or "tumor4"
  ))

#JA:#JA:calculate normalized average frequency                                                                                               
fig8c_S9.6_graph_input<- fig8c_S9.6_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Normal" ~ average_frequency / sum(fig8c_S9.6_graph_input$average_frequency[fig8c_S9.6_graph_input$Identifier == "Normal"]),
    Identifier == "T1" ~ average_frequency / sum(fig8c_S9.6_graph_input$average_frequency[fig8c_S9.6_graph_input$Identifier == "T1"]),
    Identifier == "T2" ~ average_frequency / sum(fig8c_S9.6_graph_input$average_frequency[fig8c_S9.6_graph_input$Identifier == "T2"]),
    Identifier == "T3" ~ average_frequency / sum(fig8c_S9.6_graph_input$average_frequency[fig8c_S9.6_graph_input$Identifier == "T3"]),
    Identifier == "T4" ~ average_frequency / sum(fig8c_S9.6_graph_input$average_frequency[fig8c_S9.6_graph_input$Identifier == "T4"]),
    TRUE ~ NA_real_#JA:JA: Handles cases where Identifier do not match to either "Adjacent",  "tumor2", "tumor3" or "tumor4"
  ))

fig8c_S9.6_graph_input$Identifier <- factor(fig8c_S9.6_graph_input$Identifier, levels = c("Normal", "T1", "T2", "T3", "T4"))
fwrite(fig8c_S9.6_graph_input, "BR720_fig8c_tumor1_vs_2_vs_3_vs_4_vs_normal_S9.6_geom_point.csv")

max_value <- round(max(fig8c_S9.6_graph_input$Bin_Midpoint), 1)



fig8c_S9.6_graph<- ggplot()+
  geom_density(
 #JA: density plot for "S9.6 in normal"
    data = fig8c_S9.6_input %>% filter(Identifier == "Normal"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_S9.6_graph_input$norm_average_frequency[fig8c_S9.6_graph_input$Identifier == "Normal"]) * diff(bin_breaks_75)[1],
      color = "Normal"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for T1
  geom_density(
    data = fig8c_S9.6_input %>% filter(Identifier == "T1"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_S9.6_graph_input$norm_average_frequency[fig8c_S9.6_graph_input$Identifier == "T1"]) * diff(bin_breaks_75)[1],
      color = "T1"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for tumor2
  geom_density(
    data = fig8c_S9.6_input %>% filter(Identifier == "T2"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_S9.6_graph_input$norm_average_frequency[fig8c_S9.6_graph_input$Identifier == "T2"]) * diff(bin_breaks_75)[1],
      color = "T2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for tumor3
  geom_density(
    data = fig8c_S9.6_input %>% filter(Identifier == "T3"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_S9.6_graph_input$norm_average_frequency[fig8c_S9.6_graph_input$Identifier == "T3"]) * diff(bin_breaks_75)[1],
      color = "T3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
#JA:JA:plot for tumor4
  geom_density(
    data = fig8c_S9.6_input %>% filter(Identifier == "T4"),
    aes(
   #JA: x is values before binning
   #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8c_S9.6_graph_input$norm_average_frequency[fig8c_S9.6_graph_input$Identifier == "T4"]) * diff(bin_breaks_75)[1],
      color = "T4"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  geom_point(
    data = fig8c_S9.6_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Normal" = "blue", "T1" = "darkolivegreen", "T2" = "red", "T3" = "black", "T4" = "darkorange")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different tumor size in ICNST",
    subtitle ="S9.6 Only",
    x = "DAB intensity mean per nucleus",
    y = "Normalized Frequency",
    color = "Legends"
  ) +
  
#JA:JA:scale_x_log10()+
#JA:JA:scale_y_log10()+
  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),#JA:facet wrap is not useful because of lot of tma_cores
        legend.position = "top")

#JA:change bw to 0.05
ggsave("BR720_fig8c_ICNST_tumor1_vs_2_vs_3_vs_4_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff",
       plot = fig8c_S9.6_graph, width=15,height=10, dpi=600)

#JA:change bw to 0.02
ggsave("BR720_fig8c_ICNST_tumor1_vs_2_vs_3_vs_4_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff",
       plot = fig8c_S9.6_graph, width=15,height=10, dpi=600)


#JA:Statistics:
anova_result <- aov(`DAB: Nucleus: Mean` ~ Identifier, data = fig8c_S9.6_input)
summary (anova_result)
#JA:Df Sum Sq Mean Sq F value Pr(>F)    
#JA:Identifier       3    265   88.34    2545 <2e-16 ***
#JA: Residuals   270166   9377    0.03                   
#JA:---
#JA:Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_result <- TukeyHSD(anova_result)#JA:This output is from the Tukey HSD (Honest Significant Difference) test, which shows the results of pairwise comparisons between the means of the groups in your ANOVA model. Here's a breakdown of each part:

print(tukey_result)
#JA:Tukey multiple comparisons of means
#JA:95% family-wise confidence level

#JA:Fit: aov(formula = `DAB: Nucleus: Mean` ~ Identifier, data = fig8c_S9.6_input)

#JA:$Identifier
#JA:diff          lwr          upr     p adj
#JA:T1-Normal  -0.008914034 -0.027121824  0.009293755 0.5899692
#JA:T2-Normal   0.082013532  0.064942238  0.099084825 0.0000000
#JA:T3-Normal   0.008149211 -0.009017416  0.025315837 0.6143753
#JA:T2-T1  0.090927566  0.084426455  0.097428677 0.0000000
#JA:T3-T1  0.017063245  0.010315768  0.023810722 0.0000000
#JA:T3-T2 -0.073864321 -0.076195398 -0.071533244 0.0000000




